function [pointCloud] = doa(doaMethod, radarData3D, pointList, cfarRD, Parameter)
    % radarData3D: [channel, Range, Doppler] 
    % doaMethod: 选择测角算法 【1:fft  2:dbf  3:1Dmusic  4:RDMusic  5:Capon】
    
    virtualAntenna = Parameter.virtualAntenna;
    lambda = Parameter.lambda;
    txNum = Parameter.txNum;
    rxNum = Parameter.rxNum;
    dopplerLen = Parameter.dopplerBin;
    dx = Parameter.dx;
    numPoints = size(pointList, 2);
%     pointCloud = zeros(numPoints, 4); % 4D点云
    pointCloud = []; % 【MUSIC或Capon算法可能在一个RD单元格子里同时解析出多个角度，所以可能会比numPoints点数多】
    for targetId = 1:numPoints
        rIdx = pointList(1, targetId);
        dIdx = pointList(2, targetId);
        rangeVal = (rIdx-1) * Parameter.c/2/Parameter.BandwidthValid;
        antVec = squeeze(radarData3D(:,rIdx,dIdx)); % 天线相位向量
        
        % TDM-MIMO 多普勒相位补偿
        binIdx = dIdx - (dopplerLen/2 + 1); % 例如：512 个 Bin，还原后的 binIdx 范围是 -256 到 +255。
        
        %【判断是否速度模糊折叠】
%         if Parameter.target(1,2) < 0 && binIdx > 0 % 仅适用于一倍折叠
%             % 检测到目标本质是靠近的（预设为负），但 FFT 错报成了正速度，说明发生了爆表折叠
%             % 比如目标本该是深浅负速度（-180 频点），但由于爆表折叠，在 FFT 谱上它荡到了正速度（+332 频点）。
%             % 让 binIdx 减去总通道数，将其拉回真实的物理频点。
%             binIdx = binIdx - dopplerLen; 
%         end
        
        velVal = binIdx * lambda / (2 * Parameter.Tc * txNum * Parameter.Chirps); % binIdx*vel_res
        
        % 【根据信号生成代码：currentTime = ((chirpId - 1) * txNum + (txId - 1)) * Tr;】
        % 【相邻两个 chirpId 之间的时间跨度是 txNum * Tr（因为中间要等所有 TX 轮流发完）。】
        % 【所以，相邻chirp相位演变要除以 txNum】
        dopplerPhase = 2 * pi * binIdx / dopplerLen / txNum; % 相邻 Chirp 间隔的相位演变
        % 补偿 antVec
        % 1:4 通道是 TX1 (基准，无需补偿)
        % 5:8 通道是 TX2 (延迟了 1 个 TDM 槽位) -> 补偿 -1 * dopplerPhase
        % 9:12 通道是 TX3 (延迟了 2 个 TDM 槽位) -> 补偿 -2 * dopplerPhase
        % 注意：具体延迟几个槽位取决于你 generateSignal 里的 txId 顺序
        for txId = 1:txNum
            compFactor = exp(-1j * dopplerPhase * (txId - 1)); 
            idx_start = (txId - 1) * rxNum + 1;
            idx_end = txId * rxNum;
            antVec(idx_start : idx_end) = antVec(idx_start : idx_end) .* compFactor;
        end
            
        if doaMethod == 1 % FFT 测角
            angleIndex = asin((-512:1:512-1)/512) * 180 / pi;
            doa_fft=fftshift(fft(antVec,1024));
            doa_abs=abs(doa_fft);
            [~,loc]=max(doa_abs);
            angle = angleIndex(loc);
            
        elseif doaMethod == 2 % DBF 数字波束形成
            thetaScan = -90:0.1:90; % 定义波束扫描范围和精度，从左侧 -90°扫到右侧 90°，步长为 0.1°。这决定了测角的显示分辨率。
            weightVec = zeros(virtualAntenna,1);
            doa_dbf = zeros(length(thetaScan),1);
            kk = 1;
            for degscan = thetaScan
                for txId = 1:txNum
                    for rxId = 1:rxNum
                        % 计算理论上的相位偏移
                        dphi = ((txId-1) * rxNum + rxId - 1) * 2 * pi / lambda * dx * sind(degscan); % 计算在当前假设扫描角度 degscan 下，第 N个虚拟通道应该具有的相位偏移。
                        % 构建导向矢量
                        weightVec((txId-1) * rxNum + rxId) = exp(1i * dphi); % 【为了抵消掉中频信号中已经存在的那个正相位偏移wx，见generateSignal.m】
                    end
                end
                doa_dbf(kk) = antVec'*weightVec; % 能量投影(内积)
                kk = kk + 1;
            end
            doa_abs = abs(doa_dbf);
            [pk,loc]=max(doa_abs); % 寻找能量最高的索引
            angle = thetaScan(loc); % 映射回角度值
    
        elseif doaMethod == 3 % 1DMUSIC
            % 普通的 MUSIC 算法能分离开不同的角度，有一个核心的大前提：
            % 不同目标发射或反射回来的信号，在统计上必须是独立的（不相关的）。
            % 空间平滑的本质是：利用“切香肠”的方式，牺牲一部分虚拟阵列的长度，来换取矩阵之秩的恢复。
            % --- Step 1: 空间平滑 ---
            subLen = 8; % 子阵长度，天线阵列孔径变小，导致角度分辨率变宽了一点
            numSub = virtualAntenna - subLen + 1; % 子阵数量 
            Rx_smooth = zeros(subLen, subLen);
            
            for k = 1:numSub
                subV = antVec(k : k + subLen - 1);
                Rx_smooth = Rx_smooth + (subV * subV'); % 计算子阵协方差矩阵[subLen, subLen]
            end
            Rxx = Rx_smooth / numSub; % 所有子阵协方差矩阵取平均
            
            % --- Step 2: 对角加载 (增强数值稳定性) ---
            % 加一个相对于矩阵迹(Trace)很小的量，解决矩阵接近奇异的问题
%             Rxx = Rxx + 2e-3 * eye(subLen);
            
            % --- Step 4: MUSIC 核心算法 ---
            thetaScan = -90:0.5:90;
            
            [V, D] = eig(Rxx);
            [~, sortIdx] = sort(diag(D), 'descend');
            En = V(:, sortIdx(2:end)); % [subLen, subLen-1] 噪声子空间，【现在是针对每个CFAR点逐个处理，视为该分辨单元内包含1个目标】
            
            P_music = zeros(length(thetaScan), 1);
            for j = 1:length(thetaScan)
                A = exp(1i * 2 * pi * dx / lambda * (0:subLen-1)' * sind(thetaScan(j)));% [subLen, 1]
                P_music(j) = 1 / abs(A' * (En * En') * A); % 计算加了括号，所以结果是标量
            end
             
            % --- Step 5: 提取峰值并转换坐标 ---
            [~, maxIdx] = max(P_music);
            angle = thetaScan(maxIdx);

        elseif doaMethod == 4 % RDMusic
            % 这里的RDMusic实际上还是一维Music，只是做了多目标多峰值提取
            % --- Step 1: 空间平滑 ---
            subLen = 8; % 子阵长度
            numSub = virtualAntenna - subLen + 1; % 子阵数量 
            Rx_smooth = zeros(subLen, subLen);
            
            for k = 1:numSub
                subV = antVec(k : k + subLen - 1);
                Rx_smooth = Rx_smooth + (subV * subV');
            end
            Rx = Rx_smooth / numSub; 
            % --- Step 2: 自适应对角加载 (增强数值稳定性) ---
            % 加一个相对于矩阵迹(Trace)很小的量，解决矩阵接近奇异的问题
            loading_factor = 0.01; % 常用经验值 1% ~ 5%
            Rx = Rx + loading_factor * trace(Rx) * eye(subLen); % [subLen, subLen]

            % --- Step 3:特征值分解
            [V, D] = eig(Rx);
            [~, sortIdx] = sort(diag(D), 'descend');
            V = V(:, sortIdx);

            % 通过特征值能量占比来判断目标数
            eigVals = diag(D);
            eigVals = eigVals(sortIdx);
            energyTh = sum(eigVals) * 0.1; 
            numTargets = sum(eigVals > energyTh); % 能量占比大于10%的视为有效目标信号
            numTargets = max(1, min(numTargets, 2)); % 强制限制单点最多解析出2个并排目标，最少不低于1个目标
            En = V(:, numTargets + 1 : end);% 噪声子空间 % [subLen, numNoise]
            
            deg = -90:0.5:90;
            subAntIdx = (0:subLen-1)';
            dphi_matrix = 2 * pi / lambda * dx * subAntIdx * sind(deg); % [subLen, deg]
            A_matrix = exp(1i * dphi_matrix);
            
            % --- MUSIC 核心数学公式 矩阵化并行计算 ---
            % P_music = 1 / (a' * En * En' * a)，令y = En' * a
            % y' * y = |y1|^2+|y2|^2+...+|yn|^2
            En_A = En' * A_matrix; % [numNoise, deg]
            denominator = sum(abs(En_A).^2, 1); % 将角度导向矢量在噪声维度的分量平方和累加
            music_spectrum = 1 ./ (denominator + eps); % 加上 eps 防止分母为 0
            music_spectrum_abs = abs(music_spectrum') ; % 转为列向量

            [~, locs] = findpeaks(music_spectrum_abs, ...
                            'MinPeakHeight', max(music_spectrum_abs) * 0.05, ...
                            'SortStr', 'descend');
            if ~isempty(locs)
                angle = deg(locs); % 吐出通过超分辨分离开的【所有】目标角度数组
            else
                angle = 0; % 兜底保护
            end
            
%             % 真正的 2DMusic【同时测出，方位角+俯仰角】 伪代码如下
%             % 导向矢量需要同时考虑 X 轴和 Y 轴的天线坐标
%             thetaScan = -90:1:90; % 方位角扫描范围
%             phiScan   = -30:1:30; % 俯仰角扫描范围
%             
%             P_music_2D = zeros(length(thetaScan), length(phiScan));
%             
%             % 进行高耗能的二维双重循环搜索（或者矩阵化扩展）
%             for t = 1:length(thetaScan)
%                 for p = 1:length(phiScan)
%                     % 构造同时包含方位 theta 和俯仰 phi 的 2D 导向矢量
%                     A = exp(1i * 2 * pi / lambda * (x_idx*sind(thetaScan(t))*cosd(phiScan(p)) + y_idx*sind(phiScan(p))));
%                     
%                     % 映射到 2D 空间谱图上
%                     P_music_2D(t, p) = 1 / abs(A' * (En * En') * A);
%                 end
%             end
%             % 最终在 P_music_2D 这张二维热力图上寻找局部极大值峰值

        elseif doaMethod == 5 % Capon
            % --- Step 1: 空间平滑 ---
            subLen = 8; % 子阵长度
            numSub = virtualAntenna - subLen + 1; % 子阵数量 
            Rx_smooth = zeros(subLen, subLen);
            
            for k = 1:numSub
                subV = antVec(k : k + subLen - 1);
                Rx_smooth = Rx_smooth + (subV * subV');
            end
            Rx = Rx_smooth / numSub; 
            
            % --- Step 2: 自适应对角加载 (增强数值稳定性) ---
            % 加一个相对于矩阵迹(Trace)很小的量，解决矩阵接近奇异的问题
            loading_factor = 0.01; % 常用经验值 1% ~ 5%
            Rx = Rx + loading_factor * trace(Rx) * eye(subLen);
        
            deg = -90:0.5:90;
            subAntIdx = (0:subLen-1)'; % 只需要子阵内部的相对索引
        %     a = zeros(1,virtualAntenna);
            kk = 1;
            doa_capon = zeros(length(deg),1);
            for degscan = deg
                dphi = 2 * pi / lambda * dx * subAntIdx * sind(degscan);
                a = exp(1i * dphi);
        
                % 数学等价于 a' * (Rx \ a)
                % Rx \ a 实际上解出了 Rx * w = a 中的 w
                doa_capon(kk) = 1/(a'*(Rx \ a));
                kk = kk + 1;
            end
            doa_abs = abs(doa_capon);
            % 【多目标下的找"峰"逻辑需调整，分辨在同一个距离和速度下，并排存在的多个目标"（比如两辆并排开的车，或并排走的两个人）】
            % 多目标下可使用findpeaks函数，寻找局部极大值峰值，[pks, locs] = findpeaks(doa_abs,'MinPeakHeight', max(doa_abs)*0.1, 'SortStr', 'descend');
            % locs 可能是一个数组，包含1个或多个并排目标的角度
            [~,loc]=max(doa_abs); 
            angle = deg(loc);
            
        end
        
        current_power = cfarRD(rIdx, dIdx); % 这个cfarRD已是db幅度
        pointCloud = [pointCloud; rangeVal, velVal, angle, current_power];

    end
end