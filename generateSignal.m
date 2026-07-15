%% 生成混频信号 -- FMCW（调频连续波）雷达信号仿真
% 根据预设的雷达参数（如带宽、斜率、天线布局等）和目标信息，生成混频后的中频（IF）信号。

function [rawData, SNR_scene_total_dB] = generateSignal(Parameter, targetnum)
    kB = Parameter.kB;       % 玻尔兹曼常数
    T = Parameter.T;             % 环境温度 (K)
    Pt = Parameter.Pt;         % 发射功率 12 dBm
    G_ant = Parameter.G_ant;  % 天线增益 15 dBi
    Loss = Parameter.Loss;    % 系统损耗 4 dB
    NoiseFigure = Parameter.NoiseFigure; % 接收机噪声系数 10 dB
    if targetnum == 1 % 点目标
        RCS = 1; 
    else % 飞机目标
        RCS = 2; % 假设飞机单点目标的 RCS = 2 平方米
    end
    
    c = Parameter.c;                 %光速
    stratFreq = Parameter.stratFreq; %起始频率
    
    Tr = Parameter.Tr;            % 单个chirp内的扫频时间
    Tc = Parameter.Tc;            % 相邻 Chirp 的起始间隔时间
    samples = Parameter.Samples;  %采样点
    Fs = Parameter.Fs;            %采样率

    rangeBin = Parameter.rangeBin;     %rangeBin
    chirps = Parameter.Chirps;         %chirp数
    dopplerBin = Parameter.dopplerBin; %dopplerBin

    slope = Parameter.Slope;           %chirp斜率
    bandwidth = Parameter.Bandwidth;   %发射信号带宽
    centerFreq = Parameter.centerFreq; %中心频率
    lambda = Parameter.lambda;
    txAntenna = Parameter.txAntenna; %发射天线
    txNum = length(txAntenna);       %发射天线数
    rxAntenna = Parameter.rxAntenna; %接收天线
    rxNum = length(rxAntenna);       %接收天线数
    dx = Parameter.dx;          %水平间距
    
    target = Parameter.target;  %目标
    targetNum = size(target,1); %目标数
    rawData = zeros(txNum*rxNum,rangeBin,dopplerBin); % 雷达立方体 
    
    % 计算中频接收机热噪声功率 (时域总噪声功率)
    Pn = kB * T * Fs * NoiseFigure; % 大自然底噪 * 雷达网眼 * 芯片低噪声放大器

    t = 0:1/Fs:Tr-(1/Fs); % 单个chirp采样时间序列，采样点数Samples = 扫频时间*采样率 
    noise_sigma = sqrt(Pn / 2); % 复噪声两路分流（保持在循环外，避免重复计算开方）
    for chirpId = 1:chirps % 第i个Chirp，慢时间
       for txId = 1:txNum % 发射天线，模拟MIMO不同发射源
            % 即使 chirp 之间有间隔，载波的中心频率 centerFreq 依然在“滴答”地累积相位。
            % 如果不加上 (chirpId-1)*Tr，那么每一个新的 chirp 都会从相同的初始相位开始，这在物理上是不准确的（除非雷达在每个脉冲之间都关闭并重启振荡器）。
            % 发射信号的全局绝对起始时刻由 Tc 决定，即使在 Idle_time 雷达不发波，高频本振源的绝对相位依然随时间 Tc 轴连续累积
           
            % 【目前的车载雷达大多采用 TDM-MIMO（时分复用） 模式，
            % TDM 模式：TX1 发一个 Chirp，TX2 再发一个，然后 TX3……在这种情况下，天线之间存在时间偏移。
            % 第 1 轮 (chirpId=1): TX1 发射 -> TX2 发射 -> TX3 发射
            % 第 2 轮 (chirpId=2): TX1 发射 -> TX2 发射 -> TX3 发射
            % ...
            % 此时，txId 引起的位移会导致所谓的"多普勒相位偏移”，在处理角度（AoA）时需要补偿。】

            % (chirpId - 1) * txNum * Tc : 慢时间（不同脉冲间的时间偏移）
            % (txId - 1) * Tc : TDM-MIMO 下发射天线切换的时间偏移
            % 注意：如果是真正的 TDM，通常是一个时隙内只开一个 TX，
            % 实际偏移时间应该是：((chirpId - 1) * txNum + (txId - 1)) * Tc
            
            startTime_tx = ((chirpId - 1) * txNum + (txId - 1)) * Tc; % 【使用的Tc！！！】
            St = exp((1i*2*pi)*(centerFreq*(t+startTime_tx)+slope/2*t.^2)); % exp((1i*2*pi)*(centerFreq*(t+(chirpId-1)*Tr)+slope/2*t.^2)); % 确保每个天线发射时都带有正确的绝对时间戳
            
            for rxId = 1:rxNum % 接收天线，每个发射信号被所有接收天线接收
                Sif_sum = zeros(1,rangeBin); 
                for targetId = 1:targetNum % 遍历所有预设目标，叠加它们的回波

                    targetRange = target(targetId,1);
                    targetSpeed = target(targetId,2); 
                    targetAngle = target(targetId,3);
                    
                    % 计算目标到达雷达的【理论接收功率】 (瓦特)
                    Pr = (Pt * G_ant * G_ant * lambda^2 * RCS) / ((4*pi)^3 * targetRange^4 * Loss);
                    A_rx = 2 * sqrt(Pr); % 反推电压振幅：中频信号振幅 = 2 * sqrt(Pr)， *2为了把混频分一半的能量补偿回来
                     
                    % 加上速度乘以时间，【目标远离雷达速度为正，靠近雷达速度为负】
                    tau = 2 * targetRange / c; % 实际上 tau 在一个 Chirp 内部也有微小变化，但通常忽略不计（Stop-and-Go 假设）
                    fd = -2 * targetSpeed / lambda; % 多普勒频移：靠近为负速度，算出正的 fd
                    % 核心：空间相位：不同天线接收同一个目标时，因为位置不同，信号到达时间有微小差异，体现为相位差 wx
                    % (txId-1) * rxNum + rxId)：当前虚拟通道
                    % dx * sind(targetAngle)：相邻天线路径差，再除以λ转成相位，（2π对应一个波长）
                    % (txId-1) * rxNum + rxId 假设了虚拟天线在物理空间上是排成一排的线性阵列
                    wx = ((txId-1) * rxNum + rxId-1) / lambda * dx * sind(targetAngle);
				    % 混频中频信号相位
				    phi_IF = 2 * pi * (slope * tau * t + centerFreq * tau - fd * startTime_tx + wx);% - wx
                    Sif_sum = Sif_sum + A_rx * exp(1i * phi_IF);
                end
                noise = noise_sigma * (randn(1, rangeBin) + 1i*randn(1, rangeBin));
                Sif_sum = Sif_sum + noise; 
                rawData((txId-1) * rxNum + rxId,:,chirpId) = Sif_sum;
            end
        end
    end
    
    
    % ================== 【SNR 计算】 ==================
    unique_target_ids = unique(target(:, 4)); 
    total_objects = length(unique_target_ids);
    
    % 初始化一个结构体数组，用来分别完整记录每一个大目标的 SNR
    SNR_scene_total_dB = struct('TargetID', [], 'Total_SNR_dB', [], 'Points_SNR_dB', []);
    
    % 2. 遍历每一个物理整体目标（Object）
    for objIdx = 1:total_objects
        current_obj_id = unique_target_ids(objIdx);
        
        % 找出属于当前目标整体的所有散射点（Rows）
        point_indices = find(target(:, 4) == current_obj_id);
        num_points = length(point_indices);
        
        obj_total_Pr = 0; % 用于累加该目标所有散射点的总接收功率
        points_snr_list = zeros(num_points, 1); % 记录当前大目标下每个散射点的 SNR
        
        % 3. 遍历该目标内部的每一个散射点（Point）
        for pIdx = 1:num_points
            row_idx = point_indices(pIdx);
            r_init = target(row_idx, 1);
            
            % 计算该散射点的单点理论接收功率
            Pr_point = (Pt * G_ant * G_ant * lambda^2 * RCS) / ((4*pi)^3 * r_init^4 * Loss);
            obj_total_Pr = obj_total_Pr + Pr_point; % 累加到整体目标功率账本中
            
            % 计算该散射点单兵作战时，经 2D-FFT 压缩后的独立检测 SNR
            SNR_point_fft = (Pr_point / Pn) * samples * chirps;
            points_snr_list(pIdx) = 10 * log10(SNR_point_fft);
            
        end
        
        % 4. 计算整个宏观目标聚合后的总信噪比
        SNR_object_total = (obj_total_Pr / Pn) * samples * chirps;
        SNR_object_total_dB = 10 * log10(SNR_object_total);
        
            
        % 5. 将当前目标的完整账本归档到结构体数组的对应格子里
        SNR_scene_total_dB(objIdx).TargetID = current_obj_id;
        SNR_scene_total_dB(objIdx).Total_SNR_dB = SNR_object_total_dB;
        SNR_scene_total_dB(objIdx).Points_SNR_dB = points_snr_list;
    end
    fprintf('==================================================================\n');
    % ====================================================================


end