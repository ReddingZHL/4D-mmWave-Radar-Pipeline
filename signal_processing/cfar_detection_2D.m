function [pointList, cfarRD] = cfar_detection_2D(accumulateRD, targetnum)
    % 输入说明：accumulateRD 传入的为非相干积分线性功率谱
    [rangeLen, dopplerLen] = size(accumulateRD);
    
    % =======================================================================
    % 1. 2D 物理空间滑窗参数自适应配置
    % =======================================================================
    if targetnum == 1 % 1-点目标参数
        r_guard = 2;  r_train = 4;  % 距离维保护、参考【半宽】
        d_guard = 2;  d_train = 4;  % 多普勒维保护、参考【半宽】
        cfar_SNR_dB = 12;           % 2D 联合恒虚警检测门限
    else % 飞机扩展目标参数
        %【针对飞机目标的调谐：2D 下适当放宽保护区，降低门限以捞出机翼等散射点迹】
        r_guard = 4;  r_train = 6; 8; 
        d_guard = 4;  d_train = 6;  8;
        cfar_SNR_dB = 5.5;      3.5;    
    end
    
    % =======================================================================
    % 2. 构建“二维空心双层方砖”卷积核 (Nested Matrix)
    % =======================================================================
    % 计算大方砖（外圈参考窗）的完整边长
    N_r_total = 2 * (r_train + r_guard) + 1;
    N_d_total = 2 * (d_train + d_guard) + 1;
    
    % 步骤 A：创建一个充满 1 的全尺寸实心大方砖
    kernel_2D = ones(N_r_total, N_d_total);
    
    % 步骤 B：计算出内圈“保护方砖”在矩阵内部的绝对坐标边界
    r_center = r_train + r_guard + 1;
    d_center = d_train + d_guard + 1;
    
    r_guard_range = (r_center - r_guard) : (r_center + r_guard);
    d_guard_range = (d_center - d_guard) : (d_center + d_guard);
    
    % 步骤 C：将中间保护区（含CUT本身）全部强行抹零，形成“方砖城墙”
    kernel_2D(r_guard_range, d_guard_range) = 0;
    
    % 步骤 D：计算这堵城墙到底由多少个纯噪声格子组成，进行均值归一化
    num_noise_cells = sum(kernel_2D(:)); 
    kernel_2D = kernel_2D / num_noise_cells; 
    
    % =======================================================================
    % 3. 2D 空间滑窗
    % =======================================================================
    % 边界处理：距离向（纵向）采用对称反射填充 'symmetric'
    
    pad_size_d = d_train + d_guard; % 拼接宽度，多普勒向（横向）由于速度具有周期折叠特性
    % 横向环形手动填充：把右边的数据贴到左边，左边的数据贴到右边
    padded_RD = [accumulateRD(:, end-pad_size_d+1:end), accumulateRD, accumulateRD(:, 1:pad_size_d)];
    padded_noise = imfilter(padded_RD, kernel_2D, 'symmetric'); % 对称反射填充
    % 裁剪回原有 2D 图尺寸
    mean_noise_floor = padded_noise(:, pad_size_d+1 : end-pad_size_d); % 还原原始 RD 图尺寸
    
    % =======================================================================
    % 4. 2D 联合信噪比对账与点迹提取
    % =======================================================================
    % 线性域直接相除，得到整张二维图每一个坐标点的真实 2D SNR Map
    snr_map_2D = accumulateRD ./ (mean_noise_floor + eps);
    
    % 转换门限
    th_linear = 10^(cfar_SNR_dB / 10);
    
    % 生成最终的 2D 检测掩膜 (1代表冲破门限的强点，0代表噪底)
    final_mask = (snr_map_2D > th_linear);
    
    % 剔除边缘极其靠边、无法形成完整保护区的死角点迹（可选防呆机制）
    final_mask(1:r_guard, :) = 0; final_mask(end-r_guard+1:end, :) = 0;
    
    % 提取最终通关的散射点坐标
    [r_hits, d_hits] = find(final_mask);
    
    % 完美的 [N x 2] 坐标矩阵（满足你后续 DOA 测角函数输入要求的 [rangeIdx, dopplerIdx]）
    pointList = [r_hits'; d_hits']; 
    
    % =======================================================================
    % 5. 谱图填充与格式转换（转回 dB 供后级显示或调试）
    % =======================================================================
    cfarRD = zeros(rangeLen, dopplerLen);
    linear_indices = find(final_mask);
    cfarRD(linear_indices) = accumulateRD(linear_indices);
    
    % 统一转回符合雷达可视化账本的 dB 格式
    cfarRD = 10 * log10(cfarRD + eps); 
end