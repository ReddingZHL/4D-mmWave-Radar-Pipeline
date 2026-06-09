% cfar_detection的改进
function [pointList, cfarRD] = cfar_detection_1D_pro(accumulateRD, targetnum)
    % 输入说明：accumulateRD 假定为线性功率谱（如果是dB谱，建议先 10.^(RD/10) 还原）
    [rangeLen, dopplerLen] = size(accumulateRD);
    
    % 1. 参数提取与配置
    if targetnum == 1 % 1-点目标
        d_guard  = 2;  d_train  = 8;  
        r_guard  = 4;  r_train  = 8;  
        d_SNR  = 10; % 多普勒维参数
        r_SNR  = 10; % 距离维参数
    else % 飞机参数
        %【仿真的飞机形状有将近一百个点，降低d_SNR和r_SNR，cfar输出的点从5个锐增到50个】
        d_guard  = 4;  d_train  = 8;  
        r_guard  = 5;  r_train  = 8;  
        d_SNR  = 5.5;3.5; % 多普勒维参数
        r_SNR  = 5.5;3.5; % 距离维参数
    end
    
    
    % 2. 使用 1D 卷积算子 【舍去多普勒维双重循环】
    % 构造多普勒方向的滑窗卷积核 -- 直接计算出噪声均值
    % 长度为 d_train + d_guard + 1 + d_guard + d_train
    total_d_len = 2 * (d_train + d_guard) + 1;
    k_doppler_left  = zeros(1, total_d_len);
    k_doppler_right = zeros(1, total_d_len);
    % 左侧训练区 
    k_doppler_left(1 : d_train) = 1 / d_train;
    % 右侧训练区
    k_doppler_right(end - d_train + 1 : end) = 1 / d_train;

    % 环形边界处理（对应环形填充，imfilter 内置 'circular'）
    leftNoise  = imfilter(accumulateRD, k_doppler_left,  'circular');
    rightNoise = imfilter(accumulateRD, k_doppler_right, 'circular');
    % 1:CA-CFAR  2:GO-CFAR    3:SO-CFAR
    noise_CA = (leftNoise + rightNoise) / 2;
    noise_GO = max(leftNoise, rightNoise);
    noise_SO = min(leftNoise, rightNoise);
    
    % 默认采用 CA-CFAR，直接矩阵相除得到全图所有点的信噪比（线性域用除法）
    snr_map_doppler = accumulateRD ./ noise_CA; % [rangeLen, dopplerLen] 
    
    % 找出通过多普勒检测的所有列（多普勒Bin）
    th_linear_d = 10^(d_SNR / 10); % 将 dB 门限转换为线性倍数
    doppler_mask = any(snr_map_doppler > th_linear_d, 1); % [1, dopplerLen] 沿着（参数1）矩阵的第1维（纵向/逐行）
    dopplerCfarList = find(doppler_mask);
    
    % 3. 距离维度的矩阵级级并行计算
    cfarRD = zeros(rangeLen, dopplerLen);
    pointList = [];
    
    if isempty(dopplerCfarList)
        return;
    end
    
    % 构造距离向的上滑窗和下滑窗卷积核 (纵向列向量)
    total_r_len = 2 * (r_train + r_guard) + 1;
    k_range_up   = zeros(total_r_len, 1);
    k_range_down = zeros(total_r_len, 1);
    % 真正的上方训练区
    k_range_up(1 : r_train) = 1 / r_train;
    % 真正的下方训练区
    k_range_down(end - r_train + 1 : end) = 1 / r_train;

    % 对称边界处理（对应对称填充，imfilter 内置 'symmetric'）
    upNoise   = imfilter(accumulateRD, k_range_up,   'symmetric');
    downNoise = imfilter(accumulateRD, k_range_down, 'symmetric');
    
    noise_range_CA = (upNoise + downNoise) / 2;
    snr_map_range = accumulateRD ./ noise_range_CA;
    th_linear_r = 10^(r_SNR / 10);
    
    % 4. 两维 Mask 交叉相乘，瞬间提取点迹
    % 只有同时通过多普勒列筛选，且满足距离维 SNR 的点才留下
    final_mask = (snr_map_range > th_linear_r); % [rangeLen, dopplerLen] 
    final_mask(:, ~doppler_mask) = 0; % 裁剪掉没通过多普勒筛选的列
    
    % 提取最终结果
    [r_hits, d_hits] = find(final_mask);
    pointList = [r_hits'; d_hits']; % 完美的 [2 x N] 坐标矩阵
    
    % 填充输出谱
    linear_indices = find(final_mask);
    cfarRD(linear_indices) = accumulateRD(linear_indices);
    
    % 如果整体流水线后面需要 dB，最后统一转回 dB
    cfarRD = 10 * log10(cfarRD + eps); 
end