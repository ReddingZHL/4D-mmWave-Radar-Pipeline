function [range_fft, rd_map, Metrics] = range_doppler_fft(rawData, Parameter)
    % rawData: [虚拟通道, 距离采样点, chirp数]
    numADCSamples = Parameter.Samples; % 256
    numChirps = Parameter.Chirps; % 512

    %% ================== 1. 理论与实际分辨率计算 ==================
    c = Parameter.c;
    bandwidth = Parameter.BandwidthValid;
    lambda = Parameter.lambda;
    txNum = Parameter.txNum;
    Tc = Parameter.Tc;
    
    % 理论分辨率
    Metrics.res_r_theory = c / (2 * bandwidth);
    Metrics.res_v_theory = lambda / (2 * txNum * Tc * numChirps);
    % 汉宁窗（Hann）的主瓣展宽因子约为 1.63（主瓣变宽，分辨率退化）
    win_broadening_factor_r = 1.63; 
    win_broadening_factor_v = 1.63;
    Metrics.res_r_actual = Metrics.res_r_theory * win_broadening_factor_r;
    Metrics.res_v_actual = Metrics.res_v_theory * win_broadening_factor_v;

    %% ================== 加窗处理链路 ==================
    % 加窗：压低距离维和速度维的旁瓣，避免强目标的泄露掩盖掉弱目标 
    win_r = hann(numADCSamples); % 定义汉宁窗[numADCSamples, 1]
    win_d = hann(numChirps); % [numChirps, 1]
    filter_2d = win_r * win_d'; % 二维滤波矩阵 [numADCSamples, numChirps]
    filter_3d = reshape(filter_2d, 1, numADCSamples, numChirps); % 三维滤波矩阵
    % 广播机制：如果两个矩阵在某个维度上的大小不一致，但其中一个是 1，那么 MATLAB 会自动将这个维度为 1 的矩阵沿着该维度进行“虚拟复制”，直到它的尺寸和另一个矩阵匹配。
    rawData_windowed = rawData .* filter_3d;
    
    % Range FFT
    range_fft = fft(rawData_windowed, numADCSamples, 2);
    % 静态杂波抑制 (可选：减去均值)
    % range_fft = range_fft - mean(range_fft, 3);
    
    % Doppler FFT
    rd_map = fftshift(fft(range_fft, numChirps, 3), 3);
    
%     % 取幅度（实际工程常取第一个或所有虚拟通道非相干平均）。先取模，再平均（防止目标因相位差被抵消）
%     rd_map_mean = squeeze(mean(abs(rd_map).^2, 1));% [numADCSamples, numChirps]
%     Metrics = calculate_blind_metrics(rd_map_mean, Parameter);
end

function Metrics = calculate_blind_metrics(rd_map_mean, Parameter)
    % rd_map_mean: 传入的2D-FFT非相干积分能量图 [rangeBin, dopplerBin]
    
    %% ================== 1. 盲测 Peak SNR ==================
    [max_val, max_idx] = max(rd_map_mean(:));
    [max_r, max_v] = ind2sub(size(rd_map_mean), max_idx);
    
    % 创建一个与RD图一样大的遮罩矩阵 (1代表目标区，0代表纯噪声区)
    mask = zeros(size(rd_map_mean));
    
    % 自适应划定一个潜在目标核心污染区 
    r_guard = 25; 
    v_guard = 15;
    r_range = max(1, max_r - r_guard) : min(size(rd_map_mean,1), max_r + r_guard);
    v_range = max(1, max_v - v_guard) : min(size(rd_map_mean,2), max_v + v_guard);
    mask(r_range, v_range) = 1;
    
    % 提取彻底远离目标的【纯净噪声区】的能量
    pure_noise_pixels = rd_map_mean(mask == 0);
    mean_noise_power = mean(pure_noise_pixels); % 统计出真实的系统频域噪底
    
    % 计算盲测 Peak SNR
    Metrics.peak_snr_dB = 10 * log10(max_val / mean_noise_power);

    
    %% ================== 2. 盲测 3dB 主瓣网格宽度 ==================
    % 提取最强能量处的距离维切片，并转为 dB 归一化格式
    slice_r = rd_map_mean(:, max_v);
    slice_r_dB = 10 * log10(slice_r / max(slice_r));
    
    % 自适应寻找跌破 -3dB 的网格点数
    % 在归一化谱线上，值大于 -3 dB 的点全算作主瓣核心
    mainlobe_bins = sum(slice_r_dB > -3); 
    
    % 将网格点数换算成物理长度 (米)
    % 距离维轴分辨率：fs * c / (2 * slope * rangeBin)
    res_bin_physics = (Parameter.Fs * Parameter.c) / (2 * Parameter.Slope * Parameter.rangeBin);
    Metrics.mainlobe_3dB_width_metres = mainlobe_bins * res_bin_physics;
    Metrics.mainlobe_3dB_bins = mainlobe_bins;
end