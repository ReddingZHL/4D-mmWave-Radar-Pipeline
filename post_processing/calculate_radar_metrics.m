function CfarDoaMetrics = calculate_radar_metrics(targets, centroids, gates, total_bins)
% CALCULATE_RADAR_METRICS 基于物理空间唯一性占位匹配计算雷达扩展目标的 Pd 和 Pfa
%
% 输入参数:
%   targets : [N_truth x 4] 矩阵，上帝视角真值 [Range, Vel, HorAngle, TargetID]
%   centroids    : [num_clusters x 5] 矩阵，已凝聚的检测质心 [X, Y, Z, V, Total_SNR_dB] 
%   gates           : 结构体，包含三个维度的物理容差门限，例如：
%                     gates.r 距离容差 (米)
%                     gates.v  速度容差 (m/s)
%                     gates.ang 角度容差 (度)
%   total_bins      : 标量，2D-FFT的总网格数（rangeBin * dopplerBin），用于计算Pfa基数
%
% 输出参数:
%   CfarDoaMetrics     : 结构体，包含 CfarDoaMetrics.Pd 、 CfarDoaMetrics.Pfa 、 CfarDoaMetrics.angle_rmse

    %% 1. 数据解析与初始化
    % 初始化输出
    CfarDoaMetrics.Cluster_Pd = 0; % 目标层面的检测概率
    CfarDoaMetrics.Cluster_Pfa = 0; % 目标层面的虚警率
    CfarDoaMetrics.Cluster_False_Alarm_Count = 0; % 没对上真实目标的假簇质心总数
    CfarDoaMetrics.angle_rmse_pure = NaN; % 仅针对成功匹配上的目标的测角误差

    % --- 阶段 1：凝聚真实目标的上帝视角质心 ---
    if isempty(targets)
        num_true_targets = 0;
        True_Centers = [];
    else
        unique_true_ids = unique(targets(:, 4)); 
        num_true_targets = length(unique_true_ids);
        True_Centers = zeros(num_true_targets, 4); % [Range, Vel, Angle, TargetID]
        
        for tIdx = 1:num_true_targets
            curr_true_id = unique_true_ids(tIdx);
            curr_mask = (targets(:, 4) == curr_true_id);
            
            % 上帝视角几何凝聚
            r_mean   = mean(targets(curr_mask, 1));
            v_mean   = mean(targets(curr_mask, 2)); % 航向径向速
            ang_mean = mean(targets(curr_mask, 3));
            
            True_Centers(tIdx, :) = [r_mean, v_mean, ang_mean, curr_true_id];
        end
    end
    
    % --- 阶段 2：检测质心 直角转极坐标 ---
    if isempty(centroids)
        num_detect_clusters = 0;
        Detect_Centers_Radar = [];
    else
        num_detect_clusters = size(centroids, 1);
        Detect_Centers_Radar = zeros(num_detect_clusters, 3); % [Range, Vel, Angle]
        
        for cIdx = 1:num_detect_clusters
            cx = centroids(cIdx, 1);
            cy = centroids(cIdx, 2);
            cv = centroids(cIdx, 4); % 第 4 列是径向速度
            
            c_range = sqrt(cx^2 + cy^2);
            c_angle = atand(cx / cy); % 保持测角定义对齐
            Detect_Centers_Radar(cIdx, :) = [c_range, cv, c_angle];
        end
    end

    % --- 阶段 3：边界防御 ---
    if num_detect_clusters == 0
        CfarDoaMetrics.Cluster_Pd = 0;
        CfarDoaMetrics.Cluster_Pfa = 0;
        return;
    end
    if num_true_targets == 0
        CfarDoaMetrics.Cluster_Pd = 0;
        CfarDoaMetrics.Cluster_Pfa = num_detect_clusters / total_bins; 
        return;
    end
    
    % --- 阶段 4：点迹级别的唯一性占位匹配 ---
    detect_cluster_used = false(num_detect_clusters, 1);
    hit_target_count = 0;
    angle_squared_errors = zeros(num_true_targets, 1);
    
    for i = 1:num_true_targets
        t_r   = True_Centers(i, 1);
        t_v   = True_Centers(i, 2);
        t_ang = True_Centers(i, 3);
        
        % 算残差
        r_diff   = abs(Detect_Centers_Radar(:, 1) - t_r);
        v_diff   = abs(Detect_Centers_Radar(:, 2) - t_v);
        ang_diff = abs(Detect_Centers_Radar(:, 3) - t_ang);
        % 先筛出所有同时满足三维门限、且未被使用的候选点（布尔向量）
        valid_matches = (r_diff <= gates.r) & ...
                        (v_diff <= gates.v) & ...
                        (ang_diff <= gates.ang) & ...
                        (~detect_cluster_used);
                    
        if any(valid_matches)
            find_idx = find(valid_matches); % 找出所有候选点的行号
            % cand_ang_diff = ang_diff(find_idx);
            score = r_diff(find_idx); % 如果只在乎空间位置，可以直接用欧氏距离 (也可综合考虑距离和角度)
            [~, min_local_idx] = min(score); % 寻觅最小值对应的相对位置
            match_idx = find_idx(min_local_idx); % 提炼出全场最优匹配点的真实行号
            
            detect_cluster_used(match_idx) = true; % 占位锁死
            hit_target_count = hit_target_count + 1;
            
            % 测角误差
            angle_squared_errors(hit_target_count) = (Detect_Centers_Radar(match_idx, 3) - t_ang)^2;
        end
    end

    % --- 阶段 5：指标计算 ---
    CfarDoaMetrics.Cluster_Pd = hit_target_count / num_true_targets; % 检出率
    
    if hit_target_count > 0
        CfarDoaMetrics.angle_rmse_pure = sqrt(mean(angle_squared_errors(1:hit_target_count))); % 测角均方误差
    end
    
    false_cluster_count = num_detect_clusters - hit_target_count;
    CfarDoaMetrics.Cluster_False_Alarm_Count = false_cluster_count;
    N_truth_points_total = size(targets, 1);
    N_noise_cells_total = total_bins - N_truth_points_total; % 分母，减去所有真实目标碎点数
    % 分子，虚警簇的个数
    CfarDoaMetrics.Cluster_Pfa = false_cluster_count / max(1, N_noise_cells_total); % 虚警率，这个值通常极小（比如 10^{-5}），它代表“单点发生误判的概率”。
end