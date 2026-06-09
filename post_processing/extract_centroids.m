function [centroids, ClusterMetrics] = extract_centroids(data_cluster,targets)
    % =====================================================================
    % 函数功能：提取当前帧点云聚类后的质心点 (Measurement Fusion)
    % 输入参数：
    %   data_cluster - [N x 6] 矩阵，聚类后的点云数据
    %                  第 1, 2, 3 列: X, Y, Z (直角坐标)
    %                  第 4 列: Velocity (多普勒速度/径向速度)
    %                  第 5 列: SNR / Power (信号能量)
    %                  第 6 列: Cluster ID (聚类标签，0或-1代表噪声)
    %   targets - [M x 4]矩阵，真实目标质心 Range, Vel, HorAngle, TargetID
    % 输出参数：
    %   centroids - [M x 5] 矩阵，凝聚后的多目标质心点迹列表 (M 为目标个数)
    %                   每一行格式为: [X, Y, Z, Velocity, Total_Power]
    %   ClusterMetrics  - 结构体数组 (长度为 M)，包含每个簇的量化评估指标：
    %                   .density_vol  - 空间绝对密度 (点数/立方米)
    %                   .compactness  - 相对紧凑度 (碎点到质心的平均欧氏距离，米)
    %                   .center_error - 凝聚质心与飞机真实几何中心的物理偏差 (米)
    % =====================================================================
    
    ClusterMetrics.avg_density    = 0.0;
    ClusterMetrics.avg_compact    = 0.0;
    ClusterMetrics.avg_center_err = NaN; % 测距测角偏差没配上时按学术惯例挂 NaN

    % 1. 安全检查：如果当前帧没有任何点云，直接返回空
    if isempty(data_cluster)
        centroids = [];
        ClusterMetrics = [];
        return;
    end
    
    % 2. 提取所有唯一的聚类标签
    all_ids = data_cluster(:, 6); % Cluster ID
    unique_ids = unique(all_ids); % 聚类数
    
    % 3. 剔除噪声标签 (DBSCAN 中通常 0 或 -1 代表无法成类的孤立噪声点)
    valid_ids = unique_ids(unique_ids > 0); 
    num_clusters = length(valid_ids); % 有效聚类数
    
    % 4. 如果没有合法的聚类目标，直接返回空
    if num_clusters == 0
        centroids = [];
        ClusterMetrics = [];
        return;
    end
    
    % 5. 合成每个有效聚类的质心点
    centroids = zeros(num_clusters, 5); % [X, Y, Z, Velocity, Total_Power]
    W_points = zeros(num_clusters, 1); % 记录每个有效聚类中的碎点数
    % 预分配指标结构体数组，单个聚类的
    TermMetrics = struct('density_vol', cell(num_clusters, 1), ...
                         'compactness', cell(num_clusters, 1), ...
                         'center_error', cell(num_clusters, 1));
    for ii = 1:num_clusters % 遍历每一个有效聚类簇，将属于同一个标签的所有碎点的坐标、速度进行算术平均，能量（SNR）累加，
        current_id = valid_ids(ii);
        points_in_cluster = data_cluster(all_ids == current_id, :); % 提取属于当前聚类 ID 的所有碎点
        N_points = size(points_in_cluster, 1); % 当前有效聚类的碎点数
        W_points(ii) = N_points; % 记录该簇碎点数

        mean_X = mean(points_in_cluster(:, 1));
        mean_Y = mean(points_in_cluster(:, 2));
        mean_Z = mean(points_in_cluster(:, 3));
        mean_V = mean(points_in_cluster(:, 4)); % 速度平均
        
        snr_db_points = points_in_cluster(:, 5); % 提取所有碎点的单点 SNR (dB)
        sum_snr_linear = sum(10.^(snr_db_points ./ 10)); % 转回再求和，【不能将db值简单相加】
        total_snr_db = 10 * log10(sum_snr_linear); % 能量和（不能将db值简单相加）
        centroids(ii, :) = [mean_X, mean_Y, mean_Z, mean_V, total_snr_db]; % 簇质心

        %% ====== 指标一：点云空间绝对密度 (Volume Density) ======
        
        % 包围单簇的三维立方体长宽高
        dx = max(points_in_cluster(:, 1)) - min(points_in_cluster(:, 1));
        dy = max(points_in_cluster(:, 2)) - min(points_in_cluster(:, 2));
        dz = max(points_in_cluster(:, 3)) - min(points_in_cluster(:, 3));
        volume = max(0.1, dx) * max(0.1, dy) * max(0.1, dz); % 立方体体积，规避分母为 0 风险
        TermMetrics(ii).density_vol = N_points / volume; % 单位: points / m^3
        
        %% ====== 指标二：点云相对紧凑度 (Compactness) ======
        % 计算该簇内所有碎点到凝聚质心中心的平均欧氏距离
        TermMetrics(ii).compactness = mean(sqrt((points_in_cluster(:, 1) - mean_X).^2 + ...
                                   (points_in_cluster(:, 2) - mean_Y).^2 + ...
                                   (points_in_cluster(:, 3) - mean_Z).^2) );
         
        %% ====== 指标三：【多目标最近邻真值偏差】 (Center Error) ======
        has_truth = (nargin > 1) && ~isempty(targets);
        if has_truth
            unique_true_ids = unique(targets(:, 4)); % 真实大目标数
            num_true_targets = length(unique_true_ids);
            % 动态建立每个大目标的直角坐标几何中心矩阵 [X Y Z]
            True_Centers = zeros(num_true_targets, 3);
            for tIdx = 1:num_true_targets % 遍历每个大目标
                curr_true_id = unique_true_ids(tIdx);
                curr_target_mask = (targets(:, 4) == curr_true_id); % 获取每个大目标的所有碎点信息
                
                % 分别计算这架特定飞机的碎点极坐标均值
                r_mean   = mean(targets(curr_target_mask, 1)); % 所有碎点距离均值
                ang_mean = mean(targets(curr_target_mask, 3)); % 所有碎点角度均值
                
                % 转换到直角坐标，得到该大目标的纯正几何中心
                True_Centers(tIdx, 1) = r_mean * sind(ang_mean); % 大目标质心的X值
                True_Centers(tIdx, 2) = r_mean * cosd(ang_mean); % 大目标质心的Y值
                True_Centers(tIdx, 3) = 0; % 2D平面高度置零 大目标质心的Z值置零
            end

            % 当前簇质心与最接近的真实目标质心距离偏差
            TermMetrics(ii).center_error = min( sqrt((mean_X - True_Centers(:, 1)).^2 + ...
                                     (mean_Y - True_Centers(:, 2)).^2 + ...
                                     (mean_Z - True_Centers(:, 3)).^2) );
        else
            TermMetrics(ii).center_error = NaN;
        end
    end

    % 【多目标情况下不能简单的用所有真实目标的均值来计算】
    raw_densities   = [TermMetrics.density_vol]';
    raw_compactness = [TermMetrics.compactness]';
    raw_errors      = [TermMetrics.center_error]';
    W_norm = W_points(:) / sum(W_points); % 归一化点数权重 (大目标权重高，小噪声权重低)
    if num_clusters > 0
        ClusterMetrics.avg_density = W_norm' * raw_densities;
        ClusterMetrics.avg_compact = W_norm' * raw_compactness;
        
        % 针对中心偏差进行无效值(NaN)过滤防护
        valid_err_mask = ~isnan(raw_errors);
        if any(valid_err_mask)
            W_err = W_points(valid_err_mask);
            W_err_norm = W_err(:) / sum(W_err);
            sub_errors = raw_errors(valid_err_mask);
            
            % 标量化相乘
            ClusterMetrics.avg_center_err = W_err_norm' * sub_errors(:);
        else
            ClusterMetrics.avg_center_err = NaN;
        end
    else
        ClusterMetrics.avg_density = TermMetrics.density_vol;
        ClusterMetrics.avg_compact = TermMetrics.compactness;
        ClusterMetrics.avg_center_err = TermMetrics.center_error;
    end

end