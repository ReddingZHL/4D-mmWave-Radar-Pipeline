% 点的聚类
% 输入1：聚类算法 1：DBSCAN 2：K-Means 3：X-Means
% 输入2：数据集     data (N x 5) [x, y, z, 径向v, power, Cluster_ID]
% 输入3：聚类半径   EPS
% 输入4：最小数据量 minPoints
% 输出1：聚类数据   data_cluster
function data_cluster = points_clustering(clusterMethod, data, targetnum, eps, minpoints)
    % data: N x 5 矩阵 [x, y, z, v, power]
    if nargin <= 3
        if targetnum == 1 % 点目标
            eps = 1.5;      % 最大搜索半径
            minpoints = 5; % 最小点数
        else % 飞机目标
            eps = 3.5;      % 根据飞机大小设置，不然单个飞机会出现多个簇
            minpoints = 10;
        end
    end
    %【除了用距离，也可以加上多普勒速度进行聚类，在这里，飞机如果侧飞，不同部位的径向速度并不一致，所以暂时先不加上速度聚类】
    features = data(:, 1:3); % 提取用于聚类的特征维度 [x,y,z] 方便做 norm 距离计算
    
    if clusterMethod == 1 % DBSCAN
        % 1-----调用MATLAB内置的高效函数dbscan
        % idx 返回每个点对应的簇标签。如果点被判定为噪声/离群点，idx 会被赋予 -1
%         idx = dbscan(features, eps, minpoints); % [dataNum, 1]
        
        % 2-----调用实现的my_dbscan函数
        idx = my_dbscan(features, eps, minpoints);

        % 剔除噪声点（离群点）
        outliers_mask = (idx == -1);
        data(outliers_mask, :) = [];
        idx(outliers_mask) = [];
        data_cluster = [data, idx]; % [dataNum, 6] 最后一列 簇标签idx（聚类ID）

    elseif clusterMethod == 2 % K-Means
        % 局限性说明：K-Means 无法识别噪声点，且必须输入固定的 K 值
        % 在实际多目标雷达管线中，通常不建议直接用标准 K-Means
        % 假设场上有 3 个目标
        fixed_K = 3;
        if size(features, 1) < fixed_K
            idx = ones(size(features, 1), 1); % 点数不够分，塞进一类
        else
            % 调用matlab内置的 kmeans，'Replicates', 3 表示尝试 3 次初始化取最优，防止陷入局部极小值
            idx = kmeans(features, fixed_K, 'Replicates', 3, 'Display', 'off');
        end
        data_cluster = [data, idx];
        
    elseif clusterMethod == 3 % X-Means 
        % 自适应K值，从 K_min 开始，通过 BIC 评分自动分裂，直到 K_max
        k_min = 1;
        k_max = 8; % 限制单帧最大可能的目标数
        
        idx = run_xmeans(features, k_min, k_max, targetnum);
        data_cluster = [data, idx];

    end
%     plot_cluster_birdseye(data_cluster, clusterMethod);

end

function idx = my_dbscan(X, epsilon, minPts)
    % 输入:
    %   X       - N x 2 或 N x 3 的点云/点迹坐标矩阵
    %   epsilon - 邻域半径（距离门限）
    %   minPts  - 成为核心点所需的最小邻居数（含自身）
    % 输出:
    %   idx     - N x 1 向量，记录每个点所属的聚类ID（0代表噪声/野点）

    numPoints = size(X, 1);
    idx = zeros(numPoints, 1); % 初始化：0 代表未分类（Unvisited / Noise）
    clusterId = 0;             % 聚类标签计数器

    % 步骤 1：利用 MATLAB 矩阵矩阵秒算“两两欧氏距离矩阵” (N x N)
    % pdist2（Pairwise Distance Matrix，成对距离矩阵计算） 比双重 for 循环快几百倍
    % D = pdist2(X, Y)：矩阵D的行数等于X的点数，列数等于Y的点数
    distMatrix = pdist2(X, X); % [N, N]对称矩阵，对角线全0（自己到自己的距离是0）

    % 步骤 2：主循环，遍历每一个点
    for i = 1:numPoints
        % 如果这个点已经被某个群吞并了，直接跳过
        if idx(i) ~= 0
            continue;
        end
        
        % 寻找当前点的所有邻居索引
        neighbors = find(distMatrix(i, :) <= epsilon);
        
        % 步骤 3：密度检查
        if length(neighbors) < minPts
            % 邻居太少，暂时标记为噪声（0），等别人来吞并它
            idx(i) = 0; 
        else
            % 邻居足够，恭喜你成为核心点！开辟一个新帮派
            clusterId = clusterId + 1;
            idx(i) = clusterId;
            
            % 开始“顺藤摸瓜”：把邻居放入一个动态队列中去扩展
            queue = neighbors;
            queue(queue == i) = []; % 队列中移出自己
            
            while ~isempty(queue)
                currentPoint = queue(1); % 取出队列第一个点
                queue(1) = [];           % 出队
                
                % 如果这个点之前被误判为噪声，现在把它收编进当前帮派
                if idx(currentPoint) == 0
                    idx(currentPoint) = clusterId;
                end
                
                % 如果这个点已经被其他正常帮派收编了，跳过它
                % 边界点的归属采用"先到先得"的原则，允许微小的边界模糊，比强行去合并两个独立的物理目标安全得多
                if idx(currentPoint) ~= clusterId && idx(currentPoint) > 0
                    continue;
                end
                
                % 关键扩展：检查这个邻居是不是也是一个核心点
                currentNeighbors = find(distMatrix(currentPoint, :) <= epsilon);
                
                % 如果它也是核心点（有带小弟的能力）
                if length(currentNeighbors) >= minPts
                    % 它的所有未分类小弟，都要加入我们的考察队列
                    for k = 1:length(currentNeighbors)
                        neighborPoint = currentNeighbors(k);
                        if idx(neighborPoint) == 0 % 未访问过或曾是噪声
                            idx(neighborPoint) = clusterId; % 现场收编
                            queue = [queue, neighborPoint]; %#ok<AGROW> % 加入扩展队列
                        end
                    end
                end
            end % while 队列扩展结束
        end
    end % for 主循环结束
end


function idx = run_xmeans(X, k_min, k_max, targetnum)
    [N, M] = size(X); % N个点，M个维度(3)
    if N <= k_min
        idx = ones(N, 1); % 塞进一类
        return;
    end
    
    % 初始状态：从 k_min 开始做初始 K-means
    current_K = k_min;
    idx = kmeans(X, current_K, 'Replicates', 2, 'Display', 'off');
    
    % 开始迭代尝试分裂
    while current_K < k_max
        changed = false;
        new_idx = idx;
        offset = 0; % 用于修正分裂后的标签偏移
        
        for i = 1:current_K
            % 抠出当前第 i 个簇的所有点
            cluster_mask = (idx == i);
            X_sub = X(cluster_mask, :);
            N_sub = size(X_sub, 1);
            
            % 如果这个簇点数太少（少于4个点），就没必要再分裂了
            if N_sub < 4
                new_idx(cluster_mask) = i + offset;
                continue;
            end
            
            % 【物理几何尺寸叫停机制】
            % 如果一个簇在 X（距离）、Y（横向距离）、Z（高度）方向上的最大跨度已经小于你雷达的物理分辨率
            % （或者一个合理的物理阈值，比如 0.3米，那它在物理上绝对已经是一个不可分割的“单体”了。
            % 计算这个簇在 X, Y, Z 三个物理维度上的最大跨度（几何尺寸）
            cluster_span = max(X_sub, [], 1) - min(X_sub, [], 1);
            % 如果它在任何一个方向的跨度都小于 0.6 米（说明已经聚成一团完美点迹了）
            % 直接强行跳过，不准它参与任何 BIC 计算和分列！
            if targetnum == 1 % 点目标
                est_siz = 0.6;
            else
                est_siz = 9.0; % 1目标 7.5;% 2目标 5.5;% 4目标 3.5;% 5 目标
            end
            if all(cluster_span < est_siz) % 这个值经过几次调试得到
                new_idx(cluster_mask) = i + offset;
                continue; % 强行保命，不许切！
            end

            
            % 计算不分裂时的原始 BIC 分数
            bic_parent = calculate_bic(X_sub, 1);
            
            % 尝试对这个簇进行 2-Means 强行一分为二
            try
                sub_idx = kmeans(X_sub, 2, 'Replicates', 2, 'Display', 'off');
                % 计算分裂后的新 BIC 分数
                bic_child = calculate_bic(X_sub, sub_idx);
            catch
                bic_child = -Inf; % 如果 kmeans 报错，放弃分裂
            end
            
            % 判断标准：如果分裂后的 BIC 分数更高，说明切开更好！
            if bic_child > bic_parent
                % 接受分裂：将原先的标签 i，变成两个互不干扰的新标签
                sub_idx_adjusted = sub_idx;
                sub_idx_adjusted(sub_idx == 1) = i + offset;
                sub_idx_adjusted(sub_idx == 2) = current_K + offset + 1;
                
                new_idx(cluster_mask) = sub_idx_adjusted;
                offset = offset + 1; % 增加了一个新簇
                changed = true;
            else
                % 拒绝分裂：保持原样，仅做标签修正
                new_idx(cluster_mask) = i + offset;
            end
        end
        
        % 更新总体局势
        idx = new_idx;
        current_K = current_K + offset;
        
        % 如果遍历了所有簇，没有任何一个簇想继续分裂，或者达到了最大 K 限制，就退出
        if ~changed || current_K >= k_max
            break;
        end
    end
    
    % 规整标签，使其连续（例如把 1, 3, 7 变成 1, 2, 3）
    [~, ~, idx] = unique(idx);
end

%% ================== 📐 数学工具：BIC (贝叶斯信息准则) 计算器 ==================
function bic = calculate_bic(X, idx)
    % 该函数专门用于评估在 X 数据集上，当前的聚类划分 idx 分数是多少
    [N, M] = size(X);
    K = length(unique(idx));
    
    % 计算全场方差 (估计最大似然)
    variance = 0;
    for i = 1:K
        cluster_mask = (idx == i);
        X_sub = X(cluster_mask, :);
        if size(X_sub, 1) > 1
            centroid = mean(X_sub, 1);
            variance = variance + sum(sum((X_sub - centroid).^2));
        end
    end
    variance = variance / (N - K);
    
    if variance <= eps
        bic = -Inf; % 方差为0说明重合，属于无效划分
        return;
    end
    
    % 计算对数似然度 (Log-Likelihood)
    log_likelihood = 0;
    for i = 1:K
        cluster_mask = (idx == i);
        N_n = sum(cluster_mask);
        if N_n > 0
            % X-Means 论文中的标准极大似然估计公式
            log_likelihood = log_likelihood + N_n * log(N_n) - N_n * log(N) ...
                             - (N_n * M / 2) * log(2 * pi * variance) ...
                             - (N_n - 1) / 2;
        end
    end
    
    % 计算自由参数个数 (自由度)
    % 每个中心有 M 维坐标，加上一个方差项
    num_parameters = K * (M + 1);
    
    % BIC 终极公式：得分 = 似然度 - 参数复杂度惩罚项
    % 惩罚项能有效防止算法无节制地分裂成无数个极小的簇
    penalty_factor = 2.5;
    bic = log_likelihood - penalty_factor*(num_parameters / 2) * log(N);
end


function plot_cluster_birdseye(data_cluster, clusterMethod)
    % plot_cluster_birdseye - 绘制雷达点云聚类的二维鸟瞰图
    %
    % 输入参数:
    %   data_cluster : N x 6 矩阵 [x, y, z, v, power, cluster_id]
    %                  最后一列必须是聚类算法输出的有效簇标签 idx
    %   clusterMethod: 标量，1：DBSCAN 2：K-Means 3：X-Means (仅用于标题字符串切换)
    
    % --- 1. 安全防呆检查 ---
    if isempty(data_cluster)
        warning("输入的数据集为空，取消图表绘制！");
        return;
    end
    
    if size(data_cluster, 2) < 6
        error("输入数据格式不合理！矩阵至少需要6列 [x, y, z, v, power, cluster_id]。");
    end
    
    % --- 2. 提取聚类标签并去重 ---
    idx = data_cluster(:, 6); 
    unique_clusters = unique(idx); 
    num_clusters = length(unique_clusters); 
    
    % --- 3. 核心绘图引擎 ---
    if ~isempty(unique_clusters)
        % 创建独立画布，设置合适的分辨率和宽高比
        figure('Color', 'w', 'Position', [200, 200, 800, 650]); 
        hold on; box on; grid on;
        
        % 定义不同独立目标的图形样式库（防止目标太多导致标记重复）
        marker_library = {'o', 's', '^', 'd', 'v', '>', '<', 'p', 'h'};
        
        for ii = 1:num_clusters
            cluster_id = unique_clusters(ii);
            
            % 筛选出当前目标集群的所有散射点迹
            cluster_mask = (data_cluster(:, 6) == cluster_id);
            cluster_data = data_cluster(cluster_mask, :);
            
            x_val = cluster_data(:, 1);
            y_val = cluster_data(:, 2);
            v_val = cluster_data(:, 4); % 提取第四列的物理速度
            
            % 根据当前的循环序数自适应选取图形形状
            m_style = marker_library{mod(ii-1, length(marker_library)) + 1};
            
            % 使用 scatter 绘制二维投影，通过 v_val 实现动态颜色映射
            scatter(x_val, y_val, 45, v_val, m_style, 'filled', ...
                    'LineWidth', 1, 'MarkerEdgeColor', 'k', ...
                    'DisplayName', ['目标集群 ', num2str(cluster_id)]);
        end
        
        % --- 4. 视角锁定与几何形状保真对账 ---
        view(0, 90); % 强制锁定正上方绝对俯视视角
        axis equal;  % 强行保持一比一物理横纵比，防止飞机等扩展目标形状被拉伸变形
        
        % --- 5. 色彩轴 (Colorbar) 渲染配置 ---
        c = colorbar;
        c.Label.String = '雷达径向速度 (m/s)';
        c.Label.FontSize = 11;
        colormap(jet); % jet 渲染器：蓝色靠近（负速），红色远离（正速）
        
        % 动态界定速度色彩轴的显示范围，防止全场单一速度时 colorbar 报错
        v_min = min(data_cluster(:, 4));
        v_max = max(data_cluster(:, 4));
        if v_max ~= v_min
            caxis([v_min - 1, v_max + 1]);
        end
        
        % --- 6. 标签与图例张贴 ---
        method_names = {'DBSCAN', 'K-Means', 'X-Means'};
        if clusterMethod >= 1 && clusterMethod <= 3
            method_str = method_names{clusterMethod};
        else
            method_str = '未知算法';
        end
        
        title(sprintf('%s 聚类鸟瞰图 (颜色表示速度)', method_str), 'FontSize', 12);
        xlabel('X 轴 (横向距离 m)', 'FontSize', 11);
        ylabel('Y 轴 (纵向距离 m)', 'FontSize', 11);
        legend('Location', 'best');
    else
        warning("数据集中未检测到任何有效分类标签，请检查前级聚类算法！");
    end
end
