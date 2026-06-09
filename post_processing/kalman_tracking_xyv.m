function [tracked_objects, SummaryMetrics] = kalman_tracking_xyv(all_frames_centroids, Parameter, all_frames_truth)
    % =====================================================================
    % 函数功能：全维多普勒卡尔曼跟踪与航迹生命周期管理 (位置+速度 联动更新)
    % 1. 新生航迹具有多普勒速度，消除孵化时延，
    % 2. 提升航迹交叉错乱能力。多目标轨迹发生交叉、并飞、斜切错过时，仅靠距离容易发生航迹掉包，
    % =====================================================================
    % 【速度估计完全不准，可放弃使用速度估计，距离估计不错，可用距离来估计速度】
    dt = Parameter.T_frame; % 帧周期 (秒)
    numFrames = Parameter.numFrames;
    has_truth = (nargin > 2) && ~isempty(all_frames_truth);
    targetnum = 1;Parameter.targetnum;

    %% --- 门限与常数设置 ---
    if targetnum == 1 % 点目标
        MERGE_THRESH = 3.5; % 关联距离门限 (米)。两帧间位移超此距离不予关联 【根据目标大致速度和帧周期可计算帧间位移最大值，略大于即可】
        MERGE_VEL_THRESH  = 2.0;
        INIT_HIT = 2; % 孵化门限：连续 2 帧匹配上，才从暂态激活为正式航迹
        MAX_MISS = 3; %
        EVAL_GATE = 2.5;
    else % 飞机目标
        MERGE_THRESH      = 12.0;  % 空间凝聚门限 (米)：大目标的宏观物理尺寸边界
        MERGE_VEL_THRESH  = 2.0;   % 速度凝聚门限 (米/秒)：防交叉误融合
        INIT_HIT          = 4;     % 孵化打卡门限
        MAX_MISS          = 5;     % 死亡销毁门限
        EVAL_GATE         = 10.0;  % 评估真值匹配时的物理宽容波门 (米)
    end
    
    
    %% --- 卡尔曼矩阵初始化 (4D全状态模型：X, Y, Vx, Vy) ---
    A = [1 0 dt 0; 
         0 1 0  dt;
         0 0 1  0;
         0 0 0  1];
    H = eye(4); % 观测矩阵 含速度，但速度估计很不准，故放弃使用
    H_pos = H(1:2,:); % 仅取位置部分 
    q_sigma = 6.0; 3.0;
    Q = [dt^3/3   0     dt^2/2   0;     % 过程噪声协方差
         0       dt^3/3   0      dt^2/2; 
         dt^2/2   0       dt      0;      
         0       dt^2/2   0       dt] * q_sigma^2; 
    
    R = diag([0.3, 0.3, 1.5, 1.5]); % 测量噪声协方差 保持高置信度
    R_pos = R(1:2, 1:2); % 只取位置测量噪声
    P_init = eye(4) * 2.0; % 增大初始协方差，让新航迹前几帧收敛更快
    
    % 航迹池、指标管理初始化
    tracks = [];  next_id = 1; 
    frame_range_se = []; frame_radialvel_se   = []; frame_ang_se   = []; 
    track_survival_rates = zeros(numFrames, 1); ospa_distances = zeros(numFrames, 1); 
    c_ospa = 10.0; p_ospa = 2;    
    
    %% --- 跨帧时间轴演进 ---
    for frameId = 1:numFrames
        centroids = all_frames_centroids{frameId}; % [x, y, z, 径向v, power]
        % 2. 卡尔曼预测与数据关联阶段
        numMeas = size(centroids, 1);              
        numTracks = length(tracks); 
        for tIdx = 1:numTracks
            tracks(tIdx).x = A * tracks(tIdx).x; % 卡尔曼预测
            tracks(tIdx).P = A * tracks(tIdx).P * A' + Q;
        end
        
        assignment = zeros(numTracks, 1); 
        meas_used = zeros(numMeas, 1);
        % 【根据速度，动态调整自适应 距离门限】
        track_gates = zeros(numTracks, 1);
        for tIdx = 1:numTracks
            if strcmp(tracks(tIdx).status, 'tentative')
                track_gates(tIdx) = 6.0; % 考虑到速度在变，略微放宽新生目标捕获带
            else
                v_est_mag = norm([tracks(tIdx).x(3), tracks(tIdx).x(4)]); % vx vy合成的速度
                track_gates(tIdx) = 1.8 + 0.1 * v_est_mag; % 自适应距离门限
            end
        end

        if numTracks > 0 && numMeas > 0 
            cost_matrix = inf(numTracks, numMeas); 
            for tIdx = 1:numTracks
                pred_state = tracks(tIdx).x; % 预测的 [X, Y, Vx, Vy]'
                for mIdx = 1:numMeas % 测量值
                    pred_range = max(0.1, sqrt(pred_state(1)^2 + pred_state(2)^2)); % 预测range
                    pred_v_radial = (pred_state(1)*pred_state(3) + pred_state(2)*pred_state(4)) / pred_range; % 预测径向vel

                    pos_err = norm(pred_state(1:2) - centroids(mIdx,1:2).'); 
                    vel_err = abs(pred_v_radial - centroids(mIdx, 4));
                    if (pos_err <= MERGE_THRESH) && (vel_err <= MERGE_VEL_THRESH)
                        cost_matrix(tIdx, mIdx) = pos_err + vel_err;
                    else
                        cost_matrix(tIdx, mIdx) = inf;
                    end

                end
            end

            % GNN 贪心择优分配
            while true
                [min_val, min_idx] = min(cost_matrix(:));
                if isinf(min_val), break; end
                
                [t_match, m_match] = ind2sub(size(cost_matrix), min_idx); 
                
                % 空间波门强截断
                mx = centroids(m_match, 1); my = centroids(m_match, 2);
                actual_pos_dist = norm(tracks(t_match).x(1:2) - [mx; my]);
                if actual_pos_dist > track_gates(t_match)
                    cost_matrix(t_match, m_match) = inf; 
                    continue; 
                end
                
                assignment(t_match) = m_match; 
                meas_used(m_match) = 1;        
                cost_matrix(t_match, :) = inf;
                cost_matrix(:, m_match) = inf;
            end
        end
        
        % 3. 卡尔曼更新阶段
        updated_tracks_idx = [];
        for tIdx = 1:numTracks
            m_idx = assignment(tIdx);
            if m_idx > 0 
                cx = centroids(m_idx, 1); cy = centroids(m_idx, 2); vr = centroids(m_idx, 4);
                c_range = max(0.1, sqrt(cx^2 + cy^2));
                % 无论 vr 是正是负，直接将其转化为绝对标量大小
                v_mag = abs(vr); 
                vx_meas = sign(vr) * v_mag * (cx / c_range);
                vy_meas = sign(vr) * v_mag * (cy / c_range);
                
                if targetnum == 1
                    z = [cx; cy; vx_meas; vy_meas];
                    S = H * tracks(tIdx).P * H' + R;
                    K = tracks(tIdx).P * H' / S; 
                    tracks(tIdx).x = tracks(tIdx).x + K * (z - H * tracks(tIdx).x);
                    tracks(tIdx).P = (eye(4) - K * H) * tracks(tIdx).P;
                    
                else
                    %-------用X Y Vx Vy 更新滤波器矩阵----------
                    z = [cx; cy; vr; 0]; 
                    S = H * tracks(tIdx).P * H' + R;
                    K = tracks(tIdx).P * H' / S; 
                    tracks(tIdx).x = tracks(tIdx).x + K * (z - H * tracks(tIdx).x);
                    tracks(tIdx).P = (eye(4) - K * H) * tracks(tIdx).P;
                    %------------------------------------------
                    %--------因速度非恒定，仅用X Y位置进行滤波矩阵更新---------
%                     z_pos = [cx; cy];
%                     % 卡尔曼更新
%                     S = H_pos * tracks(tIdx).P * H_pos' + R_pos;
%                     K = tracks(tIdx).P * H_pos' / S; 
%                     tracks(tIdx).x = tracks(tIdx).x + K * (z_pos - H_pos * tracks(tIdx).x);
%                     tracks(tIdx).P = (eye(4) - K * H_pos) * tracks(tIdx).P;
                    %------------------------------------------
                end
                tracks(tIdx).hit = tracks(tIdx).hit + 1;
                tracks(tIdx).miss = 0; 
                if strcmp(tracks(tIdx).status, 'tentative') && tracks(tIdx).hit >= INIT_HIT
                    tracks(tIdx).status = 'active';
                end
            else
                tracks(tIdx).miss = tracks(tIdx).miss + 1;
                tracks(tIdx).hit = 0; 
            end
            if tracks(tIdx).miss < MAX_MISS
                tracks(tIdx).history = [tracks(tIdx).history, tracks(tIdx).x];
                updated_tracks_idx = [updated_tracks_idx, tIdx]; 
            end
        end
        tracks = tracks(updated_tracks_idx);
        
        
        % 4. 带速度初始化的新航迹孵化
        for mIdx = 1:numMeas
            if meas_used(mIdx) == 0
                % 显式定义结构体
                new_track = struct('id', next_id, ...
                                   'x', [centroids(mIdx, 1); centroids(mIdx, 2); 0; 0], ... % 临时初始化
                                   'P', P_init, ...
                                   'status', 'tentative', ...
                                   'hit', 1, ...
                                   'miss', 0, ...
                                   'history', []); % 固定 2x1 向量
                new_track.id = next_id;
                cx = centroids(mIdx, 1);cy = centroids(mIdx, 2);vr = centroids(mIdx, 4);
                c_range = max(0.1, sqrt(cx^2 + cy^2));
%                 v_mag_init = abs(vr);
                init_vx = vr * (cx / c_range); % 设初始运动朝着/背离雷达直线方向
                init_vy = vr * (cy / c_range);
                
                new_track.x = [cx; cy; init_vx; init_vy]; 
                P_init_adaptive = diag([1.0, 1.0, 100.0, 100.0]); % 初始速度不确定度高
                new_track.P = P_init_adaptive;
                new_track.status = 'tentative'; 
                new_track.hit = 1;new_track.miss = 0;
                new_track.history = new_track.x; 
%                 new_track.first_pos = [cx; cy]; % 记录第一帧的绝对出生位置，用于第二帧进来时差分对账
                tracks = [tracks(:)', new_track];
                next_id = next_id + 1;
            end
        end
        
        %% --- 指标性能评估 ---
        % 【标准OSPA计算的是直角坐标系下的欧式距离，不包含速度。所以这里使用到的是纯距离门限】
        % OSPA 考核的是：“你的雷达在空间位置上有没有漏检目标，或者有没有虚警出幽灵目标。”
        if has_truth
            truth_frame = all_frames_truth{frameId}; % [Range Vel Angle TargetID]
            if ~isempty(truth_frame)
                unique_target_ids = unique(truth_frame(:, 4)); 
                K_true = length(unique_target_ids); 
                X_true = zeros(K_true, 1); Y_true = zeros(K_true, 1); V_true = zeros(K_true, 1);
                for tidx = 1:K_true
                    current_id = unique_target_ids(tidx);
                    id_mask = (truth_frame(:, 4) == current_id); 
                    r_sub = truth_frame(id_mask, 1);
                    v_sub = truth_frame(id_mask, 2);
                    a_sub = truth_frame(id_mask, 3);
                    % 准确计算每个目标当前帧的真实几何中心
                    X_true(tidx) = mean(r_sub .* sind(a_sub));
                    Y_true(tidx) = mean(r_sub .* cosd(a_sub));
                    V_true(tidx) = mean(v_sub); 
                end
            else
                K_true = 0; X_true = []; Y_true = []; V_true = [];
            end
            
            % 提取当前帧所有激活(Active)状态的航迹估计值
            tracks = tracks(:); 
            active_mask = false(length(tracks), 1);
            for kk = 1:length(tracks)
                if strcmp(tracks(kk).status, 'active')
                    active_mask(kk) = true;
                end
            end
            active_tracks = tracks(active_mask);
            M_active = length(active_tracks);
            
            % 统计航迹存活率
            if K_true > 0
                track_survival_rates(frameId) = min(1.0, M_active / K_true);
            else
                % 在没有目标时，实际也没捕获到目标时，给系统打满分 1.0
                track_survival_rates(frameId) = 1.0;
            end
            
            X_est = zeros(M_active, 1); Y_est = zeros(M_active, 1); 
            for t = 1:M_active
                X_est(t) = active_tracks(t).x(1);
                Y_est(t) = active_tracks(t).x(2);
            end
            
            if M_active > 0 && K_true > 0
                for t = 1:M_active
                    dist_to_truths = sqrt((X_est(t) - X_true).^2 + (Y_est(t) - Y_true).^2);
                    [min_d, matched_true_idx] = min(dist_to_truths);
                    
                    if min_d <= EVAL_GATE % 距离门限
                        r_est = norm([X_est(t), Y_est(t)]);
                        ang_est = atan2d(X_est(t), Y_est(t));
                        
                        r_true_val = sqrt(X_true(matched_true_idx)^2 + Y_true(matched_true_idx)^2);
                        ang_true_val = atan2d(X_true(matched_true_idx), Y_true(matched_true_idx));
                        tx_e = active_tracks(t).x(1); % 航迹当前估计的 X
                        ty_e = active_tracks(t).x(2); % 航迹当前估计的 Y
                        tvx = active_tracks(t).x(3); % 航迹当前估计的 Vx
                        tvy = active_tracks(t).x(4); % 航迹当前估计的 Vy
                        track_range = max(0.1, sqrt(tx_e^2 + ty_e^2));
                        v_est_radial = (tx_e * tvx + ty_e * tvy) / track_range; % 目标当前对雷达的径向运动分量
                        v_true_radial = V_true(matched_true_idx);
                        
                        frame_range_se = [frame_range_se; (r_est - r_true_val)^2];
                        frame_ang_se   = [frame_ang_se;   (ang_est - ang_true_val)^2];
                        frame_radialvel_se = [frame_radialvel_se; (v_est_radial - v_true_radial)^2];
                    end
                end
            end
            Set_True = [X_true, Y_true];
            Set_Est  = [X_est,  Y_est];
            ospa_distances(frameId) = calculate_ospa(Set_True, Set_Est, c_ospa, p_ospa);
        end
    end
    
    % 指标平滑提取
    SummaryMetrics.Range_RMSE = sqrt(mean(frame_range_se));
    SummaryMetrics.Vel_RMSE   = sqrt(mean(frame_radialvel_se));
    SummaryMetrics.Angle_RMSE = sqrt(mean(frame_ang_se));
    SummaryMetrics.Mean_Track_Survival_Rate = mean(track_survival_rates(INIT_HIT:end)); 
    SummaryMetrics.Mean_OSPA  = mean(ospa_distances);
    
    %% --- 输出封装：只向用户输出真正激活过的、高质量的合法目标航迹 ---
    valid_count = 0;
    tracked_objects = struct([]);
    for tIdx = 1:length(tracks)
        if strcmp(tracks(tIdx).status, 'active') || tracks(tIdx).hit >= INIT_HIT
            valid_count = valid_count + 1;
            tracked_objects(valid_count).id = tracks(tIdx).id;
            tracked_objects(valid_count).trajectory = tracks(tIdx).history; 
        end
    end

    % =======================================================================
    % 在函数返回前，自动打印终局物理速度分解账目
    % =======================================================================
    if ~isempty(tracked_objects)
        fprintf('\n================== 滤波状态终局速度对账控制台 ==================\n');
        for i = 1:length(tracked_objects)
            traj = tracked_objects(i).trajectory; % 4 x N 的历史矩阵
            if ~isempty(traj)
                first_est = traj(:, 1);
                last_est = traj(:, end); 
                tx_s = first_est(1); ty_s = first_est(2);
                tx_e = last_est(1);  ty_e = last_est(2); 
                L = sqrt((tx_e - tx_s)^2 +(ty_e - ty_s)^2);
                tvx = last_est(3); tvy = last_est(4);
                
                % 1. 计算绝对真实总航速 (二维合成标量，恒为正)
                total_speed = L / (length(traj) * dt);  %norm([tvx, tvy]); % 
                
                % 2. 还原带物理符号的径向运动速度 (正代表远离，负代表靠近雷达)
                r_est = max(0.1, sqrt(tx_e^2 + ty_e^2));
                radial_speed_with_sign = (tx_e*tvx + ty_e*tvy) / r_est;
                
                % 3. 规范化打印输出
                fprintf(' 目标 ID [%d] 最后一帧-> \n', tracked_objects(i).id);
                fprintf('   [绝对总航速(大小)] : %.4f m/s \n', total_speed);
                fprintf('   [径向航速(带符号)] : %.4f m/s \n', radial_speed_with_sign);
                fprintf('   [直角分量状态解算] : Vx = %.4f m/s, Vy = %.4f m/s \n', tvx, tvy);
                fprintf(' ----------------------------------------------------------------\n');
            end
        end
        fprintf('==================================================================\n\n');
    end
    % =======================================================================
end

% 指标：OSPA（Optimal Subpattern Assignment，最优子模式分配）：既考核距离误差，又严惩虚警和漏检（势误差）
function ospa_dist = calculate_ospa(X, Y, c, p)
    K = size(X, 1); M = size(Y, 1);
    if K == 0 && M == 0, ospa_dist = 0; return;
    elseif K == 0 || M == 0, ospa_dist = c; return;
    end
    if K > M
        temp = X; X = Y; Y = temp;
        K = size(X, 1); M = size(Y, 1);
    end
    dist_matrix = zeros(K, M);
    for i = 1:K
        for j = 1:M
            dist_matrix(i, j) = min(c, norm(X(i, :) - Y(j, :)));
        end
    end
    total_loc_error_p = 0;
    for i = 1:K
        [~, min_idx] = min(dist_matrix(i, :));
        total_loc_error_p = total_loc_error_p + dist_matrix(i, min_idx)^p;
        dist_matrix(:, min_idx) = inf; 
    end
    cardinality_penalty = (c^p) * (M - K);
    ospa_dist = ((total_loc_error_p + cardinality_penalty) / M)^(1/p);
end