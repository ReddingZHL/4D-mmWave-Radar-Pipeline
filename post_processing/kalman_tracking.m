function [tracked_objects, SummaryMetrics] = kalman_tracking(all_frames_centroids, Parameter, all_frames_truth)
    % =====================================================================
    % 函数功能：多目标卡尔曼跟踪与航迹生命周期管理 (GNN（Global Nearest Neighbor） 关联 + CV（constant velocity）模型)
    % 输入参数：
    %   all_frames_centroids - cell 数组 [numFrames x 1]，每元胞为当前帧 [M x 5] 质心
    %                          (M 行代表 M 个点迹，列为: [X, Y, Z, V, Power])
    %   Parameter            - 雷达参数结构体
    %   all_frames_truth     - cell 数组 [numFrames x 1]，每元胞为上帝视角物理真值 [K x 4] (Range, Speed, Angle, TargetID)
    % 输出参数：
    %   tracked_objects      - 结构体数组，包含所有演进成功的航迹历史
    %   SummaryMetrics       - 结构体，包含全生命周期的全方位多目标跟踪指标统计
    % =====================================================================
    targetnum = Parameter.targetnum;
    dt = Parameter.T_frame; % 帧周期 (秒)
    numFrames = Parameter.numFrames;
    has_truth = (nargin > 2) && ~isempty(all_frames_truth);
    
    %% --- 门限与常数设置 ---
    if targetnum == 1 % 1-点目标
        GATE_THRESH = 1; % 关联距离门限 (米)。两帧间位移超此距离不予关联 【根据目标大致速度和帧周期可计算帧间位移最大值，略大于即可】
        INIT_HIT = 2; % 孵化门限：连续 2 帧匹配上，才从暂态激活为正式航迹
        MAX_MISS = 3; % 销毁门限：连续 3 帧未匹配上，该航迹判定死亡
        EVAL_GATE = 2.5;
    else % 飞机
        GATE_THRESH = 2;
        INIT_HIT = 4;2; 
        MAX_MISS = 5;3; 
        EVAL_GATE = 10;
    end
    
    
    %% --- 卡尔曼矩阵初始化 (恒定速度 CV 模型) ---
    % 状态向量 x = [X, Y, Vx, Vy]'
    A = [1 0 dt 0; % 状态转移矩阵：X新=X+Vx*dt, Y新=Y+Vy*dt, Vx新=Vx, Vy新=Vy
         0 1 0  dt;
         0 0 1  0;
         0 0 0  1];
     
    H = [1 0 0 0;  % 观测矩阵 (直角坐标系下，测量值只有 X 和 Y)
         0 1 0 0];
    
    % 过程噪声协方差 Q (对速度扰动的物理建模，值越小航迹越平滑但对拐弯响应慢)
    q_sigma = 6.0; 5.5;3;0.5; % 【现实中的人/车行走存在微小的加速度抖动】设定加速度是零均值，方差为q_sigma的白噪声
    Q = [dt^3/3   0     dt^2/2   0;     % dt^4/4*q_sigma^2：位置项的方差  dt^3/2*q_sigma^2：位置与速度的协方差
         0       dt^3/3   0      dt^2/2; % 0：假设X轴和Y轴的运动扰动是相互独立的
         dt^2/2   0       dt      0;      % dt*q_sigma^2是速度项的方差
         0       dt^2/2   0       dt] * q_sigma^2; 
     
    % 测量噪声协方差 R (雷达测量本身的空间误差)
    R = [0.3 0;  
         0    0.3];
    
    % 初始状态协方差 P_init 衡量对当前状态估计的不确定性
    % 对角线上元素分别代表x,y,Vx,Vy对自身的方差，值越大，说明滤波器对当前这个数值越没信心
    P_init = eye(4) * 1.0; 
%     P_init = diag([1.5, 1.5, 6.0, 6.0]); % 初始不确定性（速度方向给更大不确定性）

    %% --- 航迹列表初始化 ---
    % 航迹结构体设计
    % tracks(k).id       : 独一无二的航迹号，确定每帧的同一个目标【当判定出一个无法和现有航迹对上的"野点"，且连续INIT_HIT帧稳定存在，系统就会觉得来新目标了，将next_id发给这个新点作为它的终身代表，然后next_id加1，等待下一个新点。当这个新点目标连续MAX_MISS帧漏检，这个id就会被抹去且永远不会被回收利用，下一个再进来的新人，会被派发全新的id】
    % tracks(k).x        : 卡尔曼当前状态 [4 x 1]    [X, Y, Vx, Vy]
    % tracks(k).P        : 卡尔曼当前协方差 [4 x 4]  不确定性
    % tracks(k).status   : 'tentative' (暂态) 或 'active' (激活-连续两帧稳定存在)
    % tracks(k).hit      : 连续匹配成功计数（命中），目标在所有帧的连续打卡记录
    % tracks(k).miss     : 连续漏检计数（漏检），与命中此消彼长，一旦一帧出现漏检，hit强制清零
    % tracks(k).history  : 存储卡尔曼滤波平滑后的 [X, Y, Vx, Vy] 矩阵，若目标存活了50帧，那么history最终就是一个[4,50]的矩阵 
    tracks = struct('id', {}, 'x', {}, 'P', {}, 'status', {}, 'hit', {}, 'miss', {}, 'history', {});
    next_id = 1; % 航迹 ID 计数器
    
    frame_range_se = []; % 存储所有成功关联上的点迹的 距离误差平方
    frame_vel_se   = []; % 存储 速度误差平方
    frame_ang_se   = []; % 存储 角度误差平方
    track_survival_rates = zeros(numFrames, 1); % 存活率历史
    ospa_distances = zeros(numFrames, 1); % OSPA距离历史
    
    % OSPA 核心惩罚因子设置
    c_ospa = 10.0; % 距离截止惩罚常数 (米)，目标跟丢或多跟一个直接罚 10 米
    p_ospa = 2;    % OSPA阶数，对过大误差具有更敏感的惩罚

    %% --- 核心逻辑：跨帧时间轴演进 ---
    for frameId = 1:numFrames
        centroids = all_frames_centroids{frameId}; % 提取当前帧的所有聚类簇质心
        numMeas = size(centroids, 1);              % 当前帧观测到的目标数（聚类 簇数）
        numTracks = length(tracks);                % 现有存活航迹数，只有出现不匹配原有存活的任何一条航迹的新点，这个numTracks才会增加

        % 【第1帧，系统只孵化新目标航迹】
        % 【第2帧，卡尔曼滤波器利用上一帧记录的位置和A矩阵，往未来推算一步A*tracks(tIdx).x，
        %  这个预测位置将与第二帧雷达测到的点迹进行GNN距离配对】
        
        % 1. 卡尔曼预测阶段 (Prediction)  所有存活航迹根据 CV 模型往未来推算一步
        for tIdx = 1:numTracks
            tracks(tIdx).x = A * tracks(tIdx).x;
            tracks(tIdx).P = A * tracks(tIdx).P * A' + Q;
        end
        
        if targetnum == 1 % 1-点目标
            % 2. 数据关联阶段 (Data Association - GNN 最近邻算法)
            assignment = zeros(numTracks, 1); % 记录当前帧的每一个存活航迹和雷达在当前帧测到的第几个测量点关联成功
            meas_used = zeros(numMeas, 1);    % 标记当前雷达测到的每个聚类质心，是否与存活航迹关联成功了。【用途-判官：将未被关联的测量点，孵化新航迹】
            
            if numTracks > 0 && numMeas > 0 % 要有存活航迹，且还有新测量点，才可进行接下来的匹配
                % 计算距离代价矩阵 (Distance Cost Matrix)
                cost_matrix = inf(numTracks, numMeas); % inf为正无穷大
                for tIdx = 1:numTracks
                    pred_pos = H * tracks(tIdx).x; % 预测的 [X, Y]'
                    for mIdx = 1:numMeas
                        meas_pos = centroids(mIdx, 1:2)'; % 实际的 [X, Y]'
                        % 计算欧氏距离
                        cost_matrix(tIdx, mIdx) = norm(pred_pos - meas_pos);
                    end
                end
                
                % 贪心匹配算法 (谁离得近先抱走谁)
                while true
                    [min_val, min_idx] = min(cost_matrix(:));
                    if min_val > GATE_THRESH || isinf(min_val) % isinf(min_val)：应对GATE_THRESH设置为inf的情况
                        break; % 最小的距离都大于门限，停止关联
                    end
                    
                    % 解算出对应的航迹和测量点索引，一维索引向二维索引坐标转换
                    % matlab是列优先的存储方式  size(cost_matrix)：行数
                    [t_match, m_match] = ind2sub(size(cost_matrix), min_idx); % [余数-代表航迹，整数-列代表测量点数] 和矩阵cost_matrix行列含义对应上
                    
                    assignment(t_match) = m_match; % 航迹 t_match 成功牵手 第m_match 个观测聚类质心
                    meas_used(m_match) = 1;        % 测量点被占用了
                    
                    % 将已匹配的行和列设为无穷大，防止重复匹配
                    cost_matrix(t_match, :) = inf;
                    cost_matrix(:, m_match) = inf;
                end
            end
    
            % 3. 卡尔曼【更新】与航迹生命周期【状态机】更新
            updated_tracks_idx = [];
            for tIdx = 1:numTracks
                m_idx = assignment(tIdx);
                
                if m_idx > 0 % 测量点与航迹匹配成功，m_idx就不为0
                    % --- 情况 A: 航迹匹配成功，执行卡尔曼量化更新 (Update) ---
                    z = centroids(m_idx, 1:2)'; % 提取当前的 [X, Y] 作为观测输入
                    
                    % 卡尔曼增益计算
                    S = H * tracks(tIdx).P * H' + R;
                    K = tracks(tIdx).P * H' / S; % 矩阵右除：B/A = B*A^-1  矩阵左除：A\B = A^-1*B
                    
                    % 更新状态和协方差
                    tracks(tIdx).x = tracks(tIdx).x + K * (z - H * tracks(tIdx).x);
                    tracks(tIdx).P = (eye(4) - K * H) * tracks(tIdx).P;
                    
                    % 更新生命周期参数
                    tracks(tIdx).hit = tracks(tIdx).hit + 1;
                    tracks(tIdx).miss = 0; % 漏检清零
                    
                    % 满足 M/N 条件，转为正式激活航迹
                    if strcmp(tracks(tIdx).status, 'tentative') && tracks(tIdx).hit >= INIT_HIT
                        tracks(tIdx).status = 'active';
                    end
                else
                    % --- 情况 B: 航迹出现漏检 (脱靶) ---
                    tracks(tIdx).miss = tracks(tIdx).miss + 1;
                    tracks(tIdx).hit = 0; % 连续命中清零
                end
                
                % 将未死亡的航迹保存下来
                if tracks(tIdx).miss < MAX_MISS
                    % 只要还活着，就把当前这一帧更新后的状态记录进轨迹历史
                    tracks(tIdx).history = [tracks(tIdx).history, tracks(tIdx).x];
                    updated_tracks_idx = [updated_tracks_idx, tIdx]; % 记录存活的航迹id
                end
            end
            tracks = tracks(updated_tracks_idx); % 仅保留存活的航迹id，剔除死航迹

        else % 飞机
            % 2.【数据关联与同源点凝聚阶段】
            % 建立结构体，存储每个存活航迹关联到的【多个测量点坐标】
            track_measurements = cell(numTracks, 1);
            meas_used = zeros(numMeas, 1);    
            
            if numTracks > 0 && numMeas > 0 
                % 遍历所有存活航迹，利用波门（GATE_THRESH）圈定属于自己的多个质心点
                for tIdx = 1:numTracks
                    pred_pos = H * tracks(tIdx).x; % 当前航迹预测的 [X, Y]'
                    matched_points_idx = [];       % 存储落入波门的测量点索引
                    
                    for mIdx = 1:numMeas
                        % 只有还没被其他航迹占用的测量点，才允许参与当前航迹的关联
                        if meas_used(mIdx) == 0
                            meas_pos = centroids(mIdx, 1:2)'; 
                            dist = norm(pred_pos - meas_pos);
                            % 只要质心落入当前航迹的宏观大波门内，归于该航迹
                            if dist <= GATE_THRESH
                                matched_points_idx = [matched_points_idx, mIdx];
                            end
                        end
                    end
                    
                    % 如果找到了落入波门的同源质心（可能是机头、机身、机尾等多个点）
                    if ~isempty(matched_points_idx)
                        % 提取这些点的 [X, Y] 坐标
                        points_coord = centroids(matched_points_idx, 1:2); % [num,2]
                        points_velo  = centroids(matched_points_idx, 4);
                        
                        % 将属于同一架飞机的多个点进行【规范融合】，作为该航迹唯一的量测输入
                        z_pos = mean(points_coord, 1)'; % [num,2]对num个点求平均
                        z_v_measured = mean(points_velo);

%                         % 【真实速度超过雷达最大不模糊速度边界时发生折叠，速度解模糊，利用"余数差"来恢复真实速度】
%                         V_max = Parameter.lambda / (4 * Parameter.Tc); % 不模糊边界 V_max
%                         V_ambiguous_span = 2 * V_max; % 一个完整的折叠周期
%                         pred_vx = tracks(tIdx).x(3); % vx
%                         pred_vy = tracks(tIdx).x(4); % vy
%                         % 计算航迹当前的预测径向速度 (将直角坐标速度投影到目标与雷达连线方向上)
%                         pred_v_radial = (z_pos(1)*pred_vx + z_pos(2)*pred_vy) / norm(z_pos);
%                         % 寻找最接近卡尔曼预测速度的折叠倍数 k
%                         % 数学原理：真实速度应该在 预测速度 附近，通过四舍五入解出 k
%                         k = round((pred_v_radial - z_v_measured) / V_ambiguous_span);
%                         % 终极解模糊：恢复出真实的物理径向速度
%                         z_v_real = z_v_measured + k * V_ambiguous_span;
                        
%                         z_v_real = z_v_measured;
                        track_measurements{tIdx} = z_pos;
%                         unit_vector = z_pos / norm(z_pos);
%                         tracks(tIdx).x(3) = z_v_real * unit_vector(1); % 更新 Vx
%                         tracks(tIdx).x(4) = z_v_real * unit_vector(2); % 更新 Vy

                        % 标记这些点已被吸收，防止后面的新航迹孵化重复使用它们
                        meas_used(matched_points_idx) = 1; 
                    else
                        track_measurements{tIdx} = []; % 没抓到点，说明该航迹本帧脱靶
                    end
                end
            end
            
            % 3. 卡尔曼【更新】与航迹生命周期【状态机】更新
            updated_tracks_idx = [];
            for tIdx = 1:numTracks
                z = track_measurements{tIdx};
                
                if ~isempty(z) 
                    % --- 情况 A: 凝聚出了合法的测量点，执行卡尔曼量化更新 (Update) ---
                    S = H * tracks(tIdx).P * H' + R;
                    K = tracks(tIdx).P * H' / S; 
                    tracks(tIdx).x = tracks(tIdx).x + K * (z - H * tracks(tIdx).x);
                    tracks(tIdx).P = (eye(4) - K * H) * tracks(tIdx).P;
                    tracks(tIdx).hit = tracks(tIdx).hit + 1;
                    tracks(tIdx).miss = 0; 
                    if strcmp(tracks(tIdx).status, 'tentative') && tracks(tIdx).hit >= INIT_HIT
                        tracks(tIdx).status = 'active';
                    end
                else
                    % --- 情况 B: 整个波门内没有点（脱靶），由滤波器盲操预测维持航迹 ---
                    tracks(tIdx).miss = tracks(tIdx).miss + 1;
                    tracks(tIdx).hit = 0; 
                end
                
                % 将未死亡的航迹保存下来
                if tracks(tIdx).miss < MAX_MISS
                    tracks(tIdx).history = [tracks(tIdx).history, tracks(tIdx).x];
                    updated_tracks_idx = [updated_tracks_idx, tIdx]; 
                end
            end
            tracks = tracks(updated_tracks_idx);
        end
        
        
        % 4. 【孵化新航迹】对没有被航迹关联上的"野点迹"【测量点】初始化为新暂态目标
        for mIdx = 1:numMeas
            if meas_used(mIdx) == 0
                new_track = struct();
                new_track.id = next_id;
                
                % 初始化运动状态 [X, Y, Vx, Vy]'
                % 多普勒径向速度 centroids(mIdx, 4) 是个强先验信息
                % 针对 2D 跟踪，我们可以直接将其分解或暂时初始化为 0（卡尔曼很快能迭代出来）
                init_X = centroids(mIdx, 1);
                init_Y = centroids(mIdx, 2);
                new_track.x = [init_X; init_Y; 0; 0]; 
                
                new_track.P = P_init;
                new_track.status = 'tentative'; % 新生目标一律先处于考察暂态
                new_track.hit = 1;
                new_track.miss = 0;
                new_track.history = new_track.x; % 记入首个点
                tracks(end + 1) = new_track;
                next_id = next_id + 1;
            end
        end
        
        % 指标1：航迹存活率：雷达系统对于视场内"上帝真值目标”的持续霸气锁死能力。它反映的是算法“跟得稳不稳、会不会无故跟丢”。
        % 指标2：物理残差RMSE：没跟丢的情况下，位置和速度【角度】测得准不准
        
        if has_truth
            truth_frame = all_frames_truth{frameId};
            if ~isempty(truth_frame)
                % 提取当前帧真值中所有不重复的真实目标 ID
                unique_target_ids = unique(truth_frame(:, 4)); 
                K_true = length(unique_target_ids); % 利用 TargetID 动态融合宏观真实目标
                
                % 初始化宏观真值飞机的直角坐标阵与物理值
                X_true = zeros(K_true, 1);
                Y_true = zeros(K_true, 1);
                V_true = zeros(K_true, 1);
                
                % 遍历每个目标，将其所有碎点的物理信息融合成一个宏观质心
                for tidx = 1:K_true
                    current_id = unique_target_ids(tidx);
                    id_mask = (truth_frame(:, 4) == current_id); % 找出属于当前 ID 的所有行
                    
                    % 提取当前目标的极坐标物理量
                    r_sub = truth_frame(id_mask, 1);
                    v_sub = truth_frame(id_mask, 2);
                    a_sub = truth_frame(id_mask, 3);
                    
                    % 单点转换到直角坐标系
                    x_sub = r_sub .* sind(a_sub);
                    y_sub = r_sub .* cosd(a_sub);
                    
                    % 取平均值，得到这个目标的真实宏观几何中心和物理量
                    X_true(tidx) = mean(x_sub);
                    Y_true(tidx) = mean(y_sub);
                    V_true(tidx) = mean(v_sub); 
                end
            else
                K_true = 0;
                X_true = []; Y_true = []; V_true = [];
            end
            
            tracks = tracks(:); 
            active_mask = false(length(tracks), 1);
            for kk = 1:length(tracks)
                
                if strcmp(tracks(kk).status, 'active')
                    active_mask(kk) = true;
                end
            end
            active_tracks = tracks(active_mask);
            M_active = length(active_tracks);
            
            % --- 1. 计算航迹存活率 (Track Survival Rate) ---
            % 航迹存活率 
            % 此时分子分母量纲绝对对齐：飞机【架数】 / 飞机【架数】
            if K_true > 0
                track_survival_rates(frameId) = min(1.0, M_active / K_true);
            else
                track_survival_rates(frameId) = 1.0;
            end
            
            % --- 准备卡尔曼估计值的直角坐标转换 ---
            X_est = zeros(M_active, 1); Y_est = zeros(M_active, 1); 
            V_est = zeros(M_active, 1); 
            for t = 1:M_active
                X_est(t) = active_tracks(t).x(1);
                Y_est(t) = active_tracks(t).x(2);
%                 V_est(t) = (X_est(t)*active_tracks(t).x(3) + Y_est(t)*active_tracks(t).x(4)) / max(0.1, norm([X_est(t), Y_est(t)]));
            end
            
            % --- 2. 各物理状态维统计 MSE 收集 ---
            if M_active > 0 && K_true > 0
                for t = 1:M_active
                    dist_to_truths = sqrt((X_est(t) - X_true).^2 + (Y_est(t) - Y_true).^2);
                    [min_d, matched_true_idx] = min(dist_to_truths);
                    
                    if min_d <= EVAL_GATE 
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
                        frame_vel_se   = [frame_vel_se;   (v_est_radial - v_true_radial)^2];
                    end
                end
            end
            
            % --- 3. 终极多目标考核：OSPA 距离解算 ---
            % 用融合后的点集传给 OSPA，彻底消除虚假的基数惩罚
            Set_True = [X_true, Y_true];
            Set_Est  = [X_est,  Y_est];
            ospa_distances(frameId) = calculate_ospa(Set_True, Set_Est, c_ospa, p_ospa);
        end

    end

    SummaryMetrics.Range_RMSE = sqrt(mean(frame_range_se));
    SummaryMetrics.Vel_RMSE   = sqrt(mean(frame_vel_se));
    SummaryMetrics.Angle_RMSE = sqrt(mean(frame_ang_se));
    SummaryMetrics.Mean_Track_Survival_Rate = mean(track_survival_rates(INIT_HIT:end)); % 从INIT_HIT帧才开始有存活航迹
    SummaryMetrics.Mean_OSPA  = mean(ospa_distances);
    
    %% --- 输出封装：只向用户输出真正激活过的、高质量的合法目标航迹 ---
    valid_count = 0;
    tracked_objects = struct([]);
    for tIdx = 1:length(tracks)
        % 过滤掉那些刚冒出来就死掉的噪点航迹，必须曾经是 'active'
        if strcmp(tracks(tIdx).status, 'active') || tracks(tIdx).hit >= INIT_HIT
            valid_count = valid_count + 1;
            tracked_objects(valid_count).id = tracks(tIdx).id;
            tracked_objects(valid_count).trajectory = tracks(tIdx).history; % [4 x 存活帧数] 矩阵
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
                % 提取最后一帧的滤波估计状态量 [X; Y; Vx; Vy]
                last_est = traj(:, end); 
                tx = last_est(1);  ty = last_est(2); 
                tvx = last_est(3); tvy = last_est(4);
                
                % 1. 计算绝对真实总航速 (二维合成标量，恒为正)
                total_speed = norm([tvx, tvy]);
                
                % 2. 还原带物理符号的径向运动速度 (正代表远离，负代表靠近雷达)
                r_est = max(0.1, sqrt(tx^2 + ty^2));
                radial_speed_with_sign = (tx*tvx + ty*tvy) / r_est;
                
                % 3. 规范化打印输出
                fprintf(' 目标 ID [%d] 最后一帧-> \n', tracked_objects(i).id);
                fprintf('   [绝对总航速(大小)] : %.4f m/s\n', total_speed);
                fprintf('   [径向航速(带符号)] : %.4f m/s\n', radial_speed_with_sign);
                fprintf('   [直角分量状态解算] : Vx = %.4f m/s, Vy = %.4f m/s\n', tvx, tvy);
                fprintf(' ----------------------------------------------------------------\n');
            end
        end
        fprintf('==================================================================\n\n');
    end
    % =======================================================================    
end

% 指标3：OSPA（Optimal Subpattern Assignment，最优子模式分配）：既考核距离误差，又严惩虚警和漏检（势误差）
% c_ospa：截止惩罚常数，当算法漏检或者虚警时，每多一个或少一个目标，系统应该给予的硬性距离惩罚，
% 也是位置误差的上限，如果两个目标距离超过 c_ospa，其位置误差也直接按 c_ospa 计算，
% 若设置 c_ospa = 10 米，一旦雷达跟丢了一个目标，OSPA 会直接在分母里记上一笔 10 米的惩罚
% p_ospa：ospa阶数，决定指标对大误差的敏感程度，p_ospa=2时，系统对偶尔出现的大跳变、大误差惩罚极重（类似于均方误差），能更灵敏地反映系统的极端不稳定状态。
% 平均ospa值越低，说明卡尔曼滤波器不仅位置估计准确，而且对目标生命周期（漏检，虚警）管理得越完美
function ospa_dist = calculate_ospa(X, Y, c, p)
    M = size(X, 1); % X: 真实值 [M x 2] (直角坐标系)
    N = size(Y, 1); % Y: 估计值 [N x 2] (直角坐标系)
    
    % 处理极端无目标场景
    if M == 0 && N == 0
        ospa_dist = 0; return;
    elseif M == 0 || N == 0
        ospa_dist = c; return;
    end
    
    % 1. 大小配对检查，永远保证把较短的集合放在左边。
    % 如果估计的目标数 n 比真值数 m 还要多（出现了虚警），公式在数学上就会把 X 和 Y 互换。
    % 确保 X 为较短或对等的集合，方便执行匈牙利或最邻近排列
    swapped = false;
    if M > N % 以 m < n（即真值数 < 估计数）为例
        temp = X; X = Y; Y = temp;
        M = size(X, 1); N = size(Y, 1);
        swapped = true;
    end
    
    % 2. 计算短集合X中的每一个点到长集合Y中每一个点的欧式距离，并用门限c阻断拦截
    % 计算两个集合两两之间的距离矩阵
    dist_matrix = zeros(M, N);
    for i = 1:M
        for j = 1:N
            dist_matrix(i, j) = min(c, norm(X(i, :) - Y(j, :)));
        end
    end
    
    % 3. 最优分配。算法从长集合Y中，挑选m个点和短集合X的m个点进行一对一的最近邻配对使得配对后的总位置误差之和最小。
    % 这一步通常由匈牙利算法、KM 算法或者贪心最近邻算法高效完成。
    % 采用贪心策略寻找最优子模式分配（Subpattern Assignment）：先到先得、到手拉黑
    total_loc_error_p = 0;
    for i = 1:M
        [~, min_idx] = min(dist_matrix(i, :)); % 让每个目标挑选距离它最近的点
        total_loc_error_p = total_loc_error_p + dist_matrix(i, min_idx)^p; % 位置误差项
        dist_matrix(:, min_idx) = inf; % 防止重复配对
    end
    
    % 4. 找到最优配对后，代入OSPA公式
    % OSPA 终极数学公式：位置误差项 + 势口基数惩罚项
    cardinality_penalty = (c^p) * (N - M); % 势误差（惩罚虚警和漏检），对于多出来的目标，给予最大值惩罚
    ospa_dist = ((total_loc_error_p + cardinality_penalty) / N)^(1/p);

end



