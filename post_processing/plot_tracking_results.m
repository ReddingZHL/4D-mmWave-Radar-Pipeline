function plot_tracking_results(tracked_objects, Parameter)
    % =====================================================================
    % 函数功能：可视化输出多目标跟踪轨迹，并标注/打印各个目标的起始点与终止点
    % =====================================================================
    if isempty(tracked_objects)
        warning('没有捕获到任何有效的目标激活航迹！');
        return;
    end

    dt = Parameter.T_frame;
    airplanes_config = Parameter.airplanes_config; % 真实目标质心的[距离, 角度, 速度, 航向角]
    track_to_config_idx = zeros(length(tracked_objects), 1);
    for idx = 1:length(tracked_objects)
        traj = tracked_objects(idx).trajectory;
        x_start = traj(1, 1); y_start = traj(2, 1);
        r_estimated_start = sqrt(x_start^2 + y_start^2);
        
        % 计算该航迹起点与仿真配置中哪架飞机的设定距离最接近
        [~, best_match] = min(abs(airplanes_config(:, 1) - r_estimated_start));
        track_to_config_idx(idx) = best_match;
    end
    
    % 创建高画质画布
    figure('Color', [1 1 1], 'Name', '2D 目标跟踪轨迹'); % [1 1 1]表示将画布背景设成纯白色
    hold on; grid on; box on; % box on表示显示坐标轴外框
    
    % 1. 绘制雷达原点（参考基准）
    plot(0, 0, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r'); % 默认'ro'是空心圆圈，加上'MarkerFaceColor', 'r'将变成实心的红色圆圈
    text(0, -2, '雷达中心', 'HorizontalAlignment', 'center', 'FontWeight', 'bold'); % 添加文本标签，设置文字在水平方向上居中对齐，'FontWeight', 'bold'设字体为加粗
    
    % 动态分配颜色（确保每个目标轨迹颜色不同）
    colors = lines(length(tracked_objects)); % lines是颜色图生成器，lines(N)会生成一个N*3的RGB矩阵
    
    fprintf('\n==================== 目标追踪物理轨迹 ====================\n');
    for idx = 1:length(tracked_objects)
        traj = tracked_objects(idx).trajectory; % 提取当前目标的卡尔曼状态历史 [4 x N]
        X_history = traj(1, :); Y_history = traj(2, :);
        
        % 提取物理首尾坐标
        x_start = X_history(1);   y_start = Y_history(1);
        x_end   = X_history(end); y_end   = Y_history(end);
        total_displacement = sqrt((x_end - x_start)^2 + (y_end - y_start)^2); % 计算总位移

        % 宏观几何反推航向角 (对齐Y轴0度顺时针定义)
        macro_heading = mod(atan2d(x_end - x_start, y_end - y_start), 360); % 确保在 0~360 度之间
        
        % 2. 利用卡尔曼平滑状态（用所有帧的vx和vy合成速度的平均值）解算平均航速
        vx_history = traj(3, :); vy_history = traj(4, :);
        % --- 计算最近 10 帧 【收敛了的平稳期】 的滑动平均航速 ---
        window_size = min(10, length(vx_history));
        % 通过位置与速度点积判定靠近/远离
        vx_win = vx_history(end-window_size+1:end);
        vy_win = vy_history(end-window_size+1:end);
        x_win  = X_history(end-window_size+1:end);
        y_win  = Y_history(end-window_size+1:end);
        speeds_win_mag = sqrt(vx_win.^2 + vy_win.^2);
        dot_product_win = x_win.*vx_win + y_win.*vy_win; % =位置矢量模 * 速度矢量模 * cos(两矢量夹角)，角度为钝角，结果是负数，表示靠近雷达
        stable_speed = mean(sign(dot_product_win) .* speeds_win_mag);

        
        % 提取配对的上帝视角真值
        cfg_id = track_to_config_idx(idx);
        true_speed = airplanes_config(cfg_id, 3);
        true_heading = mod(airplanes_config(cfg_id, 4), 360);
        
        speed_error_pct = (stable_speed - true_speed) / abs(true_speed) * 100;
        heading_error_pct = (macro_heading - true_heading) / true_heading * 100;
        
        % --- 终端打印报告 ---
        fprintf('航迹 ID [%d] -> 成功匹配仿真【目标 %d】:\n', tracked_objects(idx).id, cfg_id);
        fprintf('  -> 起止点区间   : (X:%.1f, Y:%.1f) 走向 (X:%.1f, Y:%.1f) 米\n', x_start, y_start, x_end, y_end);
        fprintf('  -> 宏观总位移   : %.2f 米 (共存活 %d 帧)\n', total_displacement, length(X_history));
        fprintf('  -> 平均航速 : 估计 = %+.4f m/s | 真值 = %+.2f m/s | 误差 = %+.4f m/s (相对误差: %+.2f%%)\n', ...
                stable_speed, true_speed, (stable_speed - true_speed), speed_error_pct);
        fprintf('  -> 航向角 : 估计 = %.2f °   | 真值 = %.2f °   | 误差 = %+.2f °   (相对误差: %+.2f%%)\n', ...
                macro_heading, true_heading, (macro_heading - true_heading), heading_error_pct);
        fprintf('------------------------------------------------------------\n');

%         text(x_end + 0.6, y_end - 0.4, sprintf('v=%.2fm/s, H=%.1f°', stable_speed, macro_heading), ...
%              'Color', [0.3 0.3 0.3], 'FontSize', 9);
        
        % --- 2. 绘制连续轨迹线 ---
        plot(X_history, Y_history, '-', 'LineWidth', 2.5, 'Color', colors(idx, :)); % 将轨迹串联成线
        
        % --- 3. 2D 图形精准标注 ---
        % 起始点：绿色正方形 (Square)
        plot(x_start, y_start, 'gs', 'MarkerSize', 8, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'k'); % 填充实心，边框黑色
        % 终止点：红色五角星 (Pentagram)
        plot(x_end, y_end, 'kp', 'MarkerSize', 11, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k');
        
        % --- 4. 图上文字标签 ---
        % 在起始点旁标注 "S" (Start)
        text(x_start + 0.6, y_start, sprintf('S%d', tracked_objects(idx).id), ...
             'Color', [0 0.6 0], 'FontWeight', 'bold', 'FontSize', 10);
        % 在终止点旁标注 "Target ID [E]" (End)
        text(x_end + 0.6, y_end, sprintf('Target %d (E)', tracked_objects(idx).id), ...
             'Color', 'r', 'FontWeight', 'bold', 'FontSize', 10);
    end
    
    fprintf('============================================================\n\n');
    
    % --- 5. 图例与美化 ---
    % 创建虚拟图例对象，防止图例把每一帧的点都抓进去导致混乱
    h_start = plot(nan, nan, 'gs', 'MarkerSize', 8, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'k');
    h_end   = plot(nan, nan, 'kp', 'MarkerSize', 11, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k');
    legend([h_start, h_end], {'起始点 (Start)', '终止点 (End)'}, 'Location', 'best'); % 'Location', 'best' 会自动寻找最不挡线的位置
    
    xlabel('X 方向距离 (米)');
    ylabel('Y 方向距离 (米)');
    title('FMCW 雷达 2D 多目标跟踪轨迹');
    axis equal; % 将X轴和Y轴的刻度比例设为完全相等

    % =====================================================================
    % 画布 2：绘制带物理符号的时域航速与航向误差演进对比图
    % =====================================================================
    figure('Color', [1 1 1], 'Name', '卡尔曼滤波器性能误差分析');
    
    % ---- 子图 1：带符号航速误差随时间演进曲线 ----
    subplot(2, 1, 1); hold on; grid on; box on;
    for idx = 1:length(tracked_objects)
        traj = tracked_objects(idx).trajectory;
        X_hist = traj(1, :);  Y_hist = traj(2, :);
        vx_hist = traj(3, :); vy_hist = traj(4, :);
        
        % 1. 计算每一帧的无符号速率模长
        speeds_mag = sqrt(vx_hist.^2 + vy_hist.^2);
        % 2. 核心：通过直角坐标点积判定该帧卡尔曼是在靠近（负）还是远离（正）
        dot_product = X_hist.*vx_hist + Y_hist.*vy_hist;
        estimated_speeds_signed = sign(dot_product) .* speeds_mag; % 完美注入符号！
        
        cfg_id = track_to_config_idx(idx);
        true_speed_signed = airplanes_config(cfg_id, 3); % 保留原始设定负速度
        
        % 计算带符号的速度误差 (估计值 - 真值)
        speed_errors = estimated_speeds_signed - true_speed_signed;
        plot(speed_errors, '-o', 'LineWidth', 2, 'Color', colors(idx, :), ...
             'MarkerSize', 4, 'DisplayName', sprintf('目标 %d (真值:%+.1fm/s)', tracked_objects(idx).id, true_speed_signed));
    end
    yline(0, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off'); 
    xlabel('跟踪帧数 (Frame)'); ylabel('航速误差 (m/s)');
    title('各目标航速估计误差时域收敛曲线 (包含靠近为负、远离为正)');
    legend('Location', 'best');
    
    % ---- 子图 2：航向角误差随时间演进曲线 ----
    subplot(2, 1, 2); hold on; grid on; box on;
    for idx = 1:length(tracked_objects)
        traj = tracked_objects(idx).trajectory;
        vx_hist = traj(3, :); vy_hist = traj(4, :);
        
        estimated_headings = mod(atan2d(vx_hist, vy_hist), 360);
        cfg_id = track_to_config_idx(idx);
        true_heading = mod(airplanes_config(cfg_id, 4), 360);
        
        heading_errors = estimated_headings - true_heading;
        heading_errors = mod(heading_errors + 180, 360) - 180; % 拉平跨0度跳变
        
        plot(heading_errors, '-^', 'LineWidth', 2, 'Color', colors(idx, :), ...
             'MarkerSize', 4, 'DisplayName', sprintf('目标 %d (真值:%.1f°)', tracked_objects(idx).id, true_heading));
    end
    yline(0, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    xlabel('跟踪帧数 (Frame)'); ylabel('航向角误差 (度)');
    title('各目标航向角估计误差时域收敛曲线');
    legend('Location', 'best'); 
    
    % --- 6. 自动保存图片逻辑 ---
    if ~exist('results', 'dir')
        mkdir('results');
    end
    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    saveas(figure(1), sprintf('results/tracking_trajectory_%s.png', timestamp)); 
    saveas(figure(2), sprintf('results/tracking_errors_signed_%s.png', timestamp)); 
    fprintf('\n【系统提示】轨迹图与【带符号】误差分析图已保存至 results/ 文件夹下。\n');
end
