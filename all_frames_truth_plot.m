%% =======================================================================
%% 多目标 质心运动轨迹动态可视化
%% =======================================================================
clear; clc; close all;
matFilePath = 'data/radar_pipeline_inputs.mat'; % radar_pipeline_inputs_2thTarget_10frames
if ~exist(matFilePath, 'file')
    error('未找到数据文件：%s，请先运行回波生成仿真！', matFilePath);
end
load(matFilePath, 'Parameter', 'all_frames_truth'); % 先加载轻量级雷达配置及目标先验信息
numFrames = Parameter.numFrames;

fprintf('========= 正在绘制目标质心实时运动轨迹... =========\n');
figure('Name', '雷达回波生成：目标质心运动轨迹对账控制台', 'Position', [100, 100, 850, 700]);
hold on; grid on; box on;

%% 1. 自动探查全场历史中出现过的最大目标总数，动态分配存储空间
max_targets = 0;
for frameId = 1:numFrames
    if ~isempty(all_frames_truth{frameId})
        max_targets = max(max_targets, max(all_frames_truth{frameId}(:, 4)));
    end
end

% 将轨迹容器升级为多列矩阵（每一列负责一个目标的独立航迹）
X_traj = nan(numFrames, max_targets);
Y_traj = nan(numFrames, max_targets);
V_radial_traj = nan(numFrames, max_targets);

colors = jet(numFrames); % 使用渐变色表示时间轴演进（蓝 -> 红）

%% 2. 双重循环解析：先抽取数据到矩阵中
for frameId=1:numFrames
    targets = all_frames_truth{frameId};
    if isempty(targets), continue; end

    for target_id = 1:max_targets
        id_mask = (targets(:, 4) == target_id);
        if any(id_mask)
            r_sub = targets(id_mask, 1);
            v_sub = targets(id_mask, 2);
            a_sub = targets(id_mask, 3);
            
            % 精准解算直角坐标中心
            X_traj(frameId, target_id) = mean(r_sub .* sind(a_sub));
            Y_traj(frameId, target_id) = mean(r_sub .* cosd(a_sub));
            V_radial_traj(frameId, target_id) = mean(v_sub);
        end
    end
end

%% 3. 目标质心轨迹绘制
for target_id = 1:max_targets
    % 剔除可能存在的无效 NaN 帧（比如有些目标在中间帧才出生或消亡）
    valid_mask = ~isnan(X_traj(:, target_id)) & ~isnan(Y_traj(:, target_id));
    if ~any(valid_mask), continue; end
    
    X_curr = X_traj(valid_mask, target_id);
    Y_curr = Y_traj(valid_mask, target_id);
    V_curr = V_radial_traj(valid_mask, target_id);
    valid_indices = find(valid_mask);
    
    % A. 绘制当前目标随时间渐变的散点
    for idx = 1:length(X_curr)
        f_id = valid_indices(idx);
        scatter(X_curr(idx), Y_curr(idx), 45, colors(f_id, :), 'filled', 'MarkerFaceAlpha', 0.8, 'HandleVisibility', 'off');
    end
    
    % B. 绘制该目标的专属轨迹连线
    plot(X_curr, Y_curr, 'k--', 'LineWidth', 1.2, 'HandleVisibility', 'off');
    
    % C. 精准标记各个目标的独立起点和终点
    plot(X_curr(1), Y_curr(1), 'go', 'MarkerSize', 5, 'LineWidth', 2, ...
        'DisplayName', sprintf('目标 %d 起点', target_id));
    plot(X_curr(end), Y_curr(end), 'rp', 'MarkerSize', 6, 'LineWidth', 2, ...
        'DisplayName', sprintf('目标 %d 终点', target_id));
    
    % D. 动态绘制物理矢量航向箭头（取该目标生存周期的中间帧）
    if length(X_curr) > 1
        mid_idx = round(length(X_curr) / 2);
        vx_dir = X_curr(mid_idx+1) - X_curr(mid_idx);
        vy_dir = Y_curr(mid_idx+1) - Y_curr(mid_idx);
        
        dir_scale = sqrt(vx_dir^2 + vy_dir^2);
        if dir_scale > 1e-4
            % 绘制局部缩放的绿色航向箭头
            quiver(X_curr(mid_idx), Y_curr(mid_idx), (vx_dir/dir_scale)*3, (vy_dir/dir_scale)*3, ...
                'MaxHeadSize', 2, 'Color', [0.1 0.6 0.1], 'LineWidth', 2, ...
                'DisplayName', sprintf('目标 %d 航向', target_id));
        end
    end
    
    % E. 物理账目终端打印
    fprintf('🎯 目标 [%d] 物理账目解析:\n', target_id);
    fprintf('   -> 轨迹运动区间: (X:%.1f, Y:%.1f) 走向 (X:%.1f, Y:%.1f) 米\n', ...
        X_curr(1), Y_curr(1), X_curr(end), Y_curr(end));
    fprintf('   -> 最后一帧多普勒注入速度 (真值): %.4f m/s\n', V_curr(end));
end

%% 4. 终极图表修饰
xlabel('X 笛卡尔坐标 (东-米)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Y 笛卡尔坐标 (北-米)', 'FontSize', 11, 'FontWeight', 'bold');
title('🎯 雷达前级回波：多目标质心真实运动轨迹 (雷达位于原点[0,0])', 'FontSize', 13, 'FontWeight', 'bold');

colormap(jet);
cb = colorbar;
ylabel(cb, '时间轴演进 (帧数颜色映射：蓝 -> 红)', 'FontSize', 10, 'FontWeight', 'bold');

axis equal;
legend('Location', 'best');
set(gca, 'FontSize', 10, 'LineWidth', 1.2);
hold off;

fprintf('========= 轨迹绘制完成！全场目标已隔离，请查看可视化看板。 =========\n');