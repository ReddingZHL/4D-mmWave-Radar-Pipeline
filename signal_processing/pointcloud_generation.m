function pointCloud = pointcloud_generation(detections)
    % detections: CFAR检测后的 [range, doppler, angle, rcs]
    detectPointsNum = size(detections, 1);
    pointCloud = zeros(detectPointsNum, 5); % [x y z v power]
    for i = 1:detectPointsNum
        range = detections(i,1); % 距离估计
        velocity = detections(i,2); % 速度估计
        az    = detections(i,3); % 角度估计
        
        % 转笛卡尔坐标
        x = range * sind(az);
        y = range * cosd(az);
        z = 0;   % 当前假设2D线性阵列，垂直高度设为0，【可扩展为俯仰】
        
        pointCloud(i, :) = [x, y, z, velocity, detections(i,4)];  % [x y z 径向v power]
    end
    
    % 绘制点云图
    X_pc = pointCloud(:, 1);     % 纵向距离 (m)
    Y_pc = pointCloud(:, 2);     % 横向距离 (m)
    Z_pc = pointCloud(:, 3);     % 高度 (m)
    V_pc = pointCloud(:, 4);     % 径向速度 (m/s)
    P_pc = pointCloud(:, 5);     % 反射功率 (dB或线性)
    % 2. 功率映射为点的大小 (防止功率数值太大或太小导致点消失或撑满屏幕)
    % 使用归一化动态映射，将点的大小限制在 20 到 150 之间
    if max(P_pc) ~= min(P_pc)
        marker_size = 20 + 130 * (P_pc - min(P_pc)) / (max(P_pc) - min(P_pc));
    else
        marker_size = 40 * ones(size(P_pc)); % 若功率全等，固定大小
    end
%     % 绘制点云图
%     figure('Color', [1 1 1]);
%     % 使用 scatter3 绘制4D点云
%     % X, Y, Z 为空间坐标；marker_size 决定大小；V_pc 决定颜色
%     scatter3(X_pc,Y_pc,  Z_pc, marker_size, V_pc, 'filled', 'MarkerEdgeColor', [0.3 0.3 0.3]);
%     % 4. 绘制雷达原点，方便对比相对位置
%     hold on;
%     plot3(0, 0, 0, 'rp', 'MarkerSize', 15, 'MarkerFaceColor', 'r');
%     text(0, 0, 0.5, ' 雷达(0,0,0)', 'Color', 'r', 'FontWeight', 'bold');
%     % 5. 润色图表
%     grid on; box on;
%     colormap(jet); % 蓝冷色代表速度小/靠近，红暖色代表速度大/远离
%     cb = colorbar;
%     ylabel(cb, '径向速度 v (m/s)', 'FontSize', 11);
%     % 坐标轴标签 (严格对照雷达直角坐标系)
%     xlabel('纵向距离 X (米)', 'FontSize', 11);
%     ylabel('横向距离 Y (米)', 'FontSize', 11);
%     zlabel('高度 Z (米)', 'FontSize', 11);
%     title('雷达 4D 空间点云成像结果 (大小:功率, 颜色:速度)', 'FontSize', 12, 'FontWeight', 'bold');
%     % 保持几何比例1:1，防止飞机被压扁或者拉长
%     axis equal;view(0,90);
    
end