%% 雷达参数设置
function parameter = generateParameter(targetnum)
    % === 定义物理常数与雷达硬件指标 ===
    parameter.kB = 1.38e-23;       % 玻尔兹曼常数
    parameter.T = 290;             % 环境温度 (K)
    parameter.Pt = 0.0158;         % 发射功率 12 dBm, 10^((12 - 30) / 10)≈0.01585 W
    parameter.G_ant = 10^(15/10);  % 天线增益 15 dBi
    parameter.Loss = 10^(4/10);    % 系统损耗 4 dB
    parameter.NoiseFigure = 10^(10/10); % 接收机噪声系数 10 dB

    parameter.c = 3e8;                                  %光速
    parameter.stratFreq = 77e9;                         %起始频率
    parameter.targetnum = targetnum;
    if targetnum == 1 % 1-点目标
        parameter.Tr =  40e-6;  % 10e-6; 60e-6;                               %【扫频时间】
    else % 飞机目标
        parameter.Tr =  60e-6;40e-6;
    end
    parameter.Idle_time = 10e-6; % 100e-6;                       %空闲时间
    parameter.Tc = parameter.Tr+parameter.Idle_time;    %Chirp之间的间隔
    parameter.Fs = 25.6e6; %10e6;                                %【采样率】
    parameter.Samples = parameter.Tr * parameter.Fs; % 采样点，信号生成：单个chirp采样点数 = 扫频时间*采样率
    
    parameter.rangeBin = parameter.Samples; %2^nextpow2(parameter.Samples);            %rangebin 【为了计算效率，通常将 rangeBin 设置为大于 Samples 的最小的 2 的幂（例如采样 200 点，FFT 用 256 点）】
    parameter.Chirps = 512; % 128;                             %【】chirp数，在一帧数据内，雷达连续发射的chirp总数
    parameter.dopplerBin = parameter.Chirps;            %dopplerbin

    parameter.Slope = 29.982e12;                            %chirp斜率
    parameter.Bandwidth = parameter.Slope * parameter.Tr ;  %发射信号带宽
    parameter.BandwidthValid = parameter.Samples/parameter.Fs*parameter.Slope;  %发射信号有效带宽
    parameter.centerFreq = parameter.stratFreq + parameter.Bandwidth / 2;       %中心频率
    parameter.lambda = parameter.c / parameter.centerFreq;  %波长（中心频率波长）
    
    parameter.numFrames = 50; % 50; %10;    % 帧数
    parameter.T_frame = 80e-3;   % 帧周期

    parameter.txAntenna = ones(1,3); %发射天线个数
    parameter.rxAntenna = ones(1,4); %接收天线个数
    parameter.txNum = length(parameter.txAntenna);
    parameter.rxNum = length(parameter.rxAntenna);
    parameter.virtualAntenna = length(parameter.txAntenna) * length(parameter.rxAntenna);
    parameter.angleBin = 180;           
    
    parameter.dz = parameter.lambda / 2; %接收天线俯仰间距
    parameter.dx = parameter.lambda / 2; %接收天线水平间距
    
    parameter.doaMethod = 3; % 测角方法选择 1:fft  2:dbf  3:1Dmusic  4:RDMusic  5:Capon
    parameter.clusterMethod = 1; % 点云聚类方法选择 1-DBSCAN 2-Kmeans 3-X-Means
    
    % ===  定义目标参数 === 【注意：给不同的目标添加唯一的targetID】
    % 设置目标速度时，需要根据雷达理论可测的最大不模糊速度，保证设定速度范围在这之内
    if targetnum == 1 % 1-点目标
        raw_points = [
            100 -4   0;    %target1 range speed horizontalangle
            20   3   45;  %target2 range speed horizontalangle
            20   4   -30;   %target2 range speed horizontalangle
        ]; 
        num_pts = size(raw_points, 1);
        % 为点目标赋予 1, 2, 3... 的独立 Target_ID 编号
        parameter.target = [raw_points, (1:num_pts)'];
    else % 飞机目标
        % =================================================================
        % 构建飞机目标的参数配置矩阵 
        % 每一行代表一架飞机目标: [中心距离(m), 目标中心与雷达连线角度(°), 总速度(m/s), 航向角(°)]
        % 中心角度（雷达角度，方位角）是以 +Y 轴为 0° 基准。向右偏为正（x 变大），向左偏为负。
        % 航向角也是以 +Y 轴为 0° 基准
        % 【侧滑角是以机头为基准】
        % =================================================================
        airplanes_config = [ % [距离, 目标中心与雷达角度, 航向速度, 航向角]
            60,   10,  -4,   225;   % 飞机 1: 中距离，迎面快速斜切
            100, -20,  -4, 110;   % 飞机 2: 远距离，左侧斜穿
            40,   35,   3,  75;   % 飞机 3: 近距离，大角度脱离 
        ];
        parameter.airplanes_config = airplanes_config; % [距离, 目标中心与雷达角度, 航向恒定速度, 航向角]
        num_airplanes = size(airplanes_config, 1);
        all_airplane_targets = []; % 用于拼接所有飞机散射点的终极矩阵
        
        % 循环生成每架飞机的点阵，并融合成一个统一的雷达观测矩阵
        for airplane_idx = 1:num_airplanes
            c_range  = airplanes_config(airplane_idx, 1);
            c_angle  = airplanes_config(airplane_idx, 2);
            v_total  = airplanes_config(airplane_idx, 3);
            h_error  = airplanes_config(airplane_idx, 4);
            
            % 调用你原有的单飞机点阵生成函数
            single_airplane = generate_airplane_target(c_range, c_angle, v_total, h_error);
            num_scatterers = size(single_airplane, 1);
            id_column = ones(num_scatterers, 1) * airplane_idx;
            single_airplane_with_id = [single_airplane, id_column];
            % 垂直拼接：把所有飞机的散射点汇聚到一起
            all_airplane_targets = [all_airplane_targets; single_airplane_with_id];
        end
        
        % 将拼接后的数百个散射点赋给统一的 target 变量，无缝对接到你的雷达回波生成器
        parameter.target = all_airplane_targets;
        
%         % 为了向后兼容或绘图，保留第一个目标的简要参数（可选）
%         parameter.center_range = airplanes_config(1, 1);
%         parameter.center_angle = airplanes_config(1, 2);
%         parameter.v_airplane   = airplanes_config(1, 3);
%         parameter.heading_error = airplanes_config(1, 4);
%         
%         % === 可视化飞机初始位置，可以解开此段注释 ===
%         current_targets = parameter.target;
%         R = current_targets(:, 1); V_radial = current_targets(:, 2); Th = current_targets(:, 3);
%         X_radar_mesh = R .* cosd(Th); Y_radar_mesh = R .* sind(Th);
%         figure(11); scatter(Y_radar_mesh, X_radar_mesh, 30, V_radial, 'filled', 'MarkerEdgeColor', 'k');
%         grid on; box on; colormap(jet); colorbar; axis equal; hold on;set(gcf, 'Color', 'white');
%         plot(0, 0, 'rp', 'MarkerSize', 12, 'MarkerFaceColor', 'r');

    end
    
end


function target_matrix = generate_airplane_target(center_range, center_angle, v_total, heading_error)
    % =====================================================================
    % center_range : 飞机中心到雷达的初始距离 (m)
    % center_angle : 飞机中心处于雷达的哪个方位角 (°) [+Y轴为0°, 向右为正]
    % v_total      : 飞机的实际航行速率，【正数代表远离雷达趋势，负数代表靠近雷达趋势】
    % heading_error: 飞机的【绝对飞行航向角】(°)。以+Y轴(正北)为0°，顺时针向右飞为正。
    %【飞机产生斜切角，机身点云就会显现？？？？？？？？？？？？？？？？？？？？？】
    % =====================================================================
    
    X_local = []; Y_local = [];
    
    % --- 1. 机身边缘（纺锤体机身） ---
    x_fuse_steps = -5:0.5:5; % 飞机全长10m
    for x = x_fuse_steps
        r_fuse = 0.4 * cos(x/6); % 约束机身宽度，头尾窄，中间宽
        X_local = [X_local, x, x];
        Y_local = [Y_local, r_fuse, -r_fuse]; % 机身对称，最宽时为0.4*2=0.8m，两侧边缘最窄时0.27*2=0.54m
%         figure(11);scatter(x, r_fuse, 30, 'r', 'filled');xlim([-6, 6]);ylim([-4.5, 4.5]);
%         hold on;scatter(x, -r_fuse, 30, 'r', 'filled');xlim([-6, 6]);ylim([-4.5, 4.5]);set(gcf, 'Color', 'white');box on;axis equal;
    end
    
    % --- 2. 后掠机翼轮廓线 ---
    y_wing_steps = [-4:0.4:-0.4, 0.4:0.4:4]; % 机翼展宽8m，跳过 (-0.4, 0.4) 这个区间，让机翼点云从机身外壳两侧向外生长
    for y = y_wing_steps
        x_front = -abs(y)*0.5 + 0.5; % 机翼前边缘，向两侧时（|y|增大），x坐标向后减小（x轴正方向朝向机头）
        wing_chord = 1.5 * (1 - abs(y)/5); % 机翼宽度/弦长。机翼靠近机身的地方很宽，越往翼尖走越窄
        x_back = x_front - wing_chord; % 机翼的后边缘坐标用前缘坐标减去当前位置的机翼宽度
        X_local = [X_local, x_front, x_back];
        Y_local = [Y_local, y, y];
%         hold on;scatter(x_front, y, 30, 'b', 'filled');xlim([-6, 6]);ylim([-4.5, 4.5]);
%         hold on;scatter(x_back, y, 30, 'b', 'filled');xlim([-6, 6]);ylim([-4.5, 4.5]);
    end
    % --- 3. 箭形尾翼轮廓线 ---
    y_tail_steps = [-1.5:0.4:-0.4, 0.4:0.4:1.5]; % 尾翼总宽度为3m，同样跳过中间区间（尾翼点云接在机身尾部两侧）
    for y = y_tail_steps
        x_front_tail = -4.5 - abs(y)*0.3; % 后掠前缘线，尾翼是从距离机尾最后一米（-4.5 米处）开始向外生长的，向两侧时（|y|增大），x坐标向后减小
        tail_chord = 0.8 * (1 - abs(y)/2); % 尾翼宽度，根部最宽，尾翼尖端最小
        x_back_tail = x_front_tail - tail_chord; % 后边缘坐标用前缘坐标减去当前位置的机翼宽度
        X_local = [X_local, x_front_tail, x_back_tail];
        Y_local = [Y_local, y, y];
%         hold on;scatter(x_front_tail, y, 30, 'g', 'filled');xlim([-6, 6]);ylim([-4.5, 4.5]);
%         hold on;scatter(x_back_tail, y, 30, 'g', 'filled');xlim([-6, 6]);ylim([-4.5, 4.5]);
    end
    
    % --- 4. 封闭翼尖与尾部边界 ---
    X_local = [X_local, -1.5, -1.5, -5.5]; Y_local = [Y_local, 4, -4, 0];
%     hold on;scatter(-1.5, 4, 30, 'c', 'filled');xlim([-6, 6]);ylim([-4.5, 4.5]);
%     hold on;scatter(-1.5, -4, 30, 'c', 'filled');xlim([-6, 6]);ylim([-4.5, 4.5]);
%     hold on;scatter(-5.5, 0, 30, 'c', 'filled');xlim([-6, 6]);ylim([-4.5, 4.5]);
    
    %% --- 5. 姿态旋转矩阵 (引入了 heading_error) ---
    % 当 heading_error = 0 时，飞机就是绝对笔直地、自杀式地对着雷达天线冲过来。此时雷达只能看到一个正前方的“机头剪影”，机翼面积被严重压缩。
    % “机身点云就会显现”的原因： > 只要 heading_error 不为 0，飞机就是横着翅膀、侧着身子在做斜切进场。从雷达的视线（LOS）看过去，
    % 原本被隐藏的侧面机身、巨大的后掠翼表面全都在雷达视野里暴露无遗，散射点云瞬间变得极其丰富！

    % 以雷达为中心原点，飞机目标的y值始终是正值的，分左右一二象限分析，左二象限，
    % 飞机雷达径向速度为负时，机头朝着x轴正向，飞机径向速度为正时，机头朝着x轴负向，
    % 右一象限，飞机雷达径向速度为负时，机头朝着x轴负向,飞机径向速度为正时，时，机头朝着x轴正向
    % （实际中很复杂，因为飞行姿态导致的目标RCS是变化的。）
%     if center_angle <= 0
%         if v_total<=0
%             rot_angle = 90 - heading_error;
%         else
%             rot_angle = heading_error - 180;
%         end
%     else
%         if v_total<=0
%             rot_angle = heading_error - 360;
%         else
%             rot_angle = 90 - heading_error;
%         end
%     end
    % 根据旋转矩阵，顺时针旋转角度为负，逆时针为正的原理设置旋转角
    rot_angle = 90 - heading_error;
%     rot_angle = 180 + center_angle + heading_error; % 设置机头朝向雷达角度，上面设置飞机绝对位置时（x轴正方向朝向机头），这里旋转180度类似调头，雷达位于（0，0），机头就迎着雷达
    R_matrix = [cosd(rot_angle), -sind(rot_angle); % 二维逆时针旋转矩阵
                sind(rot_angle),  cosd(rot_angle)];
    
    % 把旋转后的飞机整体搬运到指定的空间位置
    % 角度定义：center_angle（方位角）是以 +Y 轴为 0° 基准。向右偏为正（x 变大），向左偏为负。
    rotated_points = R_matrix * [X_local; Y_local]; % 将飞机上的所有散射点绕着自身中心旋转rot_angle角度
%     figure(22);scatter(rotated_points(1,:), rotated_points(2, :), 30, 'k', 'filled');%xlim([-8, 8]);ylim([-8, 8]);
%     hold on;plot(0, 0, 'rp', 'MarkerSize', 12, 'MarkerFaceColor', 'r');set(gcf, 'Color', 'white');box on;axis equal;
    
    X_radar = rotated_points(1, :) + center_range * sind(center_angle);
    Y_radar = rotated_points(2, :) + center_range * cosd(center_angle);
%     figure(33);scatter(X_radar, Y_radar, 30, 'k', 'filled');%xlim([-8, 8]);ylim([-8, 8]);
%     hold on;plot(0, 0, 'rp', 'MarkerSize', 12, 'MarkerFaceColor', 'r');set(gcf, 'Color', 'white');box on;axis equal;
    
    

    %% --- 6. 解算回雷达参数 ---
    num_scatterers = length(X_radar);
    target_matrix = zeros(num_scatterers, 3);
    
%     flight_dir = center_angle - heading_error; 
    for i = 1:num_scatterers
        r_i = sqrt(X_radar(i)^2 + Y_radar(i)^2);
        % 每一个具体点（如左翼尖、右翼尖、机尾）【相对于雷达原点】的即时视线夹角，
        % 因为飞机有 10 米长、8 米宽，在近距离下，左翼尖和右翼尖的 angle_i 甚至能差出好几度！
        angle_i = atand(X_radar(i) / Y_radar(i)); 
        
        % 【多普勒径向速度投影】。核心物理修正：因为飞机斜着飞，各个散射点到雷达的投影速度不同
        v_radial_pure = abs(v_total) * cosd(angle_i - heading_error); %v_total * cosd(angle_i - flight_dir); 
        if v_total <= 0
            % 如果顶层设置 v_total 为负（靠近趋势），确保投影出的速度整体呈现负值
            v_radial = -abs(v_radial_pure); 
        else
            % 如果顶层设置 v_total 为正（远离趋势），确保呈现正值
            v_radial = abs(v_radial_pure); 
        end
        target_matrix(i, :) = [r_i, v_radial, angle_i];
    end
end