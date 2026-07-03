clear; clc; close all;
% 1. 雷达参数加载 
targetnum = 2; % 1-点目标  2-飞机目标
Parameter = generateParameter(targetnum);
numFrames = Parameter.numFrames;   % 帧数
T_frame = Parameter.T_frame;     % 帧周期
% 预分配一个元胞数组，用来记录每一帧 聚类 凝聚出来的 "干净质心点迹”
all_frames_centroids = cell(numFrames, 1);
% 备份初始目标设置，因为在帧循环中我们需要动态更新它们的距离
initial_targets = Parameter.target; % range speed horizontalangle TargetID
numTargets = size(initial_targets, 1); % 多目标组成的点目标总数，根据这些点目标的TargetID可得到实际目标数

%【解决海量回波爆内存，让外层循环每生成一帧 rawData，就立刻以追加（-append）的形式写入到硬盘的 .mat 文件里，然后迅速用 clear rawData 把内存腾出来】
saveDataFolder = 'data';
if ~exist(saveDataFolder, 'dir')
    mkdir(saveDataFolder); % 如果文件夹不存在则创建
end
saveFilePath = fullfile(saveDataFolder, 'radar_pipeline_inputs.mat');
save(saveFilePath, 'Parameter', '-v7.3'); % 把 配置参数 写入磁盘


fprintf('========= 开始进行多帧雷达回波信号仿真 =========\n');
all_frames_truth = cell(numFrames, 1);
% all_frames_rawData = cell(numFrames, 1); % 
all_frames_SNR = cell(1, numFrames); % 存储每一帧的SNR
%% 外层【帧】时间轴大循环
for frameId = 1:numFrames
    t_absolute = (frameId - 1) * T_frame; % 计算当前帧的宏观绝对时间墙 (秒)
    current_targets = initial_targets; % [Range, V_radial, Angle, TargetID]
    if targetnum == 1 % 点目标
        for t_idx = 1:numTargets
            R0 = initial_targets(t_idx, 1);
            v_radial = initial_targets(t_idx, 2); % 原始初始速度（带符号：靠近为负，远离为正）
            current_targets(t_idx, 1) = R0 + v_radial * t_absolute; 
            current_targets(t_idx, 3) = initial_targets(t_idx, 3); 
            current_targets(t_idx, 2) = v_radial; 
        end
    else % 飞机目标
        airplanes_config = Parameter.airplanes_config; % [距离, 目标中心与雷达角度, 航向恒定速度, 航向角]
        for t_idx = 1:numTargets
            R0   = initial_targets(t_idx, 1); 
            Th0  = initial_targets(t_idx, 3); % 初始方位角
            air_id = initial_targets(t_idx, 4);
            v_mag = abs(airplanes_config(air_id, 3));  % 该飞机的恒定总速率，直接取绝对值
            flight_dir  = airplanes_config(air_id, 4);  % 该飞机的恒定绝对航向角
            
            X0 = R0 * sind(Th0);
            Y0 = R0 * cosd(Th0);
            Xt = X0 + v_mag * sind(flight_dir) * t_absolute; % 航迹向量（起始原点）以Y轴正半轴为0°，速度为负表示朝着向量方向运动
            Yt = Y0 + v_mag * cosd(flight_dir) * t_absolute;
            current_targets(t_idx, 1) = sqrt(Xt^2 + Yt^2);     % 更新 实时变化的距离
            current_targets(t_idx, 3) = atan2d(Xt, Yt);        % 更新 实时变化的中心角度;
            
            unit_x = sind(Th0);
            unit_y = cosd(Th0);
            v_x = v_mag * sind(flight_dir);
            v_y = v_mag * cosd(flight_dir);
            v_radial_curr = v_x * unit_x + v_y * unit_y;
            current_targets(t_idx, 2) = v_radial_curr; % 更新 实时变化的径向速度
        end
    end
    Parameter.target = current_targets; % 将移动后的新位置灌回参数体
    all_frames_truth{frameId} = current_targets;
    
    % 2. 信号生成 将该帧下的所有点目标当成一个快照，生成该帧的仿真回波信号
    [rawData,SNR_scene_total_dB] = generateSignal(Parameter, targetnum);
%     all_frames_rawData{frameId} = rawData;
    all_frames_SNR{frameId} = SNR_scene_total_dB; % 存入数组
    varName = sprintf('frame_raw_%d', frameId); % 动态生成专属于当前帧的独立变量名 (如 frame_raw_1, frame_raw_2...)
    eval([varName ' = rawData;']);
    save(saveFilePath, varName, '-append'); % 即产即销，直接将这一帧的巨额回波“直传”追加到本地硬盘中
    eval(['clear ' varName ';']); 
    clear rawData;
    
    fprintf('进度: [Frame %d/%d] 已成功生成回波\n', ...
            frameId, numFrames);
end

fprintf('========= 仿真结束，正在将海量回波数据写入本地磁盘... =========\n');

% 保存后级算法需要的所有原材料：回波数据、SNR、真实轨迹、雷达硬件指标
% 把以下 4 份数据 装进 radar_pipeline_inputs.mat 里，-v7.3 会让 MATLAB 启动高级 HDF5 压缩算法
save(saveFilePath, 'all_frames_SNR', 'all_frames_truth', '-append');
fprintf('数据已成功保存至：%s。\n', saveFilePath);


