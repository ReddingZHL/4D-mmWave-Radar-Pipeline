%% ==================== 4D毫米波雷达感知完整Pipeline ====================
clear; clc; close all;

addpath signal_processing
addpath post_processing
% 包含：all_frames_rawData, all_frames_SNR, all_frames_truth, Parameter
matFilePath = 'data/radar_pipeline_inputs_3Targets_50frames.mat';
load(matFilePath, 'Parameter', 'all_frames_truth', 'all_frames_SNR'); % 先加载轻量级雷达配置及目标先验信息

% 1. 雷达参数加载
targetnum = Parameter.targetnum; % 1-点目标  2-飞机目标
numFrames = Parameter.numFrames;   % 帧数
T_frame = Parameter.T_frame;     % 帧周期
% 预分配一个元胞数组，用来记录每一帧 聚类 凝聚出来的“干净质心点迹”
all_frames_centroids = cell(numFrames, 1);

% % 2. 生成/加载原始数据（可切换仿真或真实数据）
% dataFolderPath = "data";
% content = dir(dataFolderPath); % 获取文件夹内容列表
% % dir 函数会返回一个包含文件夹内容的结构体数组。需要注意的是，在 Windows 和 Linux 系统中，任何文件夹都默认包含两个隐藏的目录项：
% % (1) . （当前目录）
% % (2) .. （上级目录）
% % 因此，如果 dir 返回的结构体长度大于 2，则说明文件夹不为空。
% if length(content) > 2
% %     fname = "adc_data.bin"; 
% %     rawData = DCA1000_Read_Data(fname);
% else % 数据文件夹data内容为空，生成仿真数据并存储在data文件夹中
%     rawData = generateSignal(Parameter); % [虚拟通道数, 距离采样点, chirp数]
% end

% % 2.1 从rawData提取水平和俯仰阵列组
% % 在实际处理 IWR1843 的数据流时，从 rawData 中根据 txId 提取对应的接收数据。
% % 水平阵列组： 提取由 TX1 产生的 4 个接收通道和 TX3 产生的 4 个接收通道。
% tx1_rx = rawData(1:4, :, :);
% tx3_rx = rawData(9:12, :, :);
% horizontal_array_Rawdata = cat(1, tx1_rx, tx3_rx);
% % 俯仰阵列组： 提取由 TX2 产生的 4 个接收通道。
% tx2_rx = rawData(5:8, :, :);
% vertical_array_Rawdata = tx2_rx;

fprintf('========= 开始进行多帧雷达仿真信号点云处理 =========\n');

%% 外层【帧】时间轴大循环
for frameId = 1:numFrames
    current_targets = all_frames_truth{frameId}; % 提取当前帧真实位置
    Parameter.target = current_targets;
    % 2. 当前帧回波信号
    varName = sprintf('frame_raw_%d', frameId);
    fileData = load(matFilePath, varName); % 只从硬盘读取这 1 帧的数据 同时load所有matlab会爆内存
    rawData = fileData.(varName); % [虚拟通道数 (txNum*rxNum), rangeBin, chirps]
    clear fileData;
    SNR_scene_total_dB = all_frames_SNR{frameId};
    
    % 3. 信号处理链路
    % 3.1 range doppler 2dfft
    [range_fft, rd_map, RDMetrics] = range_doppler_fft(rawData, Parameter);          % Range + Doppler FFT
    rd_map_mean = squeeze(mean(abs(rd_map).^2, 1)); % 非相干积分二维能量图 (取模平方，通道平均)
    r_max = (Parameter.Fs * Parameter.c) / (2 * Parameter.Slope); % 最大理论距离
    range_axis = linspace(0, r_max, Parameter.rangeBin);
    % 计算速度轴时，必须使用整个 TDM 循环的总时间，TX1 的数据流在0, 3Tr, 6Tr, 9Tr...采样
    v_max = Parameter.lambda / (4 * Parameter.txNum * Parameter.Tc);
    v_res = Parameter.lambda / (2 * Parameter.Chirps * Parameter.txNum * Parameter.Tc);
    doppler_axis = linspace(-v_max, v_max - v_res, Parameter.dopplerBin); % 生成速度轴 (从 -v_max 到 v_max)
%     figure(1);imagesc(range(Parameter.dopplerBin), range_axis, squeeze(mean(abs(range_fft), 1)));
%     set(gca, 'YDir', 'normal');
%     xlabel('多普勒bin');ylabel('距离 (m)');title('Range fft');
%     figure(2);imagesc(doppler_axis, range_axis, rd_map_mean);
%     % 设置坐标轴方向（imagesc 默认 Y 轴是反的，距离 0 在上方）
%     set(gca, 'YDir', 'normal');
%     xlabel('速度 (m/s)');ylabel('距离 (m)');title('RD Map');
    
    
    % 3.2 cfar检测
%     [pointList, cfarRD] = cfar_detection_1D(db(rd_map_mean),target); % pointList记录了检测目标点的横纵维度值：纵轴为距离维度范围[1,256]，横轴为多普勒维度范围[1,512]
%     [pointList, cfarRD] = cfar_detection_1D_pro(rd_map_mean, targetnum); % 送入CFAR前不必加db，直接就是功率谱
    [pointList, cfarRD] = cfar_detection_2D(rd_map_mean, targetnum);
    if isempty(pointList)
        fprintf('第 %d/%d 帧：未检测到任何 CFAR 点迹。\n', frameId, numFrames);
        all_frames_centroids{frameId} = [];
        continue;
    end
%     figure(3);
%     imagesc(doppler_axis',range_axis,cfarRD);
%     set(gca, 'YDir', 'normal');
%     xlabel('速度(m/s)'); ylabel('距离(m)'); title("CFAR RD");set(gcf, 'Color', 'white');
    
    
    % 3.3 方位角度估计1D-MUSIC
%     % ① 正前方的 3dB 物理波束宽度 (瑞利限分辨率)
%     ang_res_boresight_rad = Parameter.lambda / (Parameter.virtualAntenna * Parameter.dx);
%     ang_res_boresight_deg = rad2deg(ang_res_boresight_rad) * 0.886; % 引入主瓣因子近似
%     ang_res = ang_res_boresight_deg; % [正前方] 测角分辨率
%     theta_edge = 45; 
%     ang_res_edge_deg = ang_res_boresight_deg / cosd(theta_edge); % 算大角度边缘（比如偏向 45 度时）的退化分辨率
%     max_fov_deg = 2 * asind(Parameter.lambda / (2 * Parameter.dx)); % 最大不模糊测角范围 (FoV)

    
    % 【记得，传入的是RD图（经过加窗FFT，旁瓣会得到抑制），而不是原始数据！！！】
    % 【仿真的飞机形状有将近一百个点，降低d_SNR和r_SNR，cfar输出的点从5个锐增到50个】
    % 【增加Tr，相当于增加带宽，提高分辨率，cfar输出点增加到80个】
    % 【给飞机加个倾斜角度，而不是直直朝着雷达飞，这样机身能出来，因为速度维的很多点分离了】
    doaMethod = 3; % parameter.doaMethod 选择测角算法 【1:fft  2:dbf  3:1Dmusic  4:RDMusic  5:Capon】
    detections_4d = doa(doaMethod, rd_map, pointList, cfarRD, Parameter); % 输出 [range, vel, angle, power]
    
    % 4. 点云后处理
    % 4.1 点云生成 极坐标转笛卡尔空间坐标
    pointCloud = pointcloud_generation(detections_4d); % 输出 [x y z 径向v power]
    
    % 4.2 点云聚类
    % 【把聚类方法从X-Means改成DBSCAN后，就不存在质心偏离飞机目标的问题了，现在所有质心都在飞机目标上】
    clusterMethod = 1; % parameter.clusterMethod 选择点云聚类方法 1-DBSCAN 2-Kmeans（计划的是3目标） 3-X-Means
    clusters = points_clustering(clusterMethod, pointCloud, targetnum); % [x, y, z, 径向v, power, Cluster_ID]
    
    % 4.3 点迹滤波，提炼点云聚类质心
    [centroids,ClusterMetrics] = extract_centroids(clusters, Parameter.target); % [x, y, z, 径向v, power]
    all_frames_centroids{frameId} = centroids;% 存储进多帧大仓库中
    % 【绘制点云及起聚类质心，有的质心不在飞机上，可设置卡尔曼滤波 "波门"，
    % 滤波器在预测下一帧飞机位置时，会以预测点为中心，画一个半径为 R 的圆（或椭圆）叫做跟踪波门。
    % 只有落在波门内部的质心，才允许用来更新卡尔曼状态】
    hold on;
    plot(centroids(:,1), centroids(:,2), 'rp', ...
     'MarkerSize', 14, ...
     'MarkerFaceColor', 'r', ...
     'MarkerEdgeColor', 'k', 'LineWidth', 1);
    
    % 【聚类后，速度解模糊，解模糊方法有多种，一种是假设检验，假设目标折叠了 -1, 0, 1 次，分别算出三个对应的真实速度相位，
    % 用这三个相位分别对天线向量 antVec 进行 TDM 补偿，把补偿后的三个不同向量，分别送进测角算法
    % （如 DBF 或 MUSIC）中计算能量谱，只有当假设的折叠次数与目标真实的物理运动完全一致时，
    % 天线阵列的孔径相位才会被完美修复，此时测角谱的能量峰值最高】
    % 【发射端时序解模糊法（工业前向雷达最常用的硬件级解法），奇/偶数 Chirp 序列采用 Tc_1/Tc_2，
    % 速度将折叠到不同的频点，只需要将这两点代入中国剩余定理，就可推出真实物理速度，但是这样就不能直接做2DFFT了，得拆分成奇/偶数矩阵，
    % 还有TDM-MIMO通道的多普勒引起的相位偏移也变了，在做doa前的相位补偿时，需要判断当前点迹来自奇数还是偶数时序动态切换补偿因子】
    
    eval_gates.r   = 5;   % 距离门限 米
    eval_gates.v   = 1;  % 速度门限 m/s
    eval_gates.ang = 1.5;3;  % 角度门限 度
    total_radar_cells = Parameter.rangeBin * Parameter.dopplerBin;
    CfarDoaMetrics = calculate_radar_metrics(Parameter.target, centroids, eval_gates, total_radar_cells);
    

    current_frame_snr_struct = all_frames_SNR{frameId}; 
    num_objs_in_frame = length(current_frame_snr_struct);
    % 拼接每个完整目标的 SNR 字符串 (例如: "ID1:+25.3dB | ID2:+18.1dB")
    snr_str_cells = cell(1, num_objs_in_frame);
    for obj_i = 1:num_objs_in_frame
        snr_str_cells{obj_i} = sprintf('Obj%d : %+5.1f dB', ...
            current_frame_snr_struct(obj_i).TargetID, ...
            current_frame_snr_struct(obj_i).Total_SNR_dB);
    end
    all_objs_snr_str = strjoin(snr_str_cells, ' | '); % 用竖线隔开
    
    fprintf('\n================================== [ 第 %2d/%2d 帧 ] ==================================\n', frameId, numFrames);
    fprintf('   ├── [时域特征]  目标群SNR: [ %s ] \n', all_objs_snr_str);
    fprintf('   ├── [频域特性]  理论距离分辨率: %5.2f 米 | 加窗距离分辨率: %5.2f 米 | 理论速度分辨率: %5.4f 米/秒 | 加窗速度分辨率: %5.4f 米 \n', ...
            RDMetrics.res_r_theory, RDMetrics.res_r_actual, RDMetrics.res_v_theory, RDMetrics.res_v_actual);
    fprintf('   ├── [聚类全场]  簇数: %d | 均密: %5.1f 点/m³ | 均紧: %4.2f 米 | 均偏: %5.4f 米 \n', ...
            size(centroids, 1), ClusterMetrics.avg_density, ClusterMetrics.avg_compact, ClusterMetrics.avg_center_err);
    fprintf('   └── [信号检测]  目标检出率: %6.2f%% | 虚警率: %10.4e | 虚警目标数: %d | 角度均方误差: %5.3f° \n', ...
            CfarDoaMetrics.Cluster_Pd * 100, CfarDoaMetrics.Cluster_Pfa, CfarDoaMetrics.Cluster_False_Alarm_Count, CfarDoaMetrics.angle_rmse_pure);
    drawnow;
end

fprintf('========= 所有帧信号模拟完成，进入卡尔曼多目标追踪器 =========\n');

%% 5. 全局卡尔曼数据关联与航迹生命周期管理
% 5.1 传入整段历史点迹，自动孵化新目标、平滑老目标、销毁死目标
[tracked_objects, SummaryMetrics] = kalman_tracking(all_frames_centroids, Parameter, all_frames_truth);% 多目标跟踪
% 卡尔曼滤波加入检测到的速度信息，跟踪轨迹更加平滑，速度估计误差大幅度下降
% [tracked_objects, SummaryMetrics] = kalman_tracking_xyv(all_frames_centroids, Parameter, all_frames_truth);% 多目标跟踪

fprintf('===========================================================\n');
fprintf('  [状态平滑度] 距离维 RMSE : %5.4f 米\n', SummaryMetrics.Range_RMSE);
fprintf('              径向速度维 RMSE : %5.4f m/s\n', SummaryMetrics.Vel_RMSE);
fprintf('              角度维 RMSE : %5.4f °\n', SummaryMetrics.Angle_RMSE);
fprintf('  [连续性评估] 航迹整体存活率: %6.2f%%\n', SummaryMetrics.Mean_Track_Survival_Rate * 100);
fprintf('  [大一统评估] 全场平均 OSPA : %5.4f 米 (阶数p=2,惩罚c=10)\n', SummaryMetrics.Mean_OSPA);
fprintf('===========================================================\n');

% 5.2 绘制连续运动轨迹
plot_tracking_results(tracked_objects, Parameter);
% disp('4D雷达感知Pipeline运行完成！');


