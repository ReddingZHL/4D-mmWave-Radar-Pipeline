function [pointList,cfarRD] = cfar_detection_1D(accumulateRD, target)

    [rangeLen,dopplerLen] = size(accumulateRD);
    %% doppler维度搜索
    dopplerCfarMethod = 1; % 1:CA-CFAR  2:GO-CFAR    3:SO-CFAR
    if target == 1 % 1-点目标
        dopplerSNR = 10; % 指目标必须比周围噪声高出多少分贝（dB）才会被判定为5发现目标1
        dopplerWinGuardLen = 2; % 紧邻 CUT 的区域（WinGuardLen）。目的是防止目标的能量泄露到参考单元中，导致门限被错误抬高。
        dopplerWinTrainLen = 8; % 最外圈的区域（WinTrainLen）。雷达通过这些单元的数据来估算当前的背景噪声水平。
    
    else % 飞机参数
        dopplerSNR  = 10;5; % 多普勒维参数
        dopplerWinGuardLen = 2; % 紧邻 CUT 的区域（WinGuardLen）。目的是防止目标的能量泄露到参考单元中，导致门限被错误抬高。
        dopplerWinTrainLen = 8; % 最外圈的区域（WinTrainLen）。雷达通过这些单元的数据来估算当前的背景噪声水平。
    
    end
    
    
    % 环形填充
    % 雷达通过连续的脉冲串之间的相位变化测速。相位移动量超过2pi，就会回绕到0
    dopplerLeft = accumulateRD(:,dopplerLen - dopplerWinGuardLen - dopplerWinTrainLen + 1:dopplerLen);
    dopplerRight = accumulateRD(:,1:dopplerWinGuardLen+dopplerWinTrainLen);
    dopplercfar = [dopplerLeft accumulateRD dopplerRight];

    dopplerCfarList = [];
    for rangeIdx = 1:rangeLen
        for dopplerIdx = 1:dopplerLen
            % 待检测单元 (CUT, Cell Under Test)：中心点，判断它是目标还是噪声。
            dopplerCfarIdx = dopplerIdx + dopplerWinGuardLen + dopplerWinTrainLen;
            leftCell = dopplercfar(rangeIdx,dopplerIdx:dopplerIdx+dopplerWinTrainLen-1);
            rightCell = dopplercfar(rangeIdx,dopplerCfarIdx+dopplerWinGuardLen:dopplerCfarIdx+dopplerWinGuardLen+dopplerWinTrainLen-1);
            leftNoise = mean(leftCell); % db域的平均数，算出来的背噪偏低
            rightNoise = mean(rightCell);
            if dopplerCfarMethod == 1
                noise = (leftNoise + rightNoise) / 2;
            elseif dopplerCfarMethod == 2
                % 目的：抑制边缘虚警，当一侧有强干扰时，选较大的一侧作为门限，门限变高，更不容易误报。
                noise = max(leftNoise, rightNoise);
            elseif dopplerCfarMethod == 3
                % 目的：减少漏检，当一侧有强干扰时，选较小的一侧作为门限，门限变低，更容易检测到微弱目标。
                noise = min(leftNoise, rightNoise);
            end
            indexdb = dopplercfar(rangeIdx,dopplerCfarIdx);
            targetSnr = indexdb - noise;
            if  targetSnr > dopplerSNR
                dopplerCfarList = [dopplerCfarList dopplerIdx];
%                 cfarRDdoppler(rangeIdx,dopplerIdx) = indexdb;
            end
        end
    end
    dopplerCfarList = unique(dopplerCfarList);

    %% range维度搜索
    rangeCfarMethod = 1; % 1:CA-CFAR  2:GO-CFAR    3:SO-CFAR 
    if target == 1 % 1-点目标
        rangeSNR = 10; 
        rangeWinGuardLen =  4; % 根据参数确定的距离分辨率，计算一个普通人体或车辆目标在距离轴上横跨的range bin，设定时需要考虑目标能量泄露扩展，将距离维保护单元扩大或缩小
        rangeWinTrainLen =  8;
    else % 飞机参数
        rangeSNR  = 10;5; 
        rangeWinGuardLen = 4; % 2; % 根据参数确定的距离分辨率，计算一个普通人体或车辆目标在距离轴上横跨的range bin，设定时需要考虑目标能量泄露扩展，将距离维保护单元扩大或缩小
        rangeWinTrainLen = 8; % 4;
    end
    
    
    % 对称/镜像填充，距离维度不具备多普勒维度"首尾相接"的周期性
    pad_size = [rangeWinGuardLen + rangeWinTrainLen, 0]; % 参数1：距离维度方向填充量；参数2：多普勒方向填充量
    rangecfar = padarray(accumulateRD, pad_size, 'symmetric', 'both'); % 做法是：以原数组的边界线为对称轴，将边缘数据"翻折"过来。both参数表示在数据上下两端均应用该填充
    
    rangeCfarList = [];
    cfarRD = zeros(rangeLen,dopplerLen);
    for dopplerIdx = dopplerCfarList % 减少搜索次数
        for rangeIdx = 1:rangeLen
            rangeCfarIdx = rangeIdx + rangeWinGuardLen + rangeWinTrainLen;
            upCell = rangecfar(rangeIdx:rangeIdx+rangeWinTrainLen-1,dopplerIdx);
            downCell = rangecfar(rangeCfarIdx+rangeWinGuardLen:rangeCfarIdx+rangeWinGuardLen+rangeWinTrainLen-1,dopplerIdx);
            upNoise = mean(upCell);
            downNoise = mean(downCell);
            if rangeCfarMethod == 1
                noise = (upNoise + downNoise) / 2;
            elseif rangeCfarMethod == 2
                noise = max(upNoise, downNoise);
            elseif rangeCfarMethod == 3
                noise = min(upNoise, downNoise);
            end
            indexdb = rangecfar(rangeCfarIdx,dopplerIdx);
            targetSnr = indexdb - noise;
            if  targetSnr > rangeSNR
                rangeCfarList = [rangeCfarList [rangeIdx;dopplerIdx]];
                cfarRD(rangeIdx,dopplerIdx) = indexdb;
            end
        end
    end
    pointList = rangeCfarList;
end
