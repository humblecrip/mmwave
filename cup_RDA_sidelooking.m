%% 杯子数据四维处理

clear all; close all; clc;

%% 1. 参数设置

fprintf('=== 四维数据处理与重排 ===\n');

cfgFilePath = 'phone.mmwave.json';
setupFilePath = 'phone.config.setup.json';
dataFilePath = 'phone_30cm.bin';

hasJson = false;
fprintf('尝试读取推荐配置文件: %s\n', cfgFilePath);
if isfile(cfgFilePath)
    try
        mmwParams = jsondecode(fileread(cfgFilePath));
        hasJson = true;
        fprintf(' 成功读取推荐配置参数。\n');
    catch jsonErr
        fprintf(' 警告: 配置文件解析失败 (%s)，将尝试自动推断参数。\n', jsonErr.message);
    end
else
    fprintf(' 警告: 未找到配置文件，将尝试自动推断参数。\n');
end

fprintf('尝试读取采集配置文件: %s\n', setupFilePath);
if isfile(setupFilePath)
    try
        captureSetup = jsondecode(fileread(setupFilePath));
        if isfield(captureSetup, 'capturedFiles') && isfield(captureSetup.capturedFiles, 'files')
            filesInfo = captureSetup.capturedFiles.files;
            if ~isempty(filesInfo)
                basePath = '';
                if isfield(captureSetup.capturedFiles, 'fileBasePath')
                    basePath = captureSetup.capturedFiles.fileBasePath;
                end
                candidateNames = {filesInfo.processedFileName};
                if all(cellfun(@isempty, candidateNames)) && isfield(filesInfo, 'rawFileName')
                    candidateNames = {filesInfo.rawFileName};
                end
                for idx = 1:numel(candidateNames)
                    if isempty(candidateNames{idx})
                        continue;
                    end
                    candidate = candidateNames{idx};
                    if ~isempty(basePath)
                        candidateFull = fullfile(basePath, candidate);
                    else
                        candidateFull = candidate;
                    end
                    if isfile(candidateFull)
                        dataFilePath = candidateFull;
                        break;
                    end
                end
            end
        end
        fprintf(' 数据文件候选路径: %s\n', dataFilePath);
    catch setupErr
        fprintf(' 警告: 采集配置文件解析失败 (%s)。\n', setupErr.message);
    end
else
    fprintf(' 提示: 未找到采集配置文件，使用默认数据文件路径。\n');
end

% 如果推荐配置可用，提取关键参数
if hasJson
    rfConfig = mmwParams.mmWaveDevices(1).rfConfig;
    profileCfg = rfConfig.rlProfiles(1).rlProfileCfg_t;
    frameCfg = rfConfig.rlFrameCfg_t;
end

% 雷达系统参数初始化
c = 3e8; % 光速
if hasJson
    fc = mmwParams.regulatoryRestrictions.frequencyRangeBegin_GHz * 1e9; % 载波频率
    lambda = c / fc;
    chirpTime = (profileCfg.rampEndTime_usec + profileCfg.idleTimeConst_usec) * 1e-6;
    samplingRate = profileCfg.digOutSampleRate * 1e3; % JSON中单位为kHz
    slope_Hz_per_s = profileCfg.freqSlopeConst_MHz_usec * 1e12; % MHz/us -> Hz/s
    bandwidth = slope_Hz_per_s * profileCfg.rampEndTime_usec * 1e-6;
    frameTime = frameCfg.framePeriodicity_msec * 1e-3;
    numSamples = profileCfg.numAdcSamples;
    numChirps = frameCfg.numLoops;
    numFrames = frameCfg.numFrames;

    rxChannelEn = sscanf(rfConfig.rlChanCfg_t.rxChannelEn, '0x%x');
    txChannelEn = sscanf(rfConfig.rlChanCfg_t.txChannelEn, '0x%x');
    numRx = sum(bitget(uint32(rxChannelEn), 1:8));
    numTx = sum(bitget(uint32(txChannelEn), 1:3));

    range_resolution = c / (2 * bandwidth);
    PRF = 1 / chirpTime;
    velocity_resolution = (PRF / numChirps) * lambda / 2;
    % 天线间距设置（AWR2243 级联）
    % 方位角方向：相邻 TX 间距为 2λ
    d_tx_azimuth = 2 * lambda;
    % 俯仰方向：RX 阵列 A 与 C 间距 16λ，C 与 B 间距 4λ
    d_rx_elev_AC = 16 * lambda;
    d_rx_elev_CB = 4 * lambda;

    % 角度分辨率（用于后续角度轴近似换算，沿方位向）
    % 注意：这里沿用原脚本的简化公式，仅将原先 5 cm 间距替换为方位向的设计间距
    angle_resolution = lambda / (2 * d_tx_azimuth);

    fprintf('从配置文件提取的系统参数:\n');
    fprintf(' 载波频率: %.2f GHz\n', fc/1e9);
    fprintf(' 信号带宽: %.1f MHz\n', bandwidth/1e6);
    fprintf(' Chirp周期: %.2f us\n', chirpTime*1e6);
    fprintf(' 采样率: %.2f Msps\n', samplingRate/1e6);
    fprintf(' 距离分辨率: %.2f cm\n', range_resolution*100);
    fprintf(' 速度分辨率: %.2f m/s\n', velocity_resolution);
    fprintf(' 天线间距(方位TX): %.3f mm (2λ)\n', d_tx_azimuth*1e3);
    fprintf(' 天线间距(俯仰A-C): %.3f mm (16λ)\n', d_rx_elev_AC*1e3);
    fprintf(' 天线间距(俯仰C-B): %.3f mm (4λ)\n', d_rx_elev_CB*1e3);
    fprintf(' 采样点数: %d\n', numSamples);
    fprintf(' Chirp数量: %d\n', numChirps);
    fprintf(' 接收天线数: %d\n', numRx);
    fprintf(' 发送天线数: %d\n', numTx);
    fprintf(' 帧数: %d\n', numFrames);
else
    fprintf('无法使用配置文件，将根据数据自动推断参数。\n');
end

%% 2. 读取原始ADC数据

fprintf('\n尝试读取原始ADC数据...\n');

if ~isfile(dataFilePath)
    % 如果采集配置路径不可用，尝试当前目录下的默认文件
    alternatePath = fullfile(pwd, 'adc_data.bin');
    if isfile(alternatePath)
        fprintf(' 使用当前目录数据文件: %s\n', alternatePath);
        dataFilePath = alternatePath;
    else
        error('未找到数据文件。请确认 %s 是否存在或更新路径。', dataFilePath);
    end
end

fprintf(' 读取数据文件: %s\n', dataFilePath);
fileID = fopen(dataFilePath, 'r');
if fileID == -1
    error('无法打开数据文件: %s', dataFilePath);
end

fseek(fileID, 0, 'eof');
fileSize = ftell(fileID);
fseek(fileID, 0, 'bof');

rawData = fread(fileID, [fileSize/2, 1], 'int16');
fclose(fileID);

fprintf(' 文件大小: %.2f MB\n', fileSize/1024/1024);
fprintf(' 数据点数: %d\n', length(rawData));

% 如果没有JSON，尝试自动推断关键参数
if ~hasJson
fprintf('\n未提供JSON配置，开始自动推断参数...\n');
L = length(rawData);
candidateNumRx = [8 4 3 2 1];
candidateNumChirps = [256 192 128 96 64 48 32];
candidateNumSamples = [1024 768 512 384 256 192 128 96 64];
best = [];
for rx = candidateNumRx
for ch = candidateNumChirps
for ns = candidateNumSamples
den = ns * ch * rx * 2; % I/Q 交错
if mod(L, den) == 0
frames = L / den;
% 约束：帧数不为零且不过大
if frames >= 1 && frames <= 2000
score = [frames, ns, ch, rx];
best = [best; score]; %#ok<AGROW>
end
end
end
end
end
if isempty(best)
error('无法自动推断 (numSamples,numChirps,numRx)。请提供配置文件或手动参数。');
end
% 按帧数最大、采样点多排序，优先选择信息量大的组合
[~, idx] = sortrows(best, [-1 -2 -3 -4]);
sel = best(idx(1), :);
numFrames = sel(1);
numSamples = sel(2);
numChirps = sel(3);
numRx = sel(4);
% 设置典型物理参数以便估算分辨率（可按实际设备调整）
fc = 77e9; % 77 GHz 车载毫米波常用
lambda = c / fc;
chirpTime = 60e-6; % 60 us 典型
frameTime = 50e-3; % 50 ms 典型
samplingRate = 5e6; % 5 Msps 典型
bandwidth = 4e9; % 4 GHz 典型
range_resolution = c / (2 * bandwidth);
velocity_resolution = lambda / (2 * frameTime);
angle_resolution = lambda / (2 * 0.05);
PRF = 1 / frameTime;
fprintf(' 自动推断结果:\n');
fprintf(' 采样点数(numSamples): %d\n', numSamples);
fprintf(' Chirp数量(numChirps): %d\n', numChirps);
fprintf(' 接收天线数(numRx): %d\n', numRx);
fprintf(' 帧数(numFrames): %d\n', numFrames);
fprintf(' 假设: fc=%.1f GHz, 带宽=%.1f GHz, 帧周期=%.0f ms\n', fc/1e9, bandwidth/1e9, frameTime*1e3);
end

%% 3. 实现四维数据处理 - 将一维数据重组为[nSample, nLoops, nRx, nFrame]格式

fprintf('\n实现四维数据处理...\n');

% 检查数据量

expectedPointsPerFrame = numSamples * numChirps * numRx * 2; % 实部和虚部

totalExpectedPoints = expectedPointsPerFrame * numFrames;

if length(rawData) < totalExpectedPoints

fprintf(' 警告：数据点数不足，调整处理帧数\n');

numFrames = floor(length(rawData) / expectedPointsPerFrame);

if numFrames == 0

error(' 数据量不足以构成一帧');

end

fprintf(' 可处理帧数: %d\n', numFrames);

else

fprintf(' 处理帧数: %d\n', numFrames);

end

% --- 建议替换的正确代码 ---



% 数据指针

current_pos = 1;

% 创建四维矩阵

dataMat = zeros(numSamples, numChirps, numRx, numFrames);

% 循环填充四维矩阵

for i_frame = 1:numFrames

for i_chirp = 1:numChirps

for i_rx = 1:numRx

% 确定当前chirp所需的数据点数

points_to_read = numSamples * 2; % 实部+虚部


% 提取单个chirp的数据

chirp_data_1d = rawData(current_pos : current_pos + points_to_read - 1);


% 将I/Q交错数据转换为复数

% 奇数位是实部(I)，偶数位是虚部(Q)

complex_data = complex(chirp_data_1d(1:2:end), chirp_data_1d(2:2:end));


% 存入四维矩阵

dataMat(:, i_chirp, i_rx, i_frame) = complex_data;


% 更新数据指针

current_pos = current_pos + points_to_read;

end

end

end

fprintf(' 四维数据重构完成!\n');

%% 4. 三维FFT处理

fprintf('\n执行三维FFT处理...\n');

% 设置FFT点数，可以设为原始尺寸或者更大以提高分辨率

rangeFftSize = numSamples;

dopplerFftSize = numChirps;

angleFftSize = numRx * 2; % 角度FFT通常需要更多点以提高角度分辨率

% 选择处理第一帧数据

frameToProcess = 1;

frameData = dataMat(:,:,:,frameToProcess);

% 1. 一维FFT (距离维) + 2. 二维FFT (多普勒维)

fprintf(' 1&2. 执行距离/多普勒向FFT...\n');

% 应用 Blackman-Harris 窗（按照 GUI 配置）
if exist('blackmanharris','file')
    windowRange = blackmanharris(numSamples, 'periodic');
    windowDoppler = blackmanharris(numChirps, 'periodic');
else
    % 兼容性备用：若无 Blackman-Harris 函数则退化为 Blackman 窗
    windowRange = blackman(numSamples, 'periodic');
    windowDoppler = blackman(numChirps, 'periodic');
end

windowRange = reshape(windowRange, [], 1, 1);
windowDoppler = reshape(windowDoppler, 1, [], 1);

frameWindowed = frameData .* windowRange;
rangeData = fft(frameWindowed, rangeFftSize, 1);

rangeDataWindowed = rangeData .* windowDoppler;
rangeDopplerData = fftshift(fft(rangeDataWindowed, dopplerFftSize, 2), 2);

% 合并 RX 通道得到幅度谱（无非相干积累）
rdCube = sum(rangeDopplerData, 3);
rdMagnitude = abs(rdCube);
maxVal = max(rdMagnitude(:));
if maxVal == 0
    warning('Range-Doppler map contains only zeros.');
    maxVal = 1;
end
rdNorm = rdMagnitude ./ maxVal;
rdB = 20 * log10(rdNorm + eps);

rangeAxis = (0:rangeFftSize-1) * range_resolution;
dopplerAxis = ((-dopplerFftSize/2):(dopplerFftSize/2-1)) * velocity_resolution;

figure;
subplot(1,2,1);
meanRangeSpectrum = mean(mean(abs(rangeData), 3), 2);
plot(rangeAxis, 20*log10(meanRangeSpectrum + eps));
xlabel('距离 / m');
ylabel('幅度 / dB');
title('距离向平均谱');
grid on;

subplot(1,2,2);
midRangeIdx = round(linspace(1, rangeFftSize, 5));
midRangeIdx = midRangeIdx(2:4);
rdmSlice = squeeze(mean(abs(rdCube(midRangeIdx, :)), 1));
plot(dopplerAxis, 20*log10(rdmSlice + eps));
xlabel('速度 / (m/s)');
ylabel('幅度 / dB');
title('多普勒向平均谱');
grid on;

figure('Name','2D FFT amplitude profile');
imagesc(dopplerAxis, rangeAxis, rdB);
set(gca, 'YDir', 'normal');
xlabel('速度 / (m/s)');
ylabel('距离 / m');
title('距离-多普勒幅度谱');
colormap(jet);
cb = colorbar;
cb.Label.String = '幅度 (dB)';
cb.Ticks = [-120, -100, -80, -60, -40];
caxis([-120, -40]);
grid on;

% 3. 三维FFT (角度维)

fprintf(' 3. 执行角度向FFT...\n');

rangeDopplerAngleData = zeros(rangeFftSize, dopplerFftSize, angleFftSize);

for rangebin = 1:rangeFftSize

for dopplerbin = 1:dopplerFftSize

% 提取当前距离-多普勒单元的所有天线数据

antennaData = squeeze(rangeDopplerData(rangebin, dopplerbin, :));

% 对角度维执行FFT (扩展到angleFftSize点)

rangeDopplerAngleData(rangebin, dopplerbin, :) = fft(antennaData, angleFftSize);

end

end

% 生成三维热图，显示距离-多普勒-角度数据

fprintf(' 生成三维热图...\n');

% 计算三维数据的幅度

dataAmplitude = abs(rangeDopplerAngleData);

% 创建距离-角度投影图 (最大值投影)

rangeAngleProj = squeeze(max(dataAmplitude, [], 2));

figure;

subplot(2,2,1);

imagesc(rangeAngleProj);

title('距离-角度投影');

xlabel('角度Bin');

ylabel('距离Bin');

colorbar;

colormap(jet);

% 创建距离-多普勒投影图 (最大值投影)

rangeDopplerProj = squeeze(max(dataAmplitude, [], 3));

subplot(2,2,2);

imagesc(rangeDopplerProj);

title('距离-多普勒投影');

xlabel('多普勒Bin');

ylabel('距离Bin');

colorbar;

% 创建多普勒-角度投影图 (最大值投影)

dopplerAngleProj = squeeze(max(dataAmplitude, [], 1));

subplot(2,2,3);

imagesc(dopplerAngleProj);

title('多普勒-角度投影');

xlabel('角度Bin');

ylabel('多普勒Bin');

colorbar;

% 寻找三维数据中的峰值

[maxVal, linearIdx] = max(dataAmplitude(:));

[maxRangeBin, maxDopplerBin, maxAngleBin] = ind2sub(size(dataAmplitude), linearIdx);

% 将Bin索引转换为物理量

rangeMeters = (maxRangeBin-1) * range_resolution;

velocityMps = (maxDopplerBin - (dopplerFftSize/2 + 1)) * velocity_resolution;

angleDegrees = (maxAngleBin - angleFftSize/2) * (angle_resolution * 180/pi) / (angleFftSize/numRx);

% 显示目标信息

subplot(2,2,4);

text(0.1, 0.8, '目标检测结果:', 'FontSize', 14);

text(0.1, 0.6, sprintf('距离: %.2f 米', rangeMeters), 'FontSize', 12);

text(0.1, 0.5, sprintf('速度: %.2f 米/秒 (%.2f 公里/小时)', velocityMps, velocityMps*3.6), 'FontSize', 12);

text(0.1, 0.4, sprintf('角度: %.2f 度', angleDegrees), 'FontSize', 12);

text(0.1, 0.2, sprintf('峰值幅度: %.2f', maxVal), 'FontSize', 12);

axis off;

% 结果打印到命令窗口

fprintf('\n检测到目标:\n');

fprintf(' 距离: %.2f 米\n', rangeMeters);

fprintf(' 速度: %.2f 米/秒 (%.2f 公里/小时)\n', velocityMps, velocityMps*3.6);

fprintf(' 角度: %.2f 度\n', angleDegrees);

fprintf(' 峰值幅度: %.2f\n', maxVal);

fprintf(' 对应Bin索引: 距离=%d, 多普勒=%d, 角度=%d\n', maxRangeBin, maxDopplerBin, maxAngleBin);

%% 5. 三维可视化 (选择性展示)

fprintf('\n生成三维可视化...\n');

% 设置阈值，只显示高于阈值的数据点，以减少数据量

threshold = 0.1 * maxVal;

[rangeIdx, dopplerIdx, angleIdx] = ind2sub(size(dataAmplitude), find(dataAmplitude > threshold));

% 将Bin索引转换为物理量

ranges = (rangeIdx-1) * range_resolution;

velocities = (dopplerIdx - (dopplerFftSize/2 + 1)) * velocity_resolution;

angles = (angleIdx - angleFftSize/2) * (angle_resolution * 180/pi) / (angleFftSize/numRx);

% 强度值

intensities = zeros(length(rangeIdx), 1);

for i = 1:length(rangeIdx)

intensities(i) = dataAmplitude(rangeIdx(i), dopplerIdx(i), angleIdx(i));

end

% 创建3D散点图

figure;

scatter3(velocities, ranges, angles, 20, intensities, 'filled');

colormap(jet);

colorbar;

xlabel('速度 (m/s)');

ylabel('距离 (m)');

zlabel('角度 (度)');

title('目标三维散点图');

grid on;

% 标记最大峰值位置

hold on;

scatter3(velocityMps, rangeMeters, angleDegrees, 100, 'r', 'filled', 'MarkerEdgeColor', 'k');

text(velocityMps, rangeMeters, angleDegrees, ' 目标', 'Color', 'r');

hold off;

fprintf('\n处理完成!\n');
