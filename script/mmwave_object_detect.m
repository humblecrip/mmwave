%% main_tdm2tx_demo.m  —— add CA-CFAR detection (Guard=8, Win=32, Thr=18 dB)
clear; close all; clc;

%% ====== 基本参数（按你的实际配置填写） ======
SampleRate     = 5e6;           % Fs, ADC采样率 [Hz]
FrequencySlope = 30e12;         % 调频斜率 [Hz/s]
Nr  = 256;                      % 每个 chirp 的采样点数
Nd  = 128;                      % 每帧 chirp 总数(包含所有Tx)
Nrx = 4;                        % 接收天线数
Ntx = 2;                        % 发射天线数 (TDM-MIMO)
fc = 77.94745e9;                % 载频 [Hz]
c  = 3e8;                       % 光速
lambda = c/fc;
Tc = 70e-6;                     % 单个chirp周期(含idle) [s]
dataFile = './reflect_right.bin';
DO_OFFSET_SEARCH = true;

%% ====== 从 JSON 读取参数（覆盖上面的默认值） ======
jsonFile = './0922v2.mmwave.json';
try
    cfg = jsondecode(fileread(jsonFile));
    dev  = cfg.mmWaveDevices(1).rfConfig;
    prof = dev.rlProfiles(1).rlProfileCfg_t;
    frm  = dev.rlFrameCfg_t;

    % 单位换算
    SampleRate     = prof.digOutSampleRate * 1e3;         % ksps → Hz
    FrequencySlope = prof.freqSlopeConst_MHz_usec * 1e12; % MHz/us → Hz/s
    Nr             = prof.numAdcSamples;                  % 采样点/Chirp
    fc             = prof.startFreqConst_GHz * 1e9;       % GHz → Hz
    lambda         = c / fc;
    Tc             = (prof.idleTimeConst_usec + prof.rampEndTime_usec) * 1e-6; % us → s

    % 通道掩码 → 天线数
    rxMask = lower(strrep(dev.rlChanCfg_t.rxChannelEn,'0x',''));
    txMask = lower(strrep(dev.rlChanCfg_t.txChannelEn,'0x',''));
    Nrx = sum(dec2bin(hex2dec(rxMask))=='1');
    Ntx = sum(dec2bin(hex2dec(txMask))=='1');

    % 配置侧：每帧理论 chirp 数（用于对比，不强制）
    chirpsPerLoop_cfg  = (frm.chirpEndIdx - frm.chirpStartIdx + 1);
    chirpsPerFrame_cfg = chirpsPerLoop_cfg * frm.numLoops;

    % 打印参数
    fprintf('==== 从 JSON 读取参数 ====\n');
    fprintf('  载频 fc: %.3f GHz\n', fc/1e9);
    fprintf('  采样率 Fs: %.3f Msps\n', SampleRate/1e6);
    fprintf('  调频斜率: %.6f MHz/us\n', prof.freqSlopeConst_MHz_usec);
    fprintf('  Nr=%d, Nrx=%d, Ntx=%d, Tc=%.2f us\n', Nr, Nrx, Ntx, Tc*1e6);
    fprintf('  每帧 chirp(配置) = %d\n', chirpsPerFrame_cfg);
catch ME
    warning('读取 JSON 失败，沿用默认参数。原因: %s', ME.message);
end

%% ====== CFAR 参数 ======
CFAR_guard  = 8;                % Guard window size
CFAR_window = 32;               % Training window (每侧)
CFAR_thr_dB = 18;               % 门限(dB) —— x > alpha * 均值
alpha = 10^(CFAR_thr_dB/10);

%% ====== 读取原始数据：返回 [Nrx x Nsamp] ======
adc = readDCA1000(dataFile);    % 默认 16bit, 4 lanes, 复数IQ

%% ====== 依据“总chirp/帧数”计算单帧chirp数，并覆盖 Nd ======
try
    if exist('frm','var') && isfield(frm,'numFrames') && frm.numFrames > 0
        totalChirps_inFile = floor(size(adc,2) / Nr);
        Nd_data = floor(totalChirps_inFile / frm.numFrames);
        if exist('chirpsPerFrame_cfg','var') && Nd_data ~= chirpsPerFrame_cfg
            fprintf('[WARN] 单帧 chirp(数据)=%d ≠ 配置=%d，按数据取值。\n', Nd_data, chirpsPerFrame_cfg);
        else
            fprintf('单帧 chirp(数据)= %d\n', Nd_data);
        end
        Nd = Nd_data;
    else
        fprintf('[INFO] JSON 未提供 numFrames，保持 Nd=%d\n', Nd);
    end
catch ME
    warning('按“总chirp/帧数”计算 Nd 失败，保持 Nd=%d。原因: %s', Nd, ME.message);
end

%% ====== 可选：粗搜最佳起点偏移（0..Nr-1 样点） ======
bestOff = 0;
if DO_OFFSET_SEARCH
    numChirps = Nd;
    Ns_need = Nr*numChirps;
    assert(size(adc,2) >= Ns_need+Nr, '数据长度不足以做offset搜索');
    bestScore = -inf;
    for off = 0:Nr-1
        adcFrame = adc(:, 1+off : off+Ns_need);
        tmp = permute( reshape(adcFrame, [Nrx, Nr, numChirps]), [2 3 1] ); % [Nr,Nd,Nrx]
        RFFT_test = fft(tmp .* reshape(hann(Nr),[Nr 1 1]), Nr, 1);
        prof = squeeze(mean(mean(abs(RFFT_test(1:floor(Nr/2),:,:)),3),2));  % 距离剖面
        score = max(prof);
        if score > bestScore
            bestScore = score; bestOff = off;
        end
    end
    fprintf('自动选择的起点偏移 = %d 样点\n', bestOff);
end

%% ====== 截取一帧并 reshape 到 [Nr, Nd, Nrx] ======
numChirps = Nd;                            % 总 chirp 数（包含所有Tx）
adcFrame  = adc(:, 1+bestOff : bestOff+Nr*numChirps);
tmp       = permute( reshape(adcFrame, [Nrx, Nr, numChirps]), [2 3 1] ); % [Nr,Nd,Nrx]

%% ====== 按 2Tx 的 TDM 顺序拆分慢时间到每个 Tx ======
Nd_perTx = numChirps / Ntx;
assert(mod(numChirps,Ntx)==0, 'Nd(%d) 不是 Ntx(%d) 的整数倍', numChirps, Ntx);
cube = zeros(Nr, Nd_perTx, Nrx, Ntx);
for tx = 1:Ntx
    chirpIdx = tx:Ntx:numChirps;     % Tx序列: Tx1,Tx2,Tx1,Tx2,...
    cube(:,:,:,tx) = tmp(:, chirpIdx, :);   % [Nr,Nd_perTx,Nrx]
end

%% ====== 1D 距离向 FFT ======
range_win = hann(Nr);
cube_win  = cube .* reshape(range_win, [Nr 1 1 1]);

NfftR = Nr;                                    % 可按需零填充
RFFT  = fft(cube_win, NfftR, 1);               % [NfftR, Nd_perTx, Nrx, Ntx]

% 只取正频半谱
Nr_pos = floor(NfftR/2);
RFFT   = RFFT(1:Nr_pos, :, :, :);              % [Nr_pos, Nd_perTx, Nrx, Ntx]

% 物理距离轴
range_axis = (0:Nr_pos-1) * (SampleRate/NfftR) * c / (2*FrequencySlope);

%% ====== 2D 多普勒/速度向 FFT（沿慢时间） ======
%RFFT_dc = RFFT - mean(RFFT, 2);                % 去静杂波
RFFT_dc = RFFT;
dop_win = hann(Nd_perTx).';
RFFT_w  = RFFT_dc .* reshape(dop_win, [1 Nd_perTx 1 1]);

NfftD = Nd_perTx;
RD = fftshift(fft(RFFT_w, NfftD, 2), 2);       % [Nr_pos, NfftD, Nrx, Ntx]

% 速度轴（注意 TDM 有效PRF = 1/(Ntx*Tc)）
PRF_perTx = 1/(Ntx*Tc);
fD_axis   = (-NfftD/2:NfftD/2-1)/NfftD * PRF_perTx;
v_axis    = (lambda/2) * fD_axis;

% 速度分辨率与不模糊速度
v_res = (lambda/2) * (PRF_perTx/NfftD);
v_max = (lambda/2) * (PRF_perTx/2);
fprintf('速度分辨率: %.4f m/s, 最大不模糊速度: ±%.4f m/s\n', v_res, v_max);

%% ====== Range–Doppler 图（非相干合成 Rx/Tx） ======
% 画幅度图（保持你原有显示）
RD_mag = squeeze(sum(sum(abs(RD), 3), 4));     % [Nr_pos, NfftD]
RD_dB  = 20*log10(RD_mag + eps);
figure(1); imagesc(v_axis, range_axis, RD_dB); axis xy;
xlabel('速度 (m/s)'); ylabel('距离 (m)'); title('Range–Doppler（非相干合成Rx/Tx）');
cb = colorbar; cb.Label.String = '幅度 (dB)'; hold on;

%% ====== —— 新增：CA-CFAR 目标检测（在功率图上做） ======
RD_pow = squeeze(sum(sum(abs(RD).^2, 3), 4));  % 功率图 [Nr_pos, NfftD]

% 每个多普勒bin做 1D CA-CFAR（沿距离维）
detMask = false(size(RD_pow));
for d = 1:NfftD
    x = RD_pow(:, d);                          % 该多普勒列的功率
    detMask(:, d) = cfar_ca_1d(x, CFAR_window, CFAR_guard, alpha);
end

% 抑制邻近多次命中，只保留局部最大
detMask = detMask & (imregionalmax(RD_pow));

% 取出检测点并叠加到 RD 图
[idx_r, idx_d] = find(detMask);
plot(v_axis(idx_d), range_axis(idx_r), 'ro', 'MarkerSize', 6, 'LineWidth', 1.2);
legend('CFAR detections');

% 打印检测到的目标 (距离, 速度)
for k = 1:numel(idx_r)
    fprintf('Detection %2d: Range = %.2f m, Velocity = %.2f m/s\n', ...
        k, range_axis(idx_r(k)), v_axis(idx_d(k)));
end




%% ====== 角度维 FFT（虚拟阵列） ======
virt_spacing = lambda/2;                       % 默认虚拟阵元间距取半波长
angleBins = Nrx * Ntx;                         % 角度FFT点数，等于虚拟阵元总数
angleWin = hann(angleBins);                    % 角度向窗函数
[angleCube, angleAxis] = compute_angle_fft_cube(RD, Nrx, Ntx, lambda, virt_spacing, angleBins, angleWin);

%% ====== 目标角度估计与3D可视化 ======
detAngles = [];
detRanges = [];
if exist('angleCube','var') && ~isempty(idx_r)
    numDet = numel(idx_r);
    detAngles = zeros(numDet, 1);
    for k = 1:numDet
        angSpec = squeeze(angleCube(idx_r(k), idx_d(k), :));
        [~, angIdx] = max(abs(angSpec));
        detAngles(k) = angleAxis(angIdx);
        fprintf('检测目标 %2d: 角度 = %.2f 度\n', k, detAngles(k));
    end
    detRanges = range_axis(idx_r);

    thetaRad = detAngles * pi/180;
    x_det = detRanges .* sin(thetaRad);         % 横向坐标
    y_det = detRanges .* cos(thetaRad);         % 前向坐标
    z_det = zeros(numDet, 1);                   % 目前假设目标位于地面，高度设为0

    fig3d = figure;
    hDet = plot3(x_det, y_det, z_det, 'o', 'LineStyle', 'none', 'MarkerSize', 8);
    set(hDet, 'MarkerFaceColor', [0 1 0], 'MarkerEdgeColor', [0 0.5 0]);
    hold on;
    hRadar = plot3(0, 0, 0, 'v', 'LineStyle', 'none', 'MarkerSize', 9);
    set(hRadar, 'MarkerFaceColor', [0 0 1], 'MarkerEdgeColor', [0 0 0.5]);
    grid on;
    xlabel('横向 x (m)');
    ylabel('前向 y (m)');
    zlabel('高度 z (m)');
    title('3D 目标检测结果');
    xRange = max(20, ceil(max(abs(x_det)) + 2));
    yRange = max(range_axis);
    xlim([-xRange, xRange]);
    ylim([0, max(25, ceil(yRange))]);
    zlim([-2, 5]);
    view([-35 35]);
    legend([hDet, hRadar], {'检测目标', '雷达位置'}, 'Location', 'best');
end



%% （保留）每个 Tx 的距离剖面 / Range–Doppler
for tx = 1:Ntx
    RP_tx = squeeze(mean(sum(abs(RFFT(:,:,:,tx)),3),2));  % [Nr_pos x 1]
    figure; plot(range_axis, 20*log10(RP_tx + eps)); grid on;
    xlabel('距离 (m)'); ylabel('幅度 (dB)');
    title(sprintf('距离剖面（Tx%d，非相干合成Rx）', tx-1));
end

for tx = 1:Ntx
    RD_tx = squeeze(sum(abs(RD(:,:,:,tx)),3));  % [Nr_pos, NfftD]
    figure; imagesc(v_axis, range_axis, 20*log10(RD_tx + eps)); axis xy;
    xlabel('速度 (m/s)'); ylabel('距离 (m)');
    title(sprintf('Range–Doppler（Tx%d，非相干合成Rx）', tx-1));
    cb = colorbar; cb.Label.String = '幅度 (dB)';
end

%% ==================== 子函数：1D CA-CFAR（均值） ====================
function detMask = cfar_ca_1d(x, win, guard, alpha)
% x      : 列向量(功率)
% win    : 每侧训练单元数
% guard  : 每侧保护单元数
% alpha  : 阈值系数（>1），等效 10^(门限_dB/10)
N = numel(x);
detMask = false(N,1);
for i = 1:N
    l1 = max(1, i - guard - win);
    l2 = max(1, i - guard - 1);
    r1 = min(N, i + guard + 1);
    r2 = min(N, i + guard + win);
    train = [x(l1:l2); x(r1:r2)];
    if numel(train) < 4, continue; end
    mu = mean(train);
    if x(i) > alpha * mu
        detMask(i) = true;
    end
end
end



function [angleCube, angleAxis] = compute_angle_fft_cube(RD, Nrx, Ntx, lambda, d_elem, NfftAngle, angleWin)
% RD        : 尺寸为 [rangeBin, dopplerBin, Nrx, Ntx] 的复数数据立方体
% Nrx/Ntx   : 真实接收与发射天线数
% lambda    : 工作波长，用于角度轴计算
% d_elem    : 相邻虚拟阵元间距，默认按半波长处理
% NfftAngle : 角度FFT点数，默认等于虚拟阵元总数
% angleWin  : 角度向窗函数（列向量），默认使用全 1
if nargin < 5 || isempty(d_elem)
    d_elem = lambda/2;
end
if nargin < 6 || isempty(NfftAngle)
    NfftAngle = Nrx * Ntx;
end
if nargin < 7 || isempty(angleWin)
    angleWin = ones(Nrx * Ntx, 1);
end

virtualElem = Nrx * Ntx;
assert(size(RD,3) == Nrx && size(RD,4) == Ntx, 'RD 尺寸与天线数量不匹配');
assert(numel(angleWin) == virtualElem, '角度窗长度需等于虚拟阵元数量');

virtCube = reshape(permute(RD, [1 2 4 3]), size(RD,1), size(RD,2), virtualElem);  % [range, doppler, virtual]
angleWin = reshape(angleWin(:).', [1 1 virtualElem]);
virtCube = virtCube .* angleWin;

angleCube = fftshift(fft(virtCube, NfftAngle, 3), 3);

spatialFreq = (-NfftAngle/2 : NfftAngle/2 - 1) / NfftAngle;
sinTheta = spatialFreq * lambda / d_elem;
sinTheta = max(min(sinTheta, 1), -1);
angleAxis = asind(sinTheta);
end
