clear; close all; clc;

%% 基本参数（按你的配置）
FrequencySlope = 30e12;      % 调频斜率 Hz/s
SampleRate     = 5e6;        % 采样率 Hz
Nr             = 256;        % 每个chirp采样点
Nd             = 128;        % 每帧chirp数
Nrx            = 4;          % 接收天线数
fc             = 77.94745e9;       % 载频
c              = 3e8;
lambda         = c/fc;
d              = lambda/2;   % 阵元间距（默认 λ/2）
Tc             = 70e-6;      % chirp周期(含idle)
Q              = 64;         % 角度向FFT点数(零填充可> Nrx)

%% 读取 4RX 复数数据：I1..I4, Q1..Q4 (int16)
fname = './reflect_verticlal.bin';
fid = fopen(fname,'rb'); assert(fid>0,'无法打开数据文件');
raw = fread(fid, Nr*Nd*Nrx*2, 'int16=>int16'); fclose(fid); % *2 for I/Q
raw = reshape(raw, 8, []);                        % 8行 = 4I + 4Q
adc = double(raw(1:4,:)) + 1i*double(raw(5:8,:)); % 4 x (Nr*Nd)

% 组成立方体 [range_sample x chirp x rx]
cube = zeros(Nr, Nd, Nrx);
for rx = 1:Nrx
    cube(:,:,rx) = reshape(adc(rx,1:Nr*Nd), Nr, Nd);
end

%% —— 1D：距离向 FFT（沿第1维）
range_win = hanning(Nr);
cube_r = bsxfun(@times, cube, range_win);   % 加窗
RFFT = fft(cube_r, Nr, 1);                  % Nr x Nd x Nrx

%% —— 2D：多普勒向 FFT（沿第2维）
RFFT = RFFT - mean(RFFT, 2);                % 去静杂波
dop_win = hanning(Nd).';
RFFT_w = bsxfun(@times, RFFT, reshape(dop_win,[1 Nd 1]));
RD = fftshift(fft(RFFT_w, Nd, 2), 2);       % Nr x Nd x Nrx

%% —— 3D：角度向 FFT（沿第3维）
ant_win = ones(Nrx,1);                      % 可改 hanning(Nrx)
RDw = bsxfun(@times, RD, reshape(ant_win,[1 1 Nrx]));
RDA = fftshift(fft(RDw, Q, 3), 3);          % Nr x Nd x Q
RDA_dB = 20*log10(abs(RDA)+eps);

%% 物理坐标轴
range_axis = (0:Nr-1) * SampleRate * c / (2*FrequencySlope*Nr); % 距离
fD_axis    = (-Nd/2:Nd/2-1)/Nd * (1/Tc);                        % 多普勒频率
v_axis     = (lambda/2) * fD_axis;                              % 速度
u          = (-Q/2:Q/2-1)/Q;                                    % 归一化空间频率
angle_axis = asind((lambda/d)*u);                               % 角度(deg)

%% —— 1D FFT 的 3D 瀑布图（Range–Time）
RT_mag = sqrt(sum(abs(RFFT).^2, 3));        % 非相干合成4路
RT_dB  = 20*log10(RT_mag + eps);
slowtime_axis = (0:Nd-1)*Tc;
figure;
mesh(range_axis, slowtime_axis, RT_dB.');    % X=距离, Y=慢时间
xlabel('距离 (m)'); ylabel('慢时间 (s)'); zlabel('幅度 (dB)');
title('1D FFT 后的三维瀑布图（Range–Time）');
xlim([0, range_axis(end)]); grid on; box on; view([-40 30]);

%% —— 2D FFT 的 3D 图（Range–Doppler）
RD_mag = sqrt(sum(abs(RD).^2, 3));          % 非相干合成4路
RD_dB  = 20*log10(RD_mag + eps);
figure;
mesh(v_axis, range_axis, RD_dB);             % X=速度, Y=距离
xlabel('速度 (m/s)'); ylabel('距离 (m)'); zlabel('幅度 (dB)');
title('2D FFT 后的三维图（Range–Doppler）');
xlim([min(v_axis), max(v_axis)]); ylim([0, range_axis(end)]);
grid on; box on; view([-40 30]);

%% —— 3D FFT 的常见切片
% A) Range–Doppler @ ~0°
abin0 = Q/2 + 1;
figure; imagesc(v_axis, range_axis, RDA_dB(:,:,abin0)); axis xy;
xlabel('速度 (m/s)'); ylabel('距离 (m)'); title('Range–Doppler (≈0° 切片)');
cb=colorbar; cb.Label.String='幅度 (dB)';

% B) Range–Angle @ fd≈0
dbin0 = Nd/2 + 1;
figure; imagesc(angle_axis, range_axis, squeeze(RDA_dB(:,dbin0,:))); axis xy;
xlabel('角度 (deg)'); ylabel('距离 (m)'); title('Range–Angle (fd≈0 切片)');
cb=colorbar; cb.Label.String='幅度 (dB)';

% C) Doppler–Angle @ 选定距离门（自动选最强）
[~, rbin_max] = max(max(max(RDA_dB,[],3),[],2));
rbin = rbin_max;
figure; imagesc(angle_axis, v_axis, squeeze(RDA_dB(rbin,:,:)).'); axis xy;
xlabel('角度 (deg)'); ylabel('速度 (m/s)');
title(sprintf('Doppler–Angle (距离门 %.2f m)', range_axis(rbin)));
cb=colorbar; cb.Label.String='幅度 (dB)';

%% 分辨率/不模糊速度参考
v_res = lambda/(2*Nd*Tc);
v_max = lambda/(4*Tc);
fprintf('速度分辨率: %.3f m/s, 最大不模糊速度: ±%.3f m/s\n', v_res, v_max);

%% 峰值目标（距离/速度/角度估计）
[pk, idx] = max(RDA_dB(:));
[rbin, dbin, abin] = ind2sub(size(RDA_dB), idx);
R_est = range_axis(rbin);
fD    = (dbin - (Nd/2+1))/Nd * (1/Tc);
v_est = (lambda/2) * fD;
u_est = (abin - (Q/2+1))/Q;
ang_est = asind((lambda/d)*u_est);
fprintf('峰值目标：R=%.2f m, v=%.2f m/s, angle=%.1f°, 幅度=%.1f dB\n', ...
        R_est, v_est, ang_est, pk);
