% 高通滤波器 MATLAB 仿真
% 截止频率 fc = 30 kHz, Q = 0.5, 双电源 ±15 V, 单位增益 Sallen-Key

clear; clc; close all;

%% 滤波器参数
fc = 30e3;          % 30 kHz
R = 2.41e3;         % 2.41 kΩ
C = 2.2e-9;         % 2.2 nF
wc = 1/(R*C);       % 角频率
Q = 0.5;

% 传递函数: H(s) = s^2 / (s^2 + (wc/Q)*s + wc^2)
num = [1 0 0];
den = [1, wc/Q, wc^2];

% 设计数字滤波器用于零状态滤波（双线性变换，保持截止频率准确）
Fs = 1e6;           % 采样频率 1 MHz
[b, a] = bilinear(num, den, Fs);   % 离散化

% 仿真时长：为了瞬态衰减，总时长取 2 ms，其中后 1 ms 作为稳态
t_total = 0:1/Fs:2e-3;    % 0 到 2 ms
N_total = length(t_total);
f_res = Fs/N_total;       % 频率分辨率 500 Hz

% 提取稳态的时间索引（舍弃前 1 ms，瞬态应已充分衰减）
steady_start = 1e-3;      % 稳态开始时间
idx_steady = t_total >= steady_start;
t_steady = t_total(idx_steady);
N_steady = length(t_steady);

%% ======== 信号 1：40 kHz 单音正弦 ========
f1 = 40e3;
x1_total = sin(2*pi*f1*t_total);
y1_total = filter(b, a, x1_total);    % 零状态滤波（包含瞬态）
x1 = x1_total(idx_steady);
y1 = y1_total(idx_steady);

plot_signals(t_steady, x1, y1, Fs, '单音正弦 40 kHz');
plot_corr_psd(x1, y1, Fs, '单音正弦 40 kHz');

%% ======== 信号 2：三音正弦 (15k + 50k + 100k) ========
f2 = [15e3, 50e3, 100e3];
x2_total = sin(2*pi*f2(1)*t_total) + sin(2*pi*f2(2)*t_total) + sin(2*pi*f2(3)*t_total);
y2_total = filter(b, a, x2_total);
x2 = x2_total(idx_steady);
y2 = y2_total(idx_steady);

plot_signals(t_steady, x2, y2, Fs, '三音正弦 15+50+100 kHz');
plot_corr_psd(x2, y2, Fs, '三音正弦 15+50+100 kHz');

%% ======== 信号 3：20 kHz 方波 ========
f3 = 20e3;
x3_total = square(2*pi*f3*t_total);   % ±1 方波
y3_total = filter(b, a, x3_total);
x3 = x3_total(idx_steady);
y3 = y3_total(idx_steady);

plot_signals(t_steady, x3, y3, Fs, '20 kHz 方波');
plot_corr_psd(x3, y3, Fs, '20 kHz 方波');

%% ======== 滤波器幅频特性 ========
figure;
[h, w] = freqz(b, a, 4096, Fs);
freqz_plot = w/1e3;
mag_dB = 20*log10(abs(h));
plot(freqz_plot, mag_dB, 'LineWidth',1.5); grid on;
xlabel('频率 (kHz)'); ylabel('幅度 (dB)');
title('数字滤波器幅频响应（逼近模拟原型）');
xline(fc/1e3, '--r', '截止频率 30 kHz');
ylim([-30 5]);

%% ******************** 自定义函数 ********************
function plot_signals(t, x, y, Fs, titlestr)
    % 绘制时域波形和单边幅度谱（峰值谱）
    figure;
    t_ms = t*1e3;
    
    subplot(2,2,1);
    plot(t_ms, x); xlabel('时间 (ms)'); ylabel('幅度');
    title(['输入时域 - ', titlestr]); grid on; xlim([0 0.5]);  % 显示前0.5 ms
    subplot(2,2,2);
    plot(t_ms, y); xlabel('时间 (ms)'); ylabel('幅度');
    title(['输出时域 - ', titlestr]); grid on; xlim([0 0.5]);
    
    N = length(x);
    X = fft(x)/N;
    Y = fft(y)/N;
    % 单边谱，乘以2（直流分量不乘2）
    X_oneside = abs(X(1:floor(N/2)+1))*2;
    Y_oneside = abs(Y(1:floor(N/2)+1))*2;
    X_oneside(1) = X_oneside(1)/2;   % 直流分量
    Y_oneside(1) = Y_oneside(1)/2;
    f_axis = Fs*(0:floor(N/2))/N / 1e3;  % kHz
    
    subplot(2,2,3);
    stem(f_axis, X_oneside, '.', 'MarkerSize', 4);
    xlabel('频率 (kHz)'); ylabel('幅度');
    title('输入频谱'); grid on; xlim([0 200]);
    subplot(2,2,4);
    stem(f_axis, Y_oneside, '.', 'MarkerSize', 4);
    xlabel('频率 (kHz)'); ylabel('幅度');
    title('输出频谱'); grid on; xlim([0 200]);
    sgtitle(['时域与频谱 - ', titlestr]);
end

function plot_corr_psd(x, y, Fs, titlestr)
    % 计算自相关和功率谱密度（修正窗长问题）
    [acfx, lags] = xcorr(x, 'biased');
    [acfy, ~] = xcorr(y, 'biased');
    lags_ms = lags/Fs * 1e3;  % ms
    
    % 根据数据长度自适应选择 pwelch 窗长
    N = length(x);
    winLen = min(512, N);   % 窗长不超过信号长度
    if winLen < 128
        winLen = N;          % 极短数据就用全数据做周期图（无平均）
    end
    noverlap = round(winLen * 0.5);
    nfft = max(1024, 2^nextpow2(winLen));
    
    % 功率谱密度（Welch法）
    [Pxx, f_psd] = pwelch(x, hamming(winLen), noverlap, nfft, Fs);
    [Pyy, ~] = pwelch(y, hamming(winLen), noverlap, nfft, Fs);
    
    figure;
    subplot(2,2,1);
    plot(lags_ms, acfx); 
    xlabel('滞后 (ms)'); ylabel('自相关');
    title(['输入自相关 - ', titlestr]); grid on;
    xlim([-0.15 0.15]);
    
    subplot(2,2,2);
    plot(lags_ms, acfy); 
    xlabel('滞后 (ms)'); ylabel('自相关');
    title(['输出自相关 - ', titlestr]); grid on;
    xlim([-0.15 0.15]);
    
    subplot(2,2,3);
    plot(f_psd/1e3, 10*log10(Pxx)); 
    xlabel('频率 (kHz)'); ylabel('PSD (dB/Hz)');
    title('输入功率谱密度'); grid on; xlim([0 200]);
    
    subplot(2,2,4);
    plot(f_psd/1e3, 10*log10(Pyy)); 
    xlabel('频率 (kHz)'); ylabel('PSD (dB/Hz)');
    title('输出功率谱密度'); grid on; xlim([0 200]);
    
    sgtitle(['自相关与功率谱密度 - ', titlestr]);
end