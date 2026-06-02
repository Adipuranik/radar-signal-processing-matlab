clc; clear; close all;

%% === RADAR SYSTEM PARAMETERS ===
c        = 3e8;          % Speed of light (m/s)
fc       = 10e9;         % Carrier frequency: 10 GHz (X-band radar)
lambda   = c / fc;       % Wavelength

B        = 5e6;          % Chirp bandwidth (Hz)
T_pulse  = 50e-6;        % Pulse duration (50 microseconds) — FIXED
fs       = 20e6;         % Sampling frequency (Hz)
PRF      = 1000;         % Pulse Repetition Frequency (Hz)
N_pulses = 64;           % Number of pulses (for Doppler processing)

% Target parameters
R_target = 750;          % Target range: 750 m — FIXED
v_target = 50;           % Target velocity: 50 m/s (toward radar)
RCS      = 1;            % Radar Cross Section (m^2)

fprintf('=== Radar System Setup ===\n');
fprintf('Carrier Frequency : %.1f GHz\n', fc/1e9);
fprintf('Target Range      : %.0f m\n',   R_target);
fprintf('Target Velocity   : %.0f m/s\n', v_target);
fprintf('Wavelength        : %.4f m\n',   lambda);


%% === GENERATE LFM CHIRP PULSE ===
t_pulse = 0 : 1/fs : T_pulse - 1/fs;
k       = B / T_pulse;              % Chirp rate (Hz/s)

% Transmitted chirp
tx = exp(1j * pi * k * t_pulse.^2);

% Time delay for echo from target
t_delay = 2 * R_target / c;
n_delay = round(t_delay * fs);

% Doppler frequency shift
f_doppler = 2 * v_target * fc / c;
fprintf('\nDoppler Frequency Shift: %.2f Hz\n', f_doppler);

N_samples = length(t_pulse);

% Build received signal
rx_clean = zeros(1, N_samples);
if n_delay < N_samples
    rx_clean(n_delay+1:end) = tx(1:end-n_delay) .* ...
        exp(1j * 2*pi * f_doppler * t_pulse(1:end-n_delay));
end

% Add complex Gaussian noise (SNR = 20 dB) — increased for clear detection
SNR_dB    = 20;
sig_pow   = mean(abs(tx).^2);
noise_pow = sig_pow / 10^(SNR_dB/10);
noise     = sqrt(noise_pow/2) * (randn(size(rx_clean)) + 1j*randn(size(rx_clean)));
rx        = rx_clean + noise;

% Plot transmitted vs received
figure('Name','Chirp Signal','Color','k');
subplot(2,1,1);
plot(t_pulse*1e6, real(tx), 'c');
title('Transmitted LFM Chirp','Color','w');
xlabel('Time (\mus)','Color','w'); ylabel('Amplitude','Color','w');
set(gca,'Color','k','XColor','w','YColor','w');

subplot(2,1,2);
plot(t_pulse*1e6, real(rx), 'y');
title('Received Signal (with noise + delay)','Color','w');
xlabel('Time (\mus)','Color','w'); ylabel('Amplitude','Color','w');
set(gca,'Color','k','XColor','w','YColor','w');


%% === MATCHED FILTER (Pulse Compression) ===
mf_kernel  = conj(fliplr(tx));
compressed = conv(rx, mf_kernel, 'same');

% Range axis
N_range    = length(compressed);
range_axis = (0:N_range-1) * c / (2*fs);  % in metres

figure('Name','Matched Filter','Color','k');
subplot(2,1,1);
plot(range_axis/1000, abs(rx), 'y');
title('Before Matched Filter (Raw Received)','Color','w');
xlabel('Range (km)','Color','w'); ylabel('Magnitude','Color','w');
set(gca,'Color','k','XColor','w','YColor','w');

subplot(2,1,2);
plot(range_axis/1000, abs(compressed), 'c');
title('After Matched Filter (Pulse Compressed)','Color','w');
xlabel('Range (km)','Color','w'); ylabel('Magnitude','Color','w');
xline(R_target/1000, 'r--', sprintf('Target %.0fm', R_target), ...
    'LabelVerticalAlignment','bottom','Color','r');
set(gca,'Color','k','XColor','w','YColor','w');

[~, peak_idx] = max(abs(compressed));
fprintf('\n=== Matched Filter Applied ===\n');
fprintf('Peak detected at range: %.0f m (true target: %.0f m)\n', ...
    range_axis(peak_idx), R_target);


%% === RANGE-DOPPLER MAP ===
slow_time   = (0:N_pulses-1) / PRF;
data_matrix = zeros(N_pulses, N_samples);

for p = 1:N_pulses
    phase_p = exp(1j * 2*pi * f_doppler * slow_time(p));

    rx_p = zeros(1, N_samples);
    if n_delay < N_samples
        rx_p(n_delay+1:end) = tx(1:end-n_delay) * phase_p;
    end

    noise_p      = sqrt(noise_pow/2) * (randn(1,N_samples) + 1j*randn(1,N_samples));
    rx_p         = rx_p + noise_p;
    comp_p       = conv(rx_p, mf_kernel, 'same');
    data_matrix(p,:) = comp_p;
end

% Doppler FFT along slow-time axis
rd_map        = fftshift(fft(data_matrix, [], 1), 1);
velocity_axis = linspace(-PRF/2, PRF/2, N_pulses) * lambda / 2;
range_axis_km = (0:N_samples-1) * c / (2*fs*1000);

figure('Name','Range-Doppler Map','Color','k');
imagesc(range_axis_km, velocity_axis, 20*log10(abs(rd_map) + eps));
colormap jet; colorbar;
title('Range-Doppler Map','Color','w');
xlabel('Range (km)','Color','w'); ylabel('Velocity (m/s)','Color','w');
clim([-20 60]);
hold on;
xline(R_target/1000, 'w--', 'LineWidth', 1.5);
yline(v_target,       'w--', 'LineWidth', 1.5);
set(gca,'Color','k','XColor','w','YColor','w');
fprintf('\nRange-Doppler Map generated with %d pulses.\n', N_pulses);


%% === CA-CFAR DETECTION ===
range_profile  = abs(compressed);
guard_cells    = 2;
training_cells = 8;
Pfa            = 1e-3;
alpha_cfar     = training_cells * (Pfa^(-1/training_cells) - 1);

cfar_threshold = zeros(size(range_profile));
detections     = zeros(size(range_profile));
half_win       = guard_cells + training_cells;

for i = half_win+1 : length(range_profile)-half_win
    train_left  = range_profile(i-half_win : i-guard_cells-1);
    train_right = range_profile(i+guard_cells+1 : i+half_win);
    noise_est   = mean([train_left, train_right]);

    cfar_threshold(i) = alpha_cfar * noise_est;
    if range_profile(i) > cfar_threshold(i)
        detections(i) = range_profile(i);
    end
end

figure('Name','CFAR Detection','Color','k');
plot(range_axis/1000, range_profile, 'b', 'LineWidth', 1.2); hold on;
plot(range_axis/1000, cfar_threshold, 'r--', 'LineWidth', 1.2);
det_idx = detections > 0;
if any(det_idx)
    stem(range_axis(det_idx)/1000, detections(det_idx), ...
        'g', 'filled', 'MarkerSize', 6);
end
xline(R_target/1000, 'w:', 'LineWidth', 1);
title('CA-CFAR Target Detection','Color','w');
xlabel('Range (km)','Color','w'); ylabel('Magnitude','Color','w');
legend('Range Profile','CFAR Threshold','Detections','True Target', ...
    'TextColor','w','Color','k');
set(gca,'Color','k','XColor','w','YColor','w');

det_ranges = range_axis(det_idx);
fprintf('\n=== CFAR Detection Results ===\n');
if ~isempty(det_ranges)
    fprintf('Target detected at: %.0f m  (true: %.0f m)\n', det_ranges(1), R_target);
else
    fprintf('No detections — try increasing SNR_dB value.\n');
end
