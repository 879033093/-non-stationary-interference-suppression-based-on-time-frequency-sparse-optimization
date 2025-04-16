close all;
clc;
clear all;

%% Radar Parameters
% %% Radar Parameters
% Parameter settings
f0 = 77e8; % Starting frequency
numADCSamples = 200; % Number of ADC samples
Doppler_Number = 64; % Number of Doppler channels
freqSlopeRate = 70.006e12; % Frequency slope rate, units: MHz/us
adcSampleRate = 4e6; % Fast-time ADC sample rate, units: Ksps
% numChirps = 128; % Number of chirps per frame
Ts = 50e-3; % Time between each frame
   % vc = 3e8; % Speed of light, units: m/s
   % lambda = vc/f0; % Wavelength, calculated for 77GHz, units: m
   % % dmax = (adcSampleRate * 1e3) * vc / (2 * freqSlopeRate * 1e12); % Maximum unambiguous distance (parameters converted to SI units)
 % d = lambda / 2; % Element spacing
  % NofPules = 1;
    % xunni_annate = 2*4;
  % Q = 180;
 % Tc = 64e-6; % Chirp total period
 % %% Algorithm Parameters
 % rangeFFTNum = 256; % Number of FFT points for range dimension
 % useFramesNum = 1200; % Number of frames to calculate once
Tc = numADCSamples / adcSampleRate;
B = Tc * freqSlopeRate;
deltaR = 3e8 / (2 * B);

%% Parameter Settings
n_chirps = 64;   % Number of pulses per frame
n_samples = 200;  % Number of samples per pulse
N = 200;    % Number of FFT points in range direction
M = 64;        % Number of FFT points in Doppler direction
fs = 4e6;    % Sampling rate
c = 3.0e8;     % Speed of light
f0 = 77e9;      % Initial frequency
lambda = c / f0;      % Radar signal wavelength
d = lambda / 2;    % Antenna array spacing
Tc = 57e-6;   % Single chirp period
Tf = 50e-3;   % Frame period
B = 3990.34e6; % Signal bandwidth in Hz
S = B / (Tc - (7e-6));
Bv = (N / fs) * S;
rangeRes = c / 2 / Bv; % Range resolution
% Ts = 0:Tc:127*Tc; % Fast time axis
Ta = 0:Tf:1024*Tf; % Total time axis
Rr = 0:rangeRes:(N - 1) * rangeRes; % Range axis
Vr = (-M / 2:M / 2 - 1) * lambda / Tc / M / 2; % Velocity axis

%% File Reading and Data Rearrangement

numADCSamples = 200; % number of ADC samples per chirp
numADCBits = 16;         % number of ADC bits per sample
numRX = 8;        % number of receivers
numTx = 2;
numLanes = 2;       % do not change. Number of lanes is always 2
isReal = 0;             % set to 1 if real-only data, 0 if complex data
chirpLoop = 2;

Xintiao_Radar = zeros(1, 10);
Xintiao_Jiange_Radar = zeros(1, 10);
Xintiao_Mov = zeros(1, 10);
Xintiao_Jiange_Mov = zeros(1, 10);

%% Read Bin Files
dataFolder = 'E:\1.Graduate\Breathing and Heartbeat Data 2025.01.08\0.5 Person 4';
videoFolder = 'D:\1.Graduate\Personal Data Measurement\Paper Replication\Algorithm for Transient Interference Suppression in Sky-wave Over-the-Horizon Radar Based on Time-Frequency Sparse Priors\Video\8';
for iii = 1  % Data ID, you can also process multiple data at once

    binFilename = fullfile(dataFolder, sprintf('%d.bin', iii));
    fid = fopen(binFilename, 'r');
    adcDataRow = fread(fid, 'int16');
    if numADCBits ~= 16
        l_max = 2^(numADCBits - 1) - 1;
        adcDataRow(adcDataRow > l_max) = adcDataRow(adcDataRow > l_max) - 2^numADCBits;
    end
    fclose(fid);

    fileSize = size(adcDataRow, 1);
    PRTnum = fix(fileSize / (numADCSamples * numRX));
    fileSize = PRTnum * numADCSamples * numRX;
    adcData = adcDataRow(1:fileSize);
    
    % Real data reshape, filesize = numADCSamples * numChirps
    if isReal
        numChirps = fileSize / numADCSamples / numRX;
        LVDS = zeros(1, fileSize);
             % Create column for each chirp
        LVDS = reshape(adcData, numADCSamples * numRX, numChirps);
              % Each row is data from one chirp
        LVDS = LVDS.';
    else
        numChirps = fileSize / 2 / numADCSamples / numRX; % Contains real and imaginary parts, divided by 2
        LVDS = zeros(1, fileSize / 2);
             % Combine real and imaginary parts into complex data
           % Read in file: 2I is followed by 2Q
        counter = 1;
        for i = 1:4:fileSize - 1
            LVDS(1, counter) = adcData(i) + sqrt(-1) * adcData(i + 2);
            LVDS(1, counter + 1) = adcData(i + 1) + sqrt(-1) * adcData(i + 3);
            counter = counter + 2;
        end
        % Create column for each chirp
        LVDS = reshape(LVDS, numADCSamples * numRX, numChirps);
        % Each row is data from one chirp
        LVDS = LVDS.';
    end
end

%% Data Rearrangement
numframe = 1024;
adcData = zeros(numRX, numChirps * numADCSamples);
for row = 1:numRX
    for i = 1:numChirps
        adcData(row, (i - 1) * numADCSamples + 1:i * numADCSamples) = LVDS(i, (row - 1) * numADCSamples + 1:row * numADCSamples);
    end
end

adcdata = adcData;

n_frame = floor(length(adcdata) / n_chirps / n_samples); % Total number of frames, result is 256
% length is used to find the length of the array
      % arr = [2.4, 5.7, 1.9]; result = floor(arr); % Result will be [2, 5, 1]

retVal = reshape(adcData(tianxian, :), numADCSamples, numChirps); % Get data from the first receiver antenna
retVal = reshape(retVal, [200, n_chirps, n_frame]); % Now retVal is still data from one antenna, but the arrangement is 200*64*1024 (samples*chirps*frames)

process_adc = zeros(numADCSamples, n_frame); % For each frame, take the first of two chirps, 200*1024

process_adc = squeeze(sum(retVal, 2)); % Non-coherent summation of all chirps in retVal
channel_data(:,:,tianxian) = process_adc';
end

   %% Import Data
   % load channel_data  %1200x200x8 (take one chirp per frame) frame samples
rangeFFTNum = 256;
useFramesNum = 1024;

%% Data Processing
loop_cnt = floor(length(channel_data(:, 1)) / useFramesNum);

k = 1;
    
use_channel_data = channel_data((k - 1) * useFramesNum + 1 : k * useFramesNum, :, :);
  % use_channel_data is the raw data of 1024*200*8
  %% Mean Cancellation and Pulse Compression
rangeProfile = MTI_PulseCompression(use_channel_data, 0, rangeFFTNum); % 1200x256x8 complex double
sum_rangProfile = sum(abs(rangeProfile(:,:,1)), 1); % sum_rangProfile size is 1*256, equivalent to data from one chirp
[~, targetIndex] = max(sum_rangProfile);  % Target distance

%% DOA Estimation
searchAngleRange = 60; % Angle search range ±searchAngleRange degrees

[~, azimuSpectrogram, Rxv] = IWR1642ODS_DOA(rangeProfile, 2, useFramesNum, searchAngleRange);
% Rxv = pinv(Rx) calculates the pseudo-inverse of the covariance matrix and computes the angle spectrum azimuSpectrogram (121*256)

%% Find Targets at Different Angles
maxAzimu = max(azimuSpectrogram, [], 2); % 121x1, max(A, [], 2) returns the maximum value for each row (maximum value for each angle)
  % [values,peaks_index] = findpeaks(maxAzimu,'minpeakheight',2e10); % findpeaks(maxAzimu);
       % The 'minpeakheight' parameter sets a threshold; only peaks above 3,000,000 are considered valid targets.
    %    Assume maxAzimu is a column vector of size 121x1
     % Use findpeaks to find all peaks
[allPeaks, allPeakIndices] = findpeaks(maxAzimu);

% Define a condition to filter out valid peaks
validPeaks = [];          % Stores indices of valid peaks
validValues = [];  % Stores values of valid peaks

for i = 1:length(allPeakIndices)
    peakIdx = allPeakIndices(i);

          % Check if the left and right four neighbors of the peak are smaller
    if peakIdx > 4 && peakIdx < length(maxAzimu) - 4
        leftNeighbors = maxAzimu(peakIdx-4:peakIdx-1);
        rightNeighbors = maxAzimu(peakIdx+1:peakIdx+4);

        % If all neighbors are smaller, the peak is valid
        if all(leftNeighbors < maxAzimu(peakIdx)) && all(rightNeighbors < maxAzimu(peakIdx))
            validPeaks = [validPeaks; peakIdx];
            validValues = [validValues; maxAzimu(peakIdx)];
        end
    end
end

     % Get the top one or two valid peaks
if length(validPeaks) >= 1
    % Sort valid peak values and select the top one
    [sortedValues, sortedIndices] = sort(validValues, 'descend');
    top2Values = sortedValues(1:1);
    top2Indices = validPeaks(sortedIndices(1:1));

       % Store results in values and peaks_index
    values = top2Values;
    peaks_index = top2Indices;
else
    disp('Not enough valid peaks found');
    values = [];
    peaks_index = [];
end

%% Plotting Range-Angle Map
figure(1);  
plot(1:length(maxAzimu),10*log10(maxAzimu));    
hold on; grid on
plot(peaks_index,10*log10(maxAzimu(peaks_index)),'bd');
xlabel('Angle (°)','FontSize',16);ylabel('Gain (dB)','FontSize',16);
       % title('MVDR Angle Estimation');
hold off

figure(2);
imagesc(-searchAngleRange:searchAngleRange,(1:rangeFFTNum)*deltaR,abs(azimuSpectrogram.'));
ylabel('Range (m)');xlabel('Angle (°)');        % title('Range-Angle Spectrogram');
ylim([0.5,2])
axis xy

%% Select the top two angles with highest amplitudes 
   % (Why not spatial fusion? Because the carotid artery may be selected, which has phase difference. See "blood pressure" paper for details.)
[values_sort,index] = sort(values,'descend');  % Sort values in descending order
 peaks_index_sort = peaks_index(index);
 peaks_index_max = peaks_index_sort(1:1);

 distance_indices = zeros(1, length(peaks_index_max)); % Stores distance index for each angle
distance_values = zeros(1, length(peaks_index_max));  % Stores corresponding maximum values

for dd = 1:length(peaks_index_max)
    [distance_values(dd), distance_indices(dd)] = max(azimuSpectrogram(peaks_index_max(dd), :));
end

% Explanation:
  % peaks_index: positions of peaks in maxAzimu, i.e., angle indices.
 % peaks_index(index): reorders original indices using sorting index.
 % peaks_index_sort is the angle index array sorted by peak value.
 % E.g., if peaks_index = [10, 20, 30], and index = [3, 1, 2], then peaks_index_sort = [30, 10, 20].

 %% Extract breathing and heartbeat from different targets
target_rangeProfile = zeros(useFramesNum,length(peaks_index_max));
 % Store the range profile data of the targets after angle filtering.
  % Each column corresponds to one target's data; typically there are only two targets.

    % rangeProfile is the result after FFT, with size 1200x256x8

for m = 1:length(peaks_index_max)

    detAngle = -searchAngleRange + peaks_index_max(m) * (searchAngleRange * 2 / length(azimuSpectrogram(:,1)));
    dashabi111(:,m) = detAngle;
    % detAngle = -60° + target angle index × (2×60° / 121)
    % Calculates angle of arrival for each target

    fai = 2 * pi * sin(detAngle / 180 * pi) * d / lambda;
    dashabi222(:,m) = fai;
    aTheta = [1,exp(-1j*1*fai),exp(-1j*2*fai),exp(-1j*3*fai),exp(-1j*4*fai),exp(-1j*5*fai),exp(-1j*6*fai),exp(-1j*7*fai)].';

    Wopt = (Rxv * aTheta) / (aTheta' * Rxv * aTheta);   
    % Classical beamforming to enhance target signal and suppress others
    % detAngle gives the angle of the target, used to compute fai and aTheta

    xt = squeeze(rangeProfile(:,distance_indices(:,m),:)); 
    % Extract 1200x8 data for each target
    target_rangeProfile(:,m) = xt * Wopt;  % 1200x2 result, one column per target
end

%% Estimate breathing and heartbeat for each target and plot
slowtime_fs = 1 / Ts;
target_num = length(target_rangeProfile(1,:));
 breathRate = zeros(1,target_num);
 heartRate  = zeros(1,target_num);
useFramesNum = 1024;

for kk = 1:target_num
     %% Extract breathing and heartbeat waveform from range profile
     target_profile = target_rangeProfile(:,kk);
     dcRemove_ag = angle(target_profile);  % Get phase
     unwrap_dcRemove_ag = unwrap(dcRemove_ag);  % Phase unwrapping
     diff_ag = unwrap_dcRemove_ag(1:end-1) - unwrap_dcRemove_ag(2:end);  % Phase difference
       agl = diff_ag;  % A 1199x1 double array of phase data

    %% Band-pass filtering to extract breathing and heartbeat signals
     agl = agl - mean(agl);  % Remove DC component

    % --- Extract breathing waveform and rate
    breath_wave = filter(RR_BPF20, agl);  % Breathing component
    % [b,a]=butter(5,[0.01,0.05]); % 5th-order Butterworth filter for 0.1-0.5Hz (breathing)
    % breath_wave = filter(b,a,agl);

    breath_fft = abs(fftshift(fft(breath_wave)));
     breath_fft = breath_fft(ceil(end/2):end);  % Use only one-sided spectrum
     [~,breath_index] = max(breath_fft);
     breath_hz = breath_index * (slowtime_fs / 2 / length(breath_fft));
    breathRate(kk) = ceil(breath_hz * 60);

    % --- Extract heartbeat waveform and rate
     heart_wave = filter(HR_BPF20, agl);  % Heartbeat component
     % [b,a]=butter(9,[0.08,0.2]); % 9th-order Butterworth filter for 0.8-2Hz (heartbeat)
    % heart_wave = filter(b,a,agl);

    heart_fft = abs(fftshift(fft(heart_wave)));
     heart_fft = heart_fft(ceil(end/2):end);
     [~,heart_index] = max(heart_fft);
     heart_hz = heart_index * (slowtime_fs / 2 / length(heart_fft));
    heartRate(kk) = ceil(heart_hz * 60);

    % --- Display waveforms and spectra
    xAxis = linspace(0,slowtime_fs/2,length(breath_fft)) * 60;
     t_axis = linspace(0,1/slowtime_fs * length(target_profile),length(target_profile)-1 );

    figure;
    subplot(141); plot(t_axis,breath_wave); title('Breathing Waveform'); xlabel('Time (s)'); ylabel('Amplitude'); grid on
    subplot(143); plot(xAxis,breath_fft); title('Breathing Spectrum'); xlabel('Frequency (bpm)'); ylabel('Amplitude'); grid on
    subplot(142); plot(t_axis,heart_wave); title('Heartbeat Waveform'); xlabel('Time (s)'); ylabel('Amplitude'); grid on
    subplot(144); plot(xAxis,heart_fft); title('Heartbeat Spectrum'); xlabel('Frequency (bpm)'); ylabel('Amplitude'); grid on
end

heart_data=agl;
fs=20;
num_chirps=1024;


% Transient Interference Removal

% Initialize parameters
fs = 20;          % Sampling rate (Hz)
T = 51.2;         % Duration (s)
N = fs * T - 1;   % Total number of samples
gamma = 0.1;      % Regularization parameter gamma
gamma1 = 0.1;     % Regularization parameter gamma1
rho = 0.1;        % Step size parameter
K = 100;          % Maximum number of iterations
epsilon = 1e-6;   % Convergence threshold

% Simulate breathing and heartbeat signal
t = (0:N-1)' / fs; 
breathing_signal = sin(2 * pi * 0.25 * t);    % Breathing component (0.25 Hz)
heartbeat_signal = 0.5 * sin(2 * pi * 1.2 * t); % Heartbeat component (1.2 Hz)
clean_signal = breathing_signal + heartbeat_signal;

% Add transient interference
transient_interference = zeros(N, 1);
transient_interference(100:120) = 2 * randn(21, 1);     % First burst
transient_interference(600:620) = 1.5 * randn(21, 1);   % Second burst

% Add noise
noise = 0.1 * randn(N, 1);

% Construct signal with interference
raw_signal = heart_data;

% Initialize iterative variables
o = zeros(N, 1);       % Transient interference component
x = zeros(N, 1);       % Clean frequency spectrum
n_var = zeros(N, 1);   % Auxiliary variable n
b = zeros(N, 1);       % Lagrange multiplier
u = zeros(N, 1);       % Intermediate variable in FFT domain
k = 1;                 % Initial iteration

% Iterative process
while k <= K
    % Gradient computation
    grad_g = o - raw_signal + sqrt(N) * ifft(x + n_var - b / rho);

    % Update transient interference component o
    step_size = 0.5;
    o = soft_threshold(o - step_size * grad_g, gamma * step_size / rho);

    % Update intermediate variable u
    u = fft(raw_signal - o) / sqrt(N);

    % Update clean frequency spectrum x
    x = soft_threshold(u - n_var + b / rho, gamma1 / rho);

    % Update auxiliary variable n_var
    n_var = (b + rho * (u - x)) / (rho + 1);

    % Update Lagrange multiplier
    b = b + rho * (u - x - n_var);

    % Convergence check
    if norm(x - fft(raw_signal - o) / sqrt(N))^2 / norm(x)^2 < epsilon
        break;
    end

    % Update iteration count
    k = k + 1;
end

% Result
pure_signal = ifft(x, 'symmetric');  % Cleaned signal

% Plot results
figure;
subplot(4, 1, 1);
plot(t, raw_signal); title('Signal with Transient Interference');
subplot(4, 1, 2);
plot(t, o); title('Extracted Transient Interference');
subplot(4, 1, 3);
plot(t, clean_signal); title('Ideal Breathing and Heartbeat Signal');
subplot(4, 1, 4);
plot(t, pure_signal); title('Cleaned Signal');

% Define soft-thresholding function
% Heartbeat frequency domain analysis
heart_data = pure_signal;

figure;
plot(t_axis, heart_data);
xlabel('Time (s)');
ylabel('Amplitude');

t = 0:51.2/1024:51.2-51.2/1024;

N1 = length(heart_data);
fshift = (-N1/2:N1/2-1)*(fs/N1); % Zero-centered frequency
f = (0:N1-1)*(fs/N1);            % Frequency at each point
heart_fre = abs(fftshift(fft(heart_data)));
heart = abs(fft(heart_data));

figure;
plot(f(1:200), heart(1:200));
xlabel('F (Hz)');
ylabel('Amplitude');
xlim([0.0 2.5])

figure;
plot(f(1:1000), heart(1:1000));
xlabel('F (Hz)');
ylabel('Amplitude');
xlim([0 2])

% Fundamental Heart Rate Frequency Extraction
% Parameter settings
fs = 20;
N = length(heart_data);
f = (0:N-1)*(fs/N);
fshift = (-N/2:N/2-1)*(fs/N);

% Spectrum computation
heart_fre = abs(fftshift(fft(heart_data)));

% Heart rate range [0.78Hz, 2Hz], divided into n sub-bands
freq_min = 0.78;
freq_max = 2.0;
n = 70;
freq_bins = linspace(freq_min, freq_max, n+1);
bin_width = freq_bins(2) - freq_bins(1);

% Harmonic energy accumulation
energy = zeros(1, n);
for i = 1:n
    f_center = (freq_bins(i) + freq_bins(i+1)) / 2;

    first_harmonic_idx = find(abs(fshift - f_center) < bin_width/2);
    second_harmonic_idx = find(abs(fshift - 2*f_center) < bin_width/2);
    third_harmonic_idx = find(abs(fshift - 3*f_center) < bin_width/2);
    fourth_harmonic_idx = find(abs(fshift - 4*f_center) < bin_width/2);

    energy(i) = sum(heart_fre(first_harmonic_idx).^2) + ...
                sum(heart_fre(second_harmonic_idx).^2) + ...
                sum(heart_fre(third_harmonic_idx).^2) + ...
                sum(heart_fre(fourth_harmonic_idx).^2);
end

% Determine fundamental frequency with maximum energy
[~, max_idx] = max(energy);
heart_rate_freq = (freq_bins(max_idx) + freq_bins(max_idx+1)) / 2;

disp(['Heart rate fundamental frequency: ', num2str(heart_rate_freq), ' Hz']);

% Comb filter with roll-off characteristics using FIR filters
% Parameter settings
fs = 20;
N = length(heart_data);
t = (0:N-1) / fs;

% Fundamental frequency and harmonics
num_harmonics = 4;
bandwidth = 0.2;

% Filter design
roll_off = 0.01;
filtered_signal = zeros(size(heart_data));

% Design bandpass filters for each harmonic
for k = 1:num_harmonics
    f_center = k * heart_rate_freq;
    f_low = max(0, f_center - bandwidth / 2);
    f_high = min(fs / 2, f_center + bandwidth / 2);
    
    freq_low = f_low / (fs / 2);
    freq_high = f_high / (fs / 2);
    
    fir_coeff = fir1(128, [freq_low, freq_high], hamming(129));
    filtered_signal = filtered_signal + filter(fir_coeff, 1, heart_data);
end

% Plot analysis
figure;
subplot(3, 1, 1);
plot(t, heart_data);
title('Original Heartbeat Signal (with breathing interference)');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(3, 1, 2);
[H, f] = freqz(fir_coeff, 1, 1024, fs);
plot(f, abs(H));
title('FIR Bandpass Filter Frequency Response for Single Harmonic');
xlabel('Frequency (Hz)');
ylabel('Amplitude');

subplot(3, 1, 3);
plot(t, filtered_signal);
title('Filtered Heartbeat Signal');
xlabel('Time (s)');
ylabel('Amplitude');

% Peak detection from filtered heartbeat signal
[peaks, locs] = findpeaks(filtered_signal, 'MinPeakHeight', mean(filtered_signal), 'MinPeakDistance', 10);

% Plot peaks
figure;
subplot(3, 1, 1);
plot(t, filtered_signal, 'b'); hold on;
plot(t(locs), peaks, 'ro');
title('Filtered Heartbeat Signal with Peaks');
xlabel('Time (s)');
ylabel('Amplitude');
legend('Heartbeat Signal', 'Peaks');
xlim([10 50]);
grid on;



