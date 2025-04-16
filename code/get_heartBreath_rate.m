function [breathRate, heartRate] = get_heartBreath_rate(target_rangeProfile, slowtime_fs)  
% slowtime_fs: sampling rate

% Calculate the number of targets
 target_num = length(target_rangeProfile(1,:));
 breathRate = zeros(1, target_num);
 heartRate  = zeros(1, target_num);
 useFramesNum = 1024;

for kk = 1:target_num
 %% Extract respiratory and heartbeat signals from the target column in the range profile
    target_profile = target_rangeProfile(:, kk);
       dcRemove_ag = angle(target_profile);              % Extract phase
      unwrap_dcRemove_ag = unwrap(dcRemove_ag);         % Phase unwrapping
      diff_ag = unwrap_dcRemove_ag(1:end-1) - unwrap_dcRemove_ag(2:end); % Phase difference
     agl = diff_ag;   % agl is a double vector of size 1199x1 representing the phase

  %% Use bandpass filters to isolate respiratory and heartbeat components
     agl = agl - mean(agl);    % Remove DC component
%---------------------------------Filter to obtain respiratory signal and rate
    breath_wave = filter(RR_BPF20, agl);              % Respiratory signal (1199x1 double)
    % [b, a] = butter(5, [0.01, 0.05]);               % Optional: 5th-order Butterworth filter for respiration
  % breath_wave = filter(b, a, agl);                % Apply optional filter

    % Compute respiratory rate
    breath_fft = abs(fftshift(fft(breath_wave)));
    breath_fft = breath_fft(ceil(end/2):end);      % Use only one side of the spectrum
    [~, breath_index] = max(breath_fft);
    breath_hz = breath_index * (slowtime_fs / 2 / length(breath_fft));
    breathRate(kk) = ceil(breath_hz * 60);    % Convert Hz to breaths per minute

    %---------------------------------Filter to obtain heartbeat signal and rate
      heart_wave = filter(HR_BPF20, agl);         % Heartbeat signal
    % [b, a] = butter(9, [0.08, 0.2]);              % Optional: 9th-order Butterworth filter for heartbeat
    % heart_wave = filter(b, a, agl);        % Apply optional filter

    % Compute heart rate
    heart_fft = abs(fftshift(fft(heart_wave)));
     heart_fft = heart_fft(ceil(end/2):end);      % Use only one side of the spectrum
    [~, heart_index] = max(heart_fft);
     heart_hz = heart_index * (slowtime_fs / 2 / length(heart_fft));
    heartRate(kk) = ceil(heart_hz * 60);             % Convert Hz to beats per minute

    %---------Visualization
    xAxis = linspace(0, slowtime_fs/2, length(breath_fft)) * 60;
    t_axis = linspace(0, 1/slowtime_fs * length(target_profile), length(target_profile) - 1);

    % figure;
    % subplot(311); plot([t_axis,60], dcRemove_ag); xlabel('Time (s)'); ylabel('Amplitude'); title('Raw Phase'); grid on
    % subplot(312); plot([t_axis,60], unwrap_dcRemove_ag); xlabel('Time (s)'); ylabel('Amplitude'); title('Unwrapped Phase'); grid on
    % subplot(313); plot(t_axis, diff_ag); xlabel('Time (s)'); ylabel('Amplitude'); title('Phase Difference'); grid on

    figure;
      subplot(141); plot(t_axis, breath_wave); title('Respiratory Signal'); xlabel('Time (s)'); ylabel('Amplitude'); grid on
       subplot(143); plot(xAxis, breath_fft); title('Respiratory Spectrum'); xlabel('Frequency (breaths/min)'); ylabel('Amplitude'); grid on
       subplot(142); plot(t_axis, heart_wave); title('Heartbeat Signal'); xlabel('Time (s)'); ylabel('Amplitude'); grid on
       subplot(144); plot(xAxis, heart_fft); title('Heartbeat Spectrum'); xlabel('Frequency (beats/min)'); ylabel('Amplitude'); grid on

end
