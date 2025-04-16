function Cancel_PulseCompre_Data = MTI_PulseCompression(framesData, cancelFlag, rangeFFTNum, pulseCancelSpace)
% Function: Perform clutter cancellation and pulse compression on one frame of radar echo data

% Inputs:   
% framesData: Raw radar echo data (before cancellation and pulse compression)
% cancelFlag: Cancellation mode flag:
%             cancelFlag = 0 → Mean cancellation
%             cancelFlag = 1 → Two-pulse cancellation
%             cancelFlag = 2 → Three-pulse cancellation
% rangeFFTNum: Number of FFT points for range pulse compression
% pulseCancelSpace: Chirp interval for two-pulse cancellation (used only when cancelFlag = 1)

% Output:
% Cancel_PulseCompre_Data: Complex data after pulse compression

%% Perform clutter cancellation and pulse compression
meanCancel_mode = 0;
twoPulseCancel_mode = 1;
threePulseCancel_mode = 2;

numChirp = size(framesData, 1);  % Number of chirps

% --------- Mean Cancellation
if (cancelFlag == meanCancel_mode)
    % Subtract the mean across chirps (DC removal)
    meanCancelData = framesData - repmat(mean(framesData, 1), [numChirp 1]);
    
    % Perform inverse FFT along range dimension (dimension 2)
    Cancel_PulseCompre_Data = ifft(meanCancelData, rangeFFTNum, 2);
end

% --------- Two-Pulse Cancellation
if (cancelFlag == twoPulseCancel_mode)
    % Subtract delayed chirps to remove static/clutter components
    twoPulseCancelData = framesData(1:end-pulseCancelSpace, :, :) - framesData(1+pulseCancelSpace:end, :, :);
    
    % Pulse compression using IFFT
    Cancel_PulseCompre_Data = ifft(twoPulseCancelData, rangeFFTNum, 2);
end

% --------- Three-Pulse Cancellation
if (cancelFlag == threePulseCancel_mode)
    % Apply second-order difference: x(n-1) - 2x(n) + x(n+1)
    threePulseCancelData = framesData(1:end-2, :, :) - 2 * framesData(2:end-1, :, :) + framesData(3:end, :, :);
    
    % Pulse compression using IFFT
    Cancel_PulseCompre_Data = ifft(threePulseCancelData, rangeFFTNum, 2);
end
