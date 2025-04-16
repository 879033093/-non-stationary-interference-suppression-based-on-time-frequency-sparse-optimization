function [pitchSpectrogram, azimuSpectrogram, Rxv] = IWR1642ODS_DOA(rangeProfile, doaFlag, useChirpNum, searchAngleRange, wopt_flag, wopt)
% Function: Perform azimuth and elevation DOA estimation for one frame of data from the IWR1642ODS radar and generate angle spectrograms

% Inputs:   
      % rangeProfile: One frame of range profile data after background subtraction
      % doaFlag: DOA mode flag: doaFlag = 0 (FFT), doaFlag = 1 (CBF), doaFlag = 2 (MVDR)
     % useChirpNum: Number of chirps used for DOA estimation
     % searchAngleRange: Search angle range for CBF or MVDR in degrees

% Outputs:
       % pitchSpectrogram: Elevation angle spectrogram
       % azimuSpectrogram: Azimuth angle spectrogram
       % Rxv: Inverse of the covariance matrix

% DOA Modes
FFT_DOA_MODE   = 0;
CBF_DOA_MODE   = 1;
MVDR_DOA_MODE  = 2;

% Only CBF and MVDR use 1D methods to compute angle spectrograms
rangeBinNum = size(rangeProfile, 2);
lambda = 5e-3;   % Wavelength
d = lambda / 2;    % Antenna spacing

if (doaFlag == CBF_DOA_MODE) || (doaFlag == MVDR_DOA_MODE)
    SEARCH_ANGLE_RANGE = searchAngleRange / 180 * pi;    % Convert search range to radians
    SEARCH_ANGLE_SPACE = 1 / 180 * pi;                     % Angular resolution in radians  
    fai = 2 * pi * sin(-SEARCH_ANGLE_RANGE:SEARCH_ANGLE_SPACE:SEARCH_ANGLE_RANGE) * d / lambda;
    
    azimuSpectrogram = zeros(length(fai), rangeBinNum);    % Azimuth angle spectrogram
    pitchSpectrogram = zeros(length(fai), rangeBinNum);  % Elevation angle spectrogram (not computed here but initialized)

    for r = 1:rangeBinNum
         % Extract 8-antenna signals for range bin r over useChirpNum chirps
        xt = squeeze(rangeProfile(1:useChirpNum, r, :))';  % xt: 8 x useChirpNum
        
        % Covariance matrix calculation
        Rx = xt * xt';   % 8x8 covariance matrix
        for an = 1:length(fai)
            % Steering vector
            aTheta = [1, exp(-1j*1*fai(an)), exp(-1j*2*fai(an)), exp(-1j*3*fai(an)), ...
                         exp(-1j*4*fai(an)), exp(-1j*5*fai(an)), exp(-1j*6*fai(an)), exp(-1j*7*fai(an))];

            if (doaFlag == CBF_DOA_MODE)  % Classical Beamforming
                azimuSpectrogram(an, r) = abs(aTheta * Rx * aTheta');
                Rxv = pinv(Rx);  
            else  % MVDR
                Rxv = pinv(Rx);  
                azimuSpectrogram(an, r) = 1 / abs(aTheta * Rxv * aTheta');
            end
        end
    end
end

% Normalize azimuth spectrogram
azimuSpectrogram = azimuSpectrogram / max(azimuSpectrogram(:));

% Optional: dB plot
% azimuSpectrogram = 10 * log10(azimuSpectrogram);
% figure;
% plot(20 * fai, 10 * log10(azimuSpectrogram)); xlabel('Angle (°)'); ylabel('dB'); title('MVDR DOA Estimation');

end
