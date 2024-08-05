% Author: Timothée Proix
% License: GPL-3.0-only

function [amp, phase, fX] = acausal4(X, muAll, option)
    samplefreq = 24;
    iCnt = 1;
    for imuAll = 1:length(muAll)
        mu = muAll{imuAll};
        for imu=1:length(mu)

            bdw = 1/3*1/mu(imu);
            cutofffreq = [1/mu(imu)-bdw, 1/mu(imu)+bdw];
            lengthFilter = min(2*mu(imu)*24, size(X, 1)/3-1);
            bPass = fir1(lengthFilter, cutofffreq/(samplefreq/2), 'bandpass');
            filtSig = filtfilt(bPass, 1, X(:, imuAll));

            phase(:, iCnt) = angle(hilbert(filtSig));
            amp(:, iCnt) = abs(hilbert(filtSig));
            fX(:, iCnt) = filtSig - mean(filtSig);
            iCnt = iCnt + 1;
        end
    end
end