function [harmonics, amplitudes] = generateHarmonicsNoLim(specEnv, fERB, fNAT, F0, nHar)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function computes the harmonics (frequencies) and corresponding
% amplitudes for a harmonic sound with given fundamental frequency and
% spectral envelope
%
% Input:
% - specEnv -> Spectral envelope
% - fERB    -> ERB-frequency vector
% - fNAT    -> Natural frequency vector
% - F0      -> Fundamental frequency
% - nHar    -> Number of harmonics
%
% Output:
% - harmonics   -> Vector with location of harmonics / Hz
% - amplitudes  -> Vector with amplitude of harmonics / a.u. (normalized to 1)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Round F0
F0 = round(F0);

% Interpolate spectral envelope from ERB to natural frequencies
interpSpecEnv = interp1(fERB, specEnv, fNAT, 'spline');

% Define harmonics for specific F0
harmonics = F0.*linspace(1, nHar, nHar);

% Exclude harmonics above 12 kHz
harmonics(harmonics>max(fNAT)) = nan;
notNANIdx = ~isnan(harmonics);
harmonics = harmonics(notNANIdx);

% % Make all harmonics vectors the same length
% if length(harmonics) > nHar
%     harmonics = harmonics(1:nHar);
% elseif length(harmonics) < nHar
%     lH = length(harmonics);
%     harmonics(lH+1:nHar) = 0;
% end

% Make amplitude spike train for harmonics
amplitudes = zeros(1,length(fNAT));
idx = ismember(fNAT, harmonics);
amplitudes(idx) = 1;

% Compute correct amplitudes for harmonics
amplitudes = interpSpecEnv.*amplitudes;     % Multiply logical amps with SE
amplitudes = amplitudes(idx);  % Extract only harmonics
amplitudes = 10.^(amplitudes./20);          % Convert to linear amplitude
% amplitudes = exp(amplitudes./10);         % Convert to linear amplitude
amplitudes = amplitudes./rms(amplitudes);   % RMS normalization

% Make all amplitudes vectors the same length
% if length(amplitudes) < 32
%     lA = length(amplitudes);
%     amplitudes(lA+1:nHar) = 0;
% end

