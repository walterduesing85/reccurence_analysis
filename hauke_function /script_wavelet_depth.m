% Script for wavelet power spectrum of Chew Bahir and other records.
%
% 1. Computes wavelet power spectrum of Chew Bahir and other records.
% 2. Stores the result in the 200th column of |data| .
% 3. Creates string array |datacorrstring|.
%
% 4 Sep 2019 - Trauth

% Interpolation to an equaly-spaced depth scale.
t = wavet(1):wavet(2):wavet(3);
x = interp1(newdata_a9(:,1),newdata_a9(:,2),t,'pchip');

% Replace data by interpolated data
clear newdata_a9
newdata_a9(:,1) = t;
newdata_a9(:,2) = x;

% Normalize data for wavelet transform. Most variables are not Gaussian but
% we are not interested in the exact statistics. We could also use medians
% and quartiles instead but this does not change the final result.
x = (x - mean(x))./std(x);

% Compute wavelet power spectrum.
fb = cwtfilterbank('Wavelet',wavename,...
    'SamplingFrequency',1/wavet(2),...
    'SignalLength',length(x),...
    'WaveletParameters',wavepar);
[wt,fr,coi] = cwt(x,'FilterBank',fb);

% Create string for graphics.
datawavestring = strcat("Wavelet power spectrum of ",titlestr_a9);


