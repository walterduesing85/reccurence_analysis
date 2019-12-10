% Script to pretreat data
%
% 1. Highpass filtering
%
% 12 Sep 2019 - Trauth

% Revert age axis for RP/RQA.
data(:,1) = -data(:,1);

% Select data to be analyzed. The array data then contains (1) the age and
% (2) the variable to be analyzed.
t = data(:,1); t = t';
x = data(:,varselectnum)./data(:,varselectdem); x = x';

% Highpass filtering to remove long-term trend, with optional display of
% magnitude response.
samplingint = abs(mean(diff(data(:,1))));
samplingfreq = 1/samplingint;
nyquistfreq = abs(0.5 * mean(diff(data(:,1)))^(-1));
if filteroption == 1
    [b,a] = butter(filterorder,filtercutoff/nyquistfreq,'high');
    [hfilt,wfilt] = freqz(b,a,1024);
    f = samplingfreq*wfilt/(2*pi);
    x = filtfilt(b,a,x);
end

% Length lenTS of time series to be used with moving windows.
lenTS = length(x);