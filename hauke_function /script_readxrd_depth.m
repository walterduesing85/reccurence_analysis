% Script to read the Chew Bahir XRD data from Verena
%
% 1. Reads data from Verena's file.
% 2. Interpolates data to evenly spaced time axis to create a new double
%    array |dataxrdage|.
%
%     1 = age (kyrs BP)
%     2 = smectite abundance group
%     3 = illite abundance group
%     4 = analcime group
%     5 = analcime intensity
%     6 = smectite inverted group
%
% 3. Creates string array |dataxrdagestring|.
%
% 26 Aug 2019 - Trauth

% Reading data from Verena's file
dataxrd = load('data_verena_xrd_2019_03_28.txt');

% Sorting.
dataxrd = sortrows(dataxrd,1);

% Removing duplicate data points.
dataxrd2 = dataxrd;
for i = 1:size(dataxrd2,1)-1
    if dataxrd2(i,1) == dataxrd2(i+1,1)
        dataxrd2(i,1) = NaN;
    end
end
dataxrd2(isnan(dataxrd2(:,1))==1,:) = [];
dataxrd = dataxrd2;
clear dataxrd2

% Replacing NaNs by zeros.
dataxrd(isnan(dataxrd)==1) = 0;

% Create string for graphics.
dataxrdstring = ["CHB Smectite Abundance Group";
                 "CHB Illite Abundance Group";
                 "CHB Analcime Group";
                 "CHB Analcime Intensity";                    
                 "CHB Smectite Inverted Group"];

clear i
