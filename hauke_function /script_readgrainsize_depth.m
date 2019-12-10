% Script to read the Chew Bahir grainsize data from Frank
%
% 1. Reads data from Frank's file.
% 2. Interpolates data to evenly spaced time axis to create a new double
%    array |datagrainsizeage|.
%
%     1 = age (kyrs BP)
%     2 = EM1
%     3 = EM2
%     4 = EM3
%     5 = EM4
%     6 = EM5
%     7 = Mean
%
% 3. Creates string array |datagrainsizestring|.
%
% 26 Aug 2019 - Trauth

% Reading data from Frank's file
datagrainsize = load('data_frank_grainsize_endmembers_2018_09_28.txt');
datagrainsizemean = load('data_frank_grainsize_mean_28_09_2018.txt');

% Sorting.
datagrainsize = sortrows(datagrainsize,1);
datagrainsizemean = sortrows(datagrainsizemean,1);

% Removing duplicate data points.
datagrainsize2 = datagrainsize;
for i = 1:size(datagrainsize2,1)-1
    if datagrainsize2(i,1) == datagrainsize2(i+1,1)
        datagrainsize2(i,1) = NaN;
    end
end
datagrainsize2(isnan(datagrainsize2(:,1))==1,:) = [];
datagrainsize = datagrainsize2;
clear datagrainsize2

datagrainsizemean2 = datagrainsizemean;
for i = 1:size(datagrainsizemean2,1)-1
    if datagrainsizemean2(i,1) == datagrainsizemean2(i+1,1)
        datagrainsizemean2(i,1) = NaN;
    end
end
datagrainsizemean2(isnan(datagrainsizemean2(:,1))==1,:) = [];
datagrainsizemean = datagrainsizemean2;
clear datagrainsizemean2

% Replacing NaNs by zeros.
datagrainsize(isnan(datagrainsize)==1) = 0;
datagrainsizemean(isnan(datagrainsizemean)==1) = 0;

% Create string for graphics.
datagrainsizestring = ["CHB Grainsize EM1";
                       "CHB Grainsize EM2";
                       "CHB Grainsize EM3";
                       "CHB Grainsize EM4";
                       "CHB Grainsize EM5"];
datagrainsizemeanstring =  "CHB Grainsize Mean";

clear i 
