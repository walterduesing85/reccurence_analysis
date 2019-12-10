% Script to read the Chew Bahir TC/TIC/TOC data from Frank
%
% 1. Reads data from Frank's file.
% 2. Interpolates data to evenly spaced time axis to create a new double
%    array |datatctictocage|.
%
%     1 = age (kyrs BP)
%     2 = TC (%)
%     3 = TIC (%)
%     4 = TOC (%)
%
% 3. Creates string array |datatctictocstring|.
%
% 26 Aug 2019 - Trauth

% Reading data from Frank's file
datatctictoc = load('data_frank_corg_2018_09_28.txt');

% Sorting.
datatctictoc = sortrows(datatctictoc,1);

% Removing duplicate data points.
datatctictoc2 = datatctictoc;
for i = 1:size(datatctictoc2,1)-1
    if datatctictoc2(i,1) == datatctictoc2(i+1,1)
        datatctictoc2(i,1) = NaN;
    end
end
datatctictoc2(isnan(datatctictoc2(:,1))==1,:) = [];
datatctictoc = datatctictoc2;
clear datatctictoc2

% Replacing NaNs by zeros.
datatctictoc(isnan(datatctictoc)==1) = 0;

% Create string for graphics.
datatctictocstring = ["CHB TC";
                      "CHB TIC";
                      "CHB TOC"];
                   
clear i 
