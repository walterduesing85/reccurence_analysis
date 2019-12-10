% Script to read the Bosumtwi data (Miller et al. 2016)
%
% 1. Reads data from file.
% 6. Interpolates data to evenly spaced time axis to create a new double
%    array |dataicdpbosumtwiage|:
%
%     1 = age (kyrs BP)
%     2 = Poaceae (%)
%
% 7. Creates string array |dataicdpbosumtwistring| with the type of data.
%
% 30 Aug 2019 - Trauth

% Read data from file.
bosumtwi_poaceae = load('data_bosumtwi_poaceae_miller.txt');

% Create string array for graphics.
dataicdpbosumtwistring = "ICDP Bosumtwi Poaceae";
                    
% Remove rows without data
bosumtwi_poaceae(isnan(bosumtwi_poaceae(:,1))==1,:) = [];

% Interpolating data to evenly spaced time axis.
dataicdpbosumtwiage(:,1) = agemodelmin : agemodelres : agemodelmax;
dataicdpbosumtwiage(:,2) = interp1(bosumtwi_poaceae(:,1),...
        bosumtwi_poaceae(:,2),dataicdpbosumtwiage(:,1),inttype);

clear bosumtwi*