%% IMPORTING, PROCESSING, ANALYZING AND DISPLAYING CHEW BAHIR DATA VS DEPTH
% This code imports all Chew Bahir data and displays them together with a
% lithologic column vs. composite depth.
%
% 5 Sep 2019 - Trauth, the first version of the script

% Clear workspace, clear command window and close figure windows
clear, clc, close all

% Preferences for time series. Please select 8 variables from the the file
% content_of_variable_data.txt to compute and display variables. Select a
% 9th variable to compute and display a wavelet power spectrum
varselect = [3 9 10 12 18 26 32 33 3];

% Preferences for wavelet power spectrum.
wavevar = varselect(9); % Index of variables for wavelet power spectrum of
                        % the variables.   
wavename = 'morse';     % Wavelet name.
wavepar = [3,30];       % Wavelet parameters.
wavet   = [1 0.05 292]; % Depth axis for interpolation, [tmin tinv tmax].
wavexlim = [0 300];     % YLim for periods in wavelet power spectrum.
waveylim = [20 150];     % YLim for periods in wavelet power spectrum.
wavedisp = 2;           % Type of display of wavelet power spectrum,
                        % 1=pcolor, 2=contourf.
waveaxis = [0 1];       % Pseudocolor axis scaling for wavelet power.
                        % spectrum. Default is [0 Inf]. Since the data are
                        % normalized prior to wavelet transformation the
                        % values of [0 1] work well for most of them.

% Importing Chew Bahir data.
script_readstratcol
script_readxrf_depth
script_readxrd_depth
script_readcolref_depth
script_readphysprop_depth
script_readgrainsize_depth
script_readtctictoc_depth

% Merging all data into a single variable |data|.
data{1}     = datacolrefpca;     % PC1-5 of Col Ref (LacCore)
data{2}     = dataxrf;           % XRF + Ones (Duluth)
data{3}     = dataxrd;           % XRD data (Foerster)
data{4}     = dataphysprop;      % Mag sus (LacCore)
data{5}     = datagrainsize;     % Grainsize endmembers (Schaebitz)
data{6}     = datagrainsizemean; % Grainsize mean (Schaebitz)
data{7}     = datatctictoc;      % TC, TIC, TOC (Schaebitz)

% Merging all the strings into a single variable |datastr|.
datastr{1}  = datacolrefstring;
datastr{2}  = dataxrfstring;
datastr{3}  = dataxrdstring;
datastr{4}  = dataphyspropstring;
datastr{5}  = datagrainsizestring;
datastr{6}  = datagrainsizemeanstring;
datastr{7}  = datatctictocstring;

% Managing data for display.
script_managedata

% Wavelet Power Spectrum.
script_wavelet_depth

% Display data with linked x-axes for zooming and browsing.
script_displayresults_depthseries_depth
script_displayresults_wavelet_depth
script_displayresults_wavelet_depth_extra

% Clearing variables that are no longer required.
%clearvars -except data datastr varselect wt

whos

