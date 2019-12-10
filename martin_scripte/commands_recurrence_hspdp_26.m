% By Martin H. Trauth and K. Hauke Kraemer 2019
%
% Version history:
%
% 12 Sep 2019 - Trauth
% 17 Sep 2019 - Kraemer
% 24 Oct 2019 - Trauth

% Clear workspace, clear command window and close figure windows
clear, clc, %close all

% Select variable to be analyzed from content_of_variable_data.txt.
varselectnum = 36;
varselectdem = 20;

% Definition of parameter values.
agemodeloption = 1;     % Choose age model: 
                        % 1=mubawa,2=oxcal,3-5 trauth 550-70, 6-8 merge
                        % oxcal+trauth 500/550/570
agemodelmin = 1;        % Minimum age (in kyrs).
agemodelmax = 570;      % Maximum age (in kyrs).
agemodelres = 0.1;      % Resolution of time axis (in kyrs).
inttype = 'pchip';      % Interpolation method.

% YLim for time series plot, potassium is 0 to 200000.
ymin = -Inf;
ymax = +Inf;

% Define time interval and spacing to be analyzed.
% Dry-wet 546, 515, 452, 410, 293, 258, 157, 133 kyr BP
% Wet-dry 299, 245 kyr BP
% Wet-dry gradual 314, 461, 228 kyr BP
% Dry-wet gradual 234 kyr BP
tmin = -agemodelmax;
tmax = -agemodelmin;
intv = agemodelres;

% Define axis settings for display of the results
xmin =  tmin;
xmax =  tmax;
dtick = 2000;

% Mask gaps, no=0, yes=1.
gapsoption = 1;            

% Saving figure, no=0, yes=1.
printres = 1;

% Optional filter to remove trend. Please use filter=1 if you would like to
% filter the data. Choose the cutoff frequency for the filter. The sampling
% interval of the data is 10 yrs corresponding to a sampling frequency of
% 0.1 kyrs^(-1). Therefore the Nyquist frequency is half of the sampling
% frequency, i.e. 0.05 kyrs^(-1). As an example, we can use a cutoff
% frequency of 0.01 kyrs^(-1) or a period of 100 kyrs, i.e. everything
% faster than 100 kyrs is preserved by the filter. Use filteroption = 1 if
% you wish to use the filter.
filteroption = 0;
filterorder  = 5;
filtercutoff = 0.1;

% Define embedding dimension m and time delay tau.
m = 2;
tau = 1;
timespan_diff = tau*(m-1);

% Define Theiler window and minimal line length.
theiler = tau;
l_min = 3;

% Define window size (w) and step size (ws).
w  = 500;
ws = 50;

% Set correction for border effects on (line_correct = 1) or off
% (line_correct = 0)
line_correct = 1;

% Select one RQA measures, which will be displayed in a separate figure
% together with RR and DET in the first window
%     Y(1) = RR     (recurrence rate)
%     Y(2) = DET    (determinism)
%     Y(3) = <L>    (mean diagonal line length)
%     Y(4) = Lmax   (maximal diagonal line length)
%     Y(5) = ENTR   (entropy of the diagonal line lengths)
%     Y(6) = LAM    (laminarity)
%     Y(7) = TT     (trapping time)
%     Y(8) = Vmax   (maximal vertical line length)
%     Y(9) = RTmax  (maximal white vertical line length)
%     Y(10) = T2    (recurrence time of 2nd type)
%     Y(11) = RTE   (recurrence time entropy, i.e., RPDE)
%     Y(12) = Clust (clustering coefficient)
%     Y(13) = Trans (transitivity)
RQA_Select = 12;

% Calculate recurrence plot using the recurrence threshold e and choose
% threshold-calculation parameter. Set parameter to define the threshold-
% calculation for the recurrence plot estimation in the next step. There
% are three options to choose from:
%   'fix'  The RP is computed under a fixed threshold e specified by input
%          parameter e.
%   'var'  The RP is computed under a fixed threshold e, which corresponds
%          to the lower e-quantile (specified by input parameter e of the 
%          distance distribution of all points in phase space.
%   'fan'  The RP is computed under a variable threshold e using a fixed
%          amount of nearest neighbours in phase space to compute the e-
%          value for each point of the phase space trajectory.
threshold_calculation = 'fan';

% Choose norm from ['euc','max']
norm = 'euc';

% If you want to receive a diagonal RP, which corrects for tangential
% motion (consisting of just diagonal lines, after Kraemer & Marwan 2019)
% tick "diagonal = 1", if you want to work with a normal RP tick
% "diagonal = 0".
diagonal = 0;

% Importing Chew Bahir age model and data.
script_readagemodel
script_readxrf_age
script_readxrd_age
script_readcolref_age
script_readphysprop_age
script_readgrainsize_age
script_readtctictoc_age

% Importing other data.
script_readorbitalforcing_age
script_readodp967all_age
script_readodp722sst_age
script_readmalawi_age
script_readco2_age
script_readbosumtwi_age

% Merging all data into a single variable |data|.
data(:,1)     = datacolrefpcaage(:,1);    % Age (kyrs BP)
data(:,2:6)   = datacolrefpcaage(:,2:6);  % PC1-5 of Col Ref (LacCore)
data(:,7:20)  = dataxrfage(:,2:15);       % XRF + Ones (Duluth)
data(:,21:24) = dataodp967age(:,2:5);     % ODP 967 data (Grant 17)
data(:,25:27) = dataicdpmalawiage(:,2:4); % ICDP Malawi data (Johnson 16)
data(:,28:30) = dataorbitalage(:,2:4);    % Orbital parameters (Laskar 04)
data(:,31:35) = dataxrdage(:,2:6);        % XRD data (Foerster)
data(:,36)    = dataphyspropage(:,2);     % Mag sus (LacCore)
data(:,37:42) = datagrainsizeage(:,2:7);  % Grainsize (Schaebitz)
data(:,43:45) = datatctictocage(:,2:4);   % TC, TIC, TOC (Schaebitz)
data(:,46)    = dataepicaco2age(:,2);     % EPICA CO2 (Bereiter 15)
data(:,47)    = dataodp722sstage(:,2);    % ODP 722 SST (Herbert 10)
data(:,48)    = dataicdpbosumtwiage(:,2); % ICDP Bosumtwi data (Miller 16)

% Merging all the strings into a single variable |datastr|.
datastr(1,:)     = "Age (kyr)";
datastr(2:6,:)   = datacolrefpcastring;
datastr(7:19,:)  = dataxrfstring;
datastr(20,:)    = "Ones";
datastr(21:24,:) = dataodp967string;
datastr(25:27,:) = dataicdpmalawistring;
datastr(28:30,:) = dataorbitalstring;
datastr(31:35,:) = dataxrdstring;
datastr(36,:)    = dataphyspropstring;
datastr(37:42,:) = datagrainsizestring;
datastr(43:45,:) = datatctictocstring;
datastr(46,:)    = dataepicaco2string;
datastr(47,:)    = dataopd722sststring;
datastr(48,:)    = dataicdpbosumtwistring;

% Data pretreatment, compute RQ, perform RQA and display results.
script_pretreatment_age
script_recurrenceplot_age
script_recurrencequantificationanalysis_age
script_maskgaps_age
script_displayresults_rqa_age


