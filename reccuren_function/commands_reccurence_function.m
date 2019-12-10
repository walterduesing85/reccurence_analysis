%This script uses the recurrence function of Kraemer et al. 2019 

clc 

close all 

clear 


dataxrf = load('data_raw_xrf_11_11_2019.txt');
% % % 

% % % 
  ages = load('data_ages_mubawa_11_11_2019.txt');
Vn1 = 4;
%   data1 = [dataxrf(:,1) dataxrf(:,Vn1)];
% % 
% data1 = load('data_magsus.txt');
% % 
  data1 = [dataxrf(:,1) log(dataxrf(:,Vn1)./dataxrf(:,14))];

%apply age model
[data_1,inv] = agemodel_2(ages,data1,650);
%%

% Input-parameters:
%
% 't'                       time vector
%
% 'y'                       corresponding values
%
% 'tau'                     time-delay (for embedding procedure)
%
% 'm'                       embedding dimension
%
% 'T'                       Theiler-window
%
% 'L'                       minimal diagonal & vertical line length
%                           
% 'w','ws'                  windowsize & stepsize for the RQA-runningwindow
%                           approach
% 'norm'                    norm for distance calculation in phasespace to
%                           'euc' (euclidic) or 'max' (maximum). Default is 
%                           euclidic.
% 'RQA-quantity'            defines which RQA-measurement is shown in the
%                           output figures
% 'window-shape'            If set to 1, the RQA is based on a 45°-rotated 
%                           window. Default is 0.
% 'line_correct'            Optional input 'line_correct' can be used in 
%                           border to correct for border effects while 
%                           counting the diagonal line length histograms 
%                           (Kraemer & Marwan, PLA, 2019).
%                           If 'line_correct = 1', the "kelo"-correction is
%                           applied,
%                           if 'line_correct = 0', the conventional method 
%                           is applied (Default)
% 'diagonal_RP'             if set to '1', a new, skeltonized recurrence 
%                           plot, which just contains diagonal lines is 
%                           computed (Default is '0'). Accounts
%                           for the effects of tangential motion and has
%                           been proposed in Kraemer & Marwan, PLA, 2019
% 'ShowAll'                 should be set to '1', if you want to receive two 
%                           additional figures containing all RQA-variance 
%                           measurements. If you don't want this, just leave
%                           it or set ShowAll = 0.
%
% 'threshold' & 'threshold-calc':
%
% Input Paramter 'threshold-calc' specifies how the threshold epsilon will
% be calculated. There are three options. Set "threshold-calc" to
%   - 'fix' The RP is computed under a fixed threshold epsilon specified by
%           input parameter 'threshold'.
%   - 'var' The RP is computed under a fixed threshold epsilon, which
%           corresponds to the lower 'threshold'-quantile (specified by input
%           parameter 'threshold') of the distance distribution of all points
%           in phasespace.
%   - 'fan' The RP is computed under a variable threshold epsilon using a
%           fixed amount of nearest neighbours in phasespace to compute the
%           epsilon-value for each point of the phasespace trajectory
%           individually (input variable 'threshold' in this case
%           corresponds to the percentage of all point to be considered as
%           "nearest neighbours"
% Default is 'fix'.    
%
% 'running-window'          (optional) specifies the type of the RQA running 
%                           window approach. If set to 0 (default), the 
%                           running window is sliding over the RP, which has
%                           been reconstructed from the whole time series. 
%                           If set to 1, the running window slides over the 
%                           time series, and each part of the time series 
%                           (of size w) will be used to construct an RP
%                           and do the RQA. 
%
%



t1        = data_1(:,1); %t
y1        = data_1(:,2);
m        = 2;
tau      = 1;
threshold = 'fix';
T        = tau;
L        = 2;
w        = 500;
ws       = 50;
RQA      = 1;
norm      = 'euc';
% threshold_calc  = 'fan';
% % diagonal_RP             =
% % line_correct            =
% window_shape            = 
% % ShowAll                 = 

varargin{1} = t1;
varargin{2} = y1;
varargin{4} = tau;
varargin{3} = m;
varargin{5} = threshold;
varargin{6} = tau;
varargin{7} = L;
varargin{8} = w;
varargin{9} = ws;
varargin{10} = RQA;
[Z1,Z2,DM] = rqaplot(varargin);


