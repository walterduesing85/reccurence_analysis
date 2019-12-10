function [Z1,Z2,DM] = rqaplot(varargin)
% RQAPLOT   Plots a measure of Recurrence Quantification Analysis
%
%   version 1.8
%   new:    - implementation of border line corrections
%
%
% Input: 
%
% Minimum input-arguments : 10
% Maximum input-arguments : 17
%
% [Z1,Z2,DM] = rqaplot(t,y,m,tau,threshold,T,L,w,ws,RQA-quantity,norm,...
%                       threshold-calc,window-shape,diagonal_RP,...
%                       line_correct,running-window,ShowAll)
%
% This function plots the measures of Recurrence Quantification Analysis
% (RQA), specified in input argument 'RQA-quantity', the corresponding time
% series as well as the corresponding RecurrencePlot (RP).
%
% Returns: Three figures and a phase-corrected vector containing all 18 
%          RQA-Measurements.
%          - Figure 1 contains the time-series with the choosen RQA-Quantity
%          - Figure 2 contains the RP with the choosen RQA-Quantity 
%          - Figure 3 contains the RP with the time-series.
%          On request the function returns 2 additional figures, containing
%          all 18 RQA-Measurements.
%
%          RQA-Measurements can also be stored in output-varibales Z1,Z2 & DM.
%          - Z1 contains all RQA-measurements of the moving window approach. 
%          - Z2 contains all RQA-measurements based on the whole RP 
%            (no moving window). 
%          - DM contains the unthresholded distance matrix.
%
% The RP will be calculated from the embedded time-series with embedding-
% parameters tau (time-delay) and m (embedding dimension) under a threshold 
% epsilon.
%
% The certain RQA-quantities are described as follows:
%
%      1 = RR       (recurrence rate)
%      2 = DET      (determinism)
%      3 = <L>      (mean diagonal line length)
%      4 = Lmax     (maximal diagonal line length)
%      5 = ENTR     (entropy of the diagonal line lengths)
%      6 = LAM      (laminarity)
%      7 = TT       (trapping time)
%      8 = Vmax     (maximal vertical line length)
%      9 = Ratio1   (ratio between T2 and RR)
%      10= T2       (recurrence time of 2nd type)
%      11= RTE      (recurrence time entropy, i.e., RPDE)
%      12= Clust    (clustering coefficient)
%      13= Trans    (transitivity)
%      14= Ratio2   (ratio between DET and RR)
%      15= Ratio3   (LAM/DET)
%      16= DIV      (1/Lmax)
%      17= RTmax    (maximal white vertical line length)
%      18= RT_DIV   (1/RTmax)
%      19= 1/RR_col (inverse of columnwise RR) (just valid for normal window shape and
%                                  RQA based on running window over global RP)
%      20= windowed 1/RR_col 
%
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
% Copyright (c) 2017
% Hauke Krämer, 
% Potsdam Institute for Climate Impact Research, Germany
% http://www.pik-potsdam.de
% Insitute of Environmental and Earth Science, University of Potsdam,
% Germany
% http://www.geo.uni-potsdam.de
% Based on work of Norbert Marwan, Potsdam Institute for Climate Impact Research, Germany
% (modified by Martin H. Trauth, Insitute of Environmental and Earth Science, University of Potsdam,
% Germany)
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or any later version.

%% Assign input
x = varargin{1}(:);
y = varargin{2}(:);
tau = varargin{4};
m = varargin{3};
epsilon = varargin{5};
T = varargin{6};
L = varargin{7};
w = varargin{8};
ws = varargin{9};
RQA = varargin{10};

methLib={'euc','max'}; % the possible norms
try
    norm = varargin{11};
    if ~isa(norm,'char') || ~ismember(norm,methLib)
       warning(['Specified norm should be one of the following possible values:',10,sprintf('''%s'' ',methLib{:})])
    end
catch
    norm = 'euc';
end

thresLib={'fix','var','fan'}; % the possible ways of threshold computation
try
    thres = varargin{12};
    if ~isa(thres,'char') || ~ismember(thres,thresLib)
       warning(['Specified way of calculating threshold should be one of the following possible values:',10,sprintf('''%s'' ',thresLib{:})])
    end
catch
    thres = 'fix';
end

try
    if varargin{13}==1
        windowshape = varargin{13};
    elseif varargin{13}==0
        windowshape = varargin{13};
    else
        windowshape = 0;
    end
catch
    windowshape = 0;
end

try
    if varargin{14}==1
        diagonal_rp = varargin{14};
    else
        diagonal_rp = 0;
    end
catch
    diagonal_rp = 0;
end

try
    if varargin{15}==1
        line_correct = varargin{15};
    else
        line_correct = 0;
    end
catch
    line_correct = 0;
end


try
    if varargin{16}==1
        runningwindow = varargin{16};
    else
        runningwindow = 0;
    end
catch
    runningwindow = 0;
end

try
    if varargin{17}==1
        ShowAll = varargin{17};
    else
        ShowAll = 0;
    end
catch
    ShowAll = 0;
end



timespan_diff = tau*(m-1);      % bind the timespan-difference covered by the embedding vectors

lenTS = length(y);              % length of the time-series

% To preserve the number of points in the moving window, when using the new
% window-shape, this windowlength of the new window, as well as the
% step-size has to be adjusted:

s = sqrt(w^2/2);
w_new = round(2*s);          % new window-size is twice the size of s
overlap = 1-(ws/w);
ws_new = round((1-overlap)*w_new);


%% check input
narginchk(10,17)
nargoutchk(0,3)
    
% check input variables x,y whether they are equal in size
if size(x,2)~=size(y,2) || size(x,1)~=size(y,1)
    error('Input vectors must have the same length.')
end

% check input variables x, y whether they are column- or line-vectors
if size(x,1)<size(x,2)
    x=x';
end

if size(y,1)<size(y,2)
    y=y';
end

% check input variabel w, whether it is an even integer
if rem(w,1)~=0 || rem(w,2)~=0
    error('Windowsize w must be an even integer.')
end

% check input variabel ws, whether it is an integer
if rem(ws,1)~=0
    error('Window-Step ws must be an integer.')
end
    
% check input variabel w_new, whether it is an uneven integer
if rem(w_new,1)~=0 || rem(w_new,2)==0
    w_new=w_new-1;
    %error('Windowsize of new window must be an uneven integer.')
end

% check input variabel ws_new, whether it is an integer
if rem(ws_new,1)~=0
    error('Window-Step of new window must be an integer.')
end

% check input variabel RQA
if (RQA == 19 && runningwindow == 1) || (RQA == 20 && runningwindow == 1)
    error('Columnwise RR just available for normal windowshape and running window over global RP')
elseif (RQA == 19 && windowshape == 1) || (RQA == 20 && windowshape == 1)
    error('Columnwise RR just available for normal windowshape and running window over global RP')
end

% check input variabel windowshape, if it is 0 or 1
if windowshape ~= 0 && rem(windowshape,1)~=1
    error('Window shape input parameter must be set to 0 (Default setting, corresponds to normal window shape) or 1 (45° rotated window)')
end


%% Calculate RQA in moving windows
if runningwindow == 0
    
    %% Embedding
    Y=embed(y,m,tau);

    %% Calculate RP
    if strcmp(thres,'fix')
        [R,DM,~] = rp(Y,epsilon,'fix',norm);
    elseif strcmp(thres,'var')
        [R,DM,~] = rp(Y,epsilon,'var',norm);
    elseif strcmp(thres,'fan')
        [R,DM,~] = rp(Y,epsilon,'fan',norm);
    end
    
    if diagonal_rp
        R = rp_diagonal(R);
    end
    
    RR = NaN(size(x,1),size(x,1));
    RR(1+round(timespan_diff/2):size(x,1)-floor(timespan_diff/2),1+round(timespan_diff/2):size(x,1)-floor(timespan_diff/2)) = R;

    % normal window shape
    if windowshape==0 || windowshape==2
        % perform RQA-measures for each window-moving time-step
        cnt = 1;                                            % counter. -> index for storing RQA-measures for each window-movement
        r_win = zeros(20,ceil((length(R)-w)/ws));

        for i = 1:ws:(length(R)-w)
            
            r_win(1:13,cnt) = rqa(R(i:(i+w),i:(i+w)),L,T,line_correct);
           
            % RT max
            r_win(17,cnt) = r_win(9,cnt);
            % ratio 1
            r_win(9,cnt) = r_win(10,cnt)/r_win(1,cnt);

            % ratio 2
            r_win(14,cnt) = r_win(2,cnt)/r_win(1,cnt);

            % ratio 3
            r_win(15,cnt) = r_win(6,cnt)/r_win(2,cnt);

            % DIV
            r_win(16,cnt) = 1/r_win(4,cnt);

            % RT-DIV
            r_win(18,cnt) = 1/r_win(17,cnt);
            cnt = cnt+1;
        end

        % Phase correction for windowed measures added by M.H. Trauth by shifting
        % the measures by half the window size and sampled at the resolution ws.

        r_win_w(1:size(r_win,1),1:size(R,1)) = NaN;
        r_win_w(:,1+w/2:ws:size(R,1)-w/2) = r_win;
        for i = 1:18
            r_win_w(i,1+w/2:size(R,1)-w/2) = ...
               fillmissing(r_win_w(i,1+w/2:size(R,1)-w/2),'linear');
        end
        
        % Calculate columnwise RR and store it in 19th column of r_win_w
        for j = 1:size(R,1)
            r_win_w(19,j) = 1/ (sum(R(:,j))/size(R,1));
%             r_win_w(19,j) = -log((sum(R(:,j))/size(R,1)));
        end
        r_win_w(19,:) = r_win_w(19,:)/max(r_win_w(19,:));
        % Perform running window over columnwise RR and store it in 20th 
        % column of r_win_w
        cnt = 1;
        for i = 1:ws:(length(R)-w)
           r_win(20,cnt) = nanmean(r_win_w(19,i:(i+w)));
           cnt = cnt+1;
        end
        % Phase correction for windowed measures by shifting the measures 
        % by half the window size and sampled at the resolution ws for 20th
        % column of r_win
        r_win_w(20,1+w/2:ws:size(R,1)-w/2) = r_win(20,:);
        r_win_w(20,1+w/2:size(R,1)-w/2) = ...
               fillmissing(r_win_w(20,1+w/2:size(R,1)-w/2),'linear');
        
        % Phase correction for embedding delay added by M.H. Trauth by shifting the
        % measures by half the embedding delay. 
        r_win_e(1:size(r_win,1),1:lenTS) = NaN;
        r_win_e(:,1+round(timespan_diff/2):lenTS-floor(timespan_diff/2)) = r_win_w;
        
        % Add gaps in the RQA measures where the data y has gaps.
        for i = 1 : size(r_win_e,1)
           r_win_e(i,isnan(y)==1) = NaN;
        end

        % Flag the NaNs by their relative location: The first NaN of a NaN block
        % is assigned the value of 1, the last NaN of the block gets the value 2,
        % and if the block consists of only one single NaN, then it is marked 3.
        index_vector_first = [];
        index_vector_last = [];
        nan_vector = isnan(y);
        flag = true;
        for i = 1 : size(r_win_e,2)
           if i == 1
             if nan_vector(i) == 0
                flag = false;
             end
           end
           if nan_vector(i) == 0 && flag == true && i == 1
             flag = false;
           elseif nan_vector(i) == 1 && flag == false 
             index_vector_first = [index_vector_first i];
             flag = true;     
           elseif nan_vector(i) == 0 && flag == true 
             index_vector_last = [index_vector_last i-1];
             flag = false;
           end
        end

        % NaN-padding before NaN block.
        for i = 1 : length(index_vector_first)
          j = index_vector_first(i);
          span = j - (w) - 2;
          gap = (w/2) + mod(span,ws);
          begin = j - gap;
  
          for k = 1 : size(r_win_e,1)
            r_win_e(k,begin:j) = NaN;
          end
        end

        % NaN-padding after NaN block.
        for i = 1 : length(index_vector_last)
          j = index_vector_last(i);
          span = j + (w);
          gap = (w/2) + mod(span,ws);
          ende = j + gap;

          for k = 1 : size(r_win_e,1)
             r_win_e(k,j:ende) = NaN;
          end
        end
 
        Z1=r_win_e;
    end
    
    
    %% RQA new-window

    if windowshape ~=0

        % perform RQA-measures for each window-moving time-step
        r_win2 = zeros(18,ceil((length(R)-w_new)/ws_new));   % allocate RQA-storing object.

        cnt2 = 1;                          % counter. -> index for storing RQA-measures
                                           % for each window-movement. 

        for startwert = 1:ws_new:length(R)-w_new+1 % Spalten-Startwert(e) für RP-Bildausschnitt 
                                                    % setzen.

            RP_Ausschnitt=zeros(w_new,w_new);  % leere Matrix mit Kantenlänge w machen. 
                                               % Hierein soll nun der RP-Ausschnitt
                                               % kommen. Danach wird diese für die RQA
                                               % verwendet.

            flag = false;
            maxcounter = round(w_new/2)-1;    % controls the maximum amount of chosen lines of old RP
            counter = 0;                    % controls the amount of chosen lines of old RP

            % Innerhalb dieses Ausschnitts jetzt Werte aus dem RP kopieren, die in 
            % der neuen Fensterform liegen.

            j2=1;                                   % Spaltenindex des RP_Auschnitt
            for j = startwert:startwert+w_new-1     % Spalten wählen
                % Zeilen wählen
                i2=round(w_new/2)-counter;          % Zeilenindex des RP_Ausschnitt
                for i = (startwert+round(w_new/2)-counter):(startwert+round(w_new/2)+counter)
                    RP_Ausschnitt(i2,j2)=R(i,j);    % Werte übernehmen
                    i2=i2+1;                        % Zeilenindex des RP_Ausschnitt anpassen
                end

                % wegen rautenform : auf maxcounter achten und dann mit dem
                % zeilenindex, welcher über counter gesteuert wird, wieder
                % heruntergehen.
                if counter == maxcounter
                    flag = true;
                end

                if flag == false
                    counter = counter + 1;
                end

                if flag == true
                    counter = counter - 1;
                end
                j2=j2+1;                % Spaltenindex des RP-Ausschnitts erhöhen
            end
            % RQA machen
            r_win2(1:18,cnt2)=rqa2(RP_Ausschnitt,L,T);


            % Spaltencounter der RQA's erhöhen
            cnt2=cnt2+1;

        end

        % Phase correction for windowed measures by shifting
        % the measures by half the window size and sampled at the resolution ws_new.

        r_win_w2(1:size(r_win2,1),1:size(R,1)) = NaN;
        r_win_w2(:,round(w_new/2):ws_new:size(R,1)-round(w_new/2)) = r_win2;

        for i = 1:18
            r_win_w2(i,round(w_new/2):size(R,1)-round(w_new/2)) = ...
               fillmissing(r_win_w2(i,round(w_new/2):size(R,1)-round(w_new/2)),'linear');
        end

        % Phase correction for embedding delay by shifting the
        % measures by half the embedding delay. 
        r_win_e2(1:size(r_win2,1),1:lenTS) = NaN;
        r_win_e2(:,1+round(timespan_diff/2):lenTS-floor(timespan_diff/2)) = r_win_w2;
        
        % Add gaps in the RQA measures where the data y has gaps.
        for i = 1 : size(r_win_e2,1)
           r_win_e2(i,isnan(y)==1) = NaN;
        end

        % Flag the NaNs by their relative location: The first NaN of a NaN block
        % is assigned the value of 1, the last NaN of the block gets the value 2,
        % and if the block consists of only one single NaN, then it is marked 3.
        index_vector_first = [];
        index_vector_last = [];
        nan_vector = isnan(y);
        flag = true;
        for i = 1 : size(r_win_e2,2)
           if i == 1
             if nan_vector(i) == 0
                flag = false;
             end
           end
           if nan_vector(i) == 0 && flag == true && i == 1
             flag = false;
           elseif nan_vector(i) == 1 && flag == false 
             index_vector_first = [index_vector_first i];
             flag = true;     
           elseif nan_vector(i) == 0 && flag == true 
             index_vector_last = [index_vector_last i-1];
             flag = false;
           end
        end

        % NaN-padding before NaN block.
        for i = 1 : length(index_vector_first)
          j = index_vector_first(i);
          span = j - (w_new) - 2;
          gap = (w_new/2) + mod(span,ws_new);
          begin = j - gap;
  
          for k = 1 : size(r_win_e2,1)
            r_win_e2(k,begin:j) = NaN;
          end
        end

        % NaN-padding after NaN block.
        for i = 1 : length(index_vector_last)
          j = index_vector_last(i);
          span = j + (w_new);
          gap = (w_new/2) + mod(span,ws_new);
          ende = j + gap;

          for k = 1 : size(r_win_e2,1)
             r_win_e2(k,j:ende) = NaN;
          end
        end        
        
        
        Z1=r_win_e2;


    end
    
% Running window over time series    
elseif runningwindow == 1
    
    if windowshape==0
        
        % preallocate vectors for storing RQA measurements 
        r_win = zeros(20,ceil((lenTS-timespan_diff-w)/ws));
        
        % preallocate cell-array for storing distance matrices
        DM = cell(1,ceil((lenTS-timespan_diff-w)/ws));
        
        % Control windowsteps
        startwert = 0;
        cnt=1;
        
        % Perform running window approach
        h = waitbar(0,'Performing RQA (running window over time series)');
        for i = 1:ws:lenTS-(w+timespan_diff) ,waitbar(i/(lenTS-(w+timespan_diff)))
            
           % Define Interval in time series
           yx = y(1+startwert:w+startwert+timespan_diff);

           % Embedding
           Y = embed(yx,m,tau);
           
           % RP
           if strcmp(thres,'fix')
                [R,P,~] = rp(Y,epsilon,'fix',norm);
           elseif strcmp(thres,'var')
                [R,P,~] = rp(Y,epsilon,'var',norm);
           elseif strcmp(thres,'fan')
                [R,P,~] = rp(Y,epsilon,'fan',norm);
           end  
           
           if diagonal_rp
                R = rp_diagonal(R);
           end
           
           DM{cnt} = P;
           % RQA
           r_win(1:13,cnt) = rqa(R,L,T,line_correct);
           % RT max
           r_win(17,cnt) = r_win(9,cnt);
           % ratio 1
           r_win(9,cnt) = r_win(10,cnt)/r_win(1,cnt);
           % ratio 2
           r_win(14,cnt) = r_win(2,cnt)/r_win(1,cnt);
           % ratio 3
           r_win(15,cnt) = r_win(6,cnt)/r_win(2,cnt);
           % DIV
           r_win(16,cnt) = 1/r_win(4,cnt);
           % RT-DIV
           r_win(18,cnt) = 1/r_win(17,cnt);
           
           % Calculate columnwise RR and store mean of it in 19th column of Z2
           helper = zeros(1,size(R,1));
           for j = 1:size(R,1)
               helper(1,j) = 1/ (nansum(R(:,j))/size(R,1));
            %             r_win_w(19,j) = -log((sum(R(:,j))/size(R,1)));
           end
           helper = helper/ max(helper);
           r_win(19,:) = nanmean(helper);
           r_win(20,:) = nanmean(helper);
           
           % Increase control parameters
           startwert = startwert + ws;
           cnt = cnt + 1;
           
        end
        close(h)
        
        % Phase correction for windowed measures added by M.H. Trauth by shifting
        % the measures by half the window size and sampled at the resolution ws.

        r_win_w(1:size(r_win,1),1:lenTS-timespan_diff) = NaN;
        r_win_w(:,1+w/2:ws:lenTS-timespan_diff-w/2) = r_win;
        for i = 1:20
            r_win_w(i,1+w/2:lenTS-timespan_diff-w/2) = ...
               fillmissing(r_win_w(i,1+w/2:lenTS-timespan_diff-w/2),'linear');
        end

        % Phase correction for embedding delay added by M.H. Trauth by shifting the
        % measures by half the embedding delay. 
        r_win_e(1:size(r_win,1),1:lenTS) = NaN;
        r_win_e(:,1+round(timespan_diff/2):lenTS-floor(timespan_diff/2)) = r_win_w;
        
        
        % Add gaps in the RQA measures where the data y has gaps.
        for i = 1 : size(r_win_e,1)
           r_win_e(i,isnan(y)==1) = NaN;
        end

        % Flag the NaNs by their relative location: The first NaN of a NaN block
        % is assigned the value of 1, the last NaN of the block gets the value 2,
        % and if the block consists of only one single NaN, then it is marked 3.
        index_vector_first = [];
        index_vector_last = [];
        nan_vector = isnan(y);
        flag = true;
        for i = 1 : size(r_win_e,2)
           if i == 1
             if nan_vector(i) == 0
                flag = false;
             end
           end
           if nan_vector(i) == 0 && flag == true && i == 1
             flag = false;
           elseif nan_vector(i) == 1 && flag == false 
             index_vector_first = [index_vector_first i];
             flag = true;     
           elseif nan_vector(i) == 0 && flag == true 
             index_vector_last = [index_vector_last i-1];
             flag = false;
           end
        end

        % NaN-padding before NaN block.
        for i = 1 : length(index_vector_first)
          j = index_vector_first(i);
          span = j - (w+timespan_diff) - 2;
          gap = floor((w+timespan_diff)/2) + mod(span,ws);
          begin = j - gap;

          for k = 1 : size(r_win_e,1)
             r_win_e(k,begin:j) = NaN;
          end
        end

        % NaN-padding after NaN block.
        for i = 1 : length(index_vector_last)
          j = index_vector_last(i);
          span = j + (w+timespan_diff);
          gap = floor((w+timespan_diff)/2) + mod(span,ws);
          ende = j + gap;

          for k = 1 : size(r_win_e,1)
             r_win_e(k,j:ende) = NaN;
          end
        end       
        
        % bind output
        Z1 = r_win_e;
        
        
    elseif windowshape ~= 0
        
        % preallocate vectors for storing RQA measurements and histograms of 
        % each moving window
        r_win2 = zeros(18,ceil((lenTS-timespan_diff-w_new)/ws_new));
        
        % preallocate cell-array for storing distance matrices
        DM = cell(1,ceil((lenTS-timespan_diff-w_new)/ws_new));
        
        % Controling windowsteps
        startwert2 = 0;  
        cnt2 = 0;
        
        h = waitbar(0,'Performing RQA (running window approach)');
        for i = 1:ws_new:lenTS-(w_new+timespan_diff),waitbar(i/(lenTS-(w_new+timespan_diff)))
               % Define Interval in time series
               yx = y(1+startwert2:w+startwert2+timespan_diff);

               % Embedding
               Y = embed(yx,m,tau);
               % RP
               if strcmp(thres,'fix')
                    [R,P,~] = rp(Y,epsilon,'fix',norm);
               elseif strcmp(thres,'var')
                    [R,P,~] = rp(Y,epsilon,'var',norm);
               elseif strcmp(thres,'fan')
                    [R,P,~] = rp(Y,epsilon,'fan',norm);
               end   
               DM{cnt2+1} = P;

               RP_Ausschnitt=zeros(w_new,w_new);  % leere Matrix mit Kantenlänge w machen. 
                                                   % Hierein soll nun der RP-Ausschnitt
                                                   % kommen. Danach wird diese für die RQA
                                                   % verwendet.

              flag = false;
              maxcounter = round(w_new/2)-1;    % controls the maximum amount of chosen lines of old RP
              counter = 0;                    % controls the amount of chosen lines of old RP

              % Innerhalb dieses Ausschnitts jetzt Werte aus dem RP kopieren, die in 
              % der neuen Fensterform liegen.

            j2=1;                                   % Spaltenindex des RP_Auschnitt
            for j = 1:w_new-1     % Spalten wählen
                % Zeilen wählen
                i2=round(w_new/2)-counter;          % Zeilenindex des RP_Ausschnitt
                for k = (round(w_new/2)-counter):(round(w_new/2)+counter)
                    RP_Ausschnitt(i2,j2)=R(k,j);    % Werte übernehmen
                    i2=i2+1;                        % Zeilenindex des RP_Ausschnitt anpassen
                end

                % wegen rautenform : auf maxcounter achten und dann mit dem
                % zeilenindex, welcher über counter gesteuert wird, wieder
                % heruntergehen.
                if counter == maxcounter
                    flag = true;
                end

                if flag == false
                    counter = counter + 1;
                end

                if flag == true
                    counter = counter - 1;
                end
                j2=j2+1;                % Spaltenindex des RP-Ausschnitts erhöhen
            end

            r_win2(:,cnt2) = rqa2(RP_Ausschnitt,L,T);

           % increase control parameter
           startwert2 = startwert2 + ws_new;
           cnt2 = cnt2 + 1;
           
        end
        close(h)
        
        % Phase correction for windowed measures added by M.H. Trauth by shifting
        % the measures by half the window size and sampled at the resolution ws.

        r_win_w2(1:size(r_win2,1),1:lenTS-timespan_diff) = NaN;
        r_win_w2(:,1+w_new/2:ws_new:lenTS-timespan_diff-w_new/2) = r_win2;
        for i = 1:18
            r_win_w2(i,1+w_new/2:lenTS-timespan_diff-w_new/2) = ...
               fillmissing(r_win_w2(i,1+w_new/2:lenTS-timespan_diff-w_new/2),'linear');
        end

        %
        % Phase correction for embedding delay added by M.H. Trauth by shifting the
        % measures by half the embedding delay. 
        r_win_e2(1:size(r_win2,1),1:lenTS) = NaN;
        r_win_e2(:,1+round(timespan_diff/2):lenTS-floor(timespan_diff/2)) = r_win_w2;
        
        
        % Add gaps in the RQA measures where the data y has gaps.
        for i = 1 : size(r_win_e2,1)
           r_win_e2(i,isnan(y)==1) = NaN;
        end

        % Flag the NaNs by their relative location: The first NaN of a NaN block
        % is assigned the value of 1, the last NaN of the block gets the value 2,
        % and if the block consists of only one single NaN, then it is marked 3.
        index_vector_first = [];
        index_vector_last = [];
        nan_vector = isnan(y);
        flag = true;
        for i = 1 : size(r_win_e2,2)
           if i == 1
             if nan_vector(i) == 0
                flag = false;
             end
           end
           if nan_vector(i) == 0 && flag == true && i == 1
             flag = false;
           elseif nan_vector(i) == 1 && flag == false 
             index_vector_first = [index_vector_first i];
             flag = true;     
           elseif nan_vector(i) == 0 && flag == true 
             index_vector_last = [index_vector_last i-1];
             flag = false;
           end
        end

        % NaN-padding before NaN block.
        for i = 1 : length(index_vector_first)
          j = index_vector_first(i);
          span = j - (w_new+timespan_diff) - 2;
          gap = floor((w_new+timespan_diff)/2) + mod(span,ws_new);
          begin = j - gap;

          for k = 1 : size(r_win_e2,1)
             r_win_e2(k,begin:j) = NaN;
          end
        end

        % NaN-padding after NaN block.
        for i = 1 : length(index_vector_last)
          j = index_vector_last(i);
          span = j + (w_new+timespan_diff);
          gap = floor((w_new+timespan_diff)/2) + mod(span,ws_new);
          ende = j + gap;

          for k = 1 : size(r_win_e2,1)
             r_win_e2(k,j:ende) = NaN;
          end
        end         
        

        % binding output-variables
        Z1=r_win_e2;
        
    end
    
end

%% RQA for the whole RP, no moving window
try
    Y=embed(y,m,tau);

    if strcmp(thres,'fix')
        [R,~,ee]=rp(Y,epsilon,'fix',norm);
    elseif strcmp(thres,'var')
        [R,~,ee]=rp(Y,epsilon,'var',norm);
    elseif strcmp(thres,'fan')
        [R,~,ee]=rp(Y,epsilon,'fan',norm);
    end

    if diagonal_rp
        R = rp_diagonal(R);
    end

    Z2(1:13)=rqa(R,L,T,line_correct);
    % RT max
    Z2(17) = Z2(9);
    % ratio 1
    Z2(9) = Z2(10)/Z2(1);
    % ratio 2
    Z2(14) = Z2(2)/Z2(1);
    % ratio 3
    Z2(15) = Z2(6)/Z2(2);
    % DIV
    Z2(16) = 1/Z2(4);
    % RT-DIV
    Z2(18) = 1/Z2(17);

    % Calculate columnwise RR and store mean of it in 19th column of Z2
    helper = zeros(1,size(R,1));
    for j = 1:size(R,1)
        helper(1,j) = 1/ (sum(R(:,j))/size(R,1));
    %             r_win_w(19,j) = -log((sum(R(:,j))/size(R,1)));
    end
    helper = helper/ max(helper);
    Z2(19) = nanmean(helper);
catch
    warning('No RQA for the entire RP possible, due to the size of the RP')
    Z2 = NaN*ones(19,1);
    ee = NaN;
end
%% Plot Time Series, RP and RQA-Quantity

% Create string array of RQS measures for figure titles
str = [...
"Recurrence rate",...
"Determinism",...
"Mean diagonal line length",...
"Maximal diagonal line length",...
"Entropy",...
"Laminarity",...
"Trapping time",...
"Maximal vertical line length",...
"Ratio RT/RR",...
"Recurrence time of 2nd type",...
"Recurrence time entropy",...
"Clustering coefficient",...
"Transitivity",...
"Ratio DET/RR",...
"Ratio LAM/DET",...
"DIV (1/Lmax)",...
"Maximal white vertical line length",...
"RT-DIV (1/RTmax)",...
"RR_{col}^{-1}",...
"windowed RR_{col}^{-1}"];

% Calculate epsilon value in case of fixed threshold depending on phase
% space diameter:
if strcmp(thres,'var')
    epsilon = round(ee,3);   
end

if windowshape==0 
    % set up first figure, containing the time-series and RQA-Quantity
    figure('Color',[1 1 1],...
    'Position',[150 400 800 600],...
    'Name','TS & spec. RQA-quantity',...
    'NumberTitle','off');

    % time-series-axis
    h(1) = axes('XLim',[x(1,1) x(length(x),1)],...
    'XDir','normal',...
    'YDir','normal',...
    'Box','on',...
    'XGrid','on',...
    'Units','Centimeters',...
    'Position',[2.5 11.5 23 8],...
    'LineWidth',1.2,...
    'FontName','Helvetica',...
    'FontSize',12); hold on

    set(get(h(1),'Title'),...
    'String','Original Time-Series',...
    'FontName','Helvetica',...
    'FontSize',12), hold on

    line(x,y,'LineWidth',1.5)

    % RQA-measurement-axis.
    h(2) = axes('XLim',[x(1,1) x(length(x),1)],...
    'XDir','normal',...
    'YDir','normal',...
    'Box','on',...
    'XGrid','on',...
    'Units','Centimeters',...
    'Position',[2.5 1.5 23 8],...
    'LineWidth',1.2,...
    'FontName','Helvetica',...
    'FontSize',12); hold on

    set(get(h(2),'Title'),...
    'String',str(RQA),...
    'FontName','Helvetica',...
    'FontSize',12), hold on
    set(get(h(2),'XLabel'),...
    'String','TIME',...
    'FontName','Helvetica',...
    'FontSize',12), hold on

    line(x,r_win_e(RQA,:),...
            'LineWidth',2)
    linkaxes(h,'x')

    if runningwindow==0
        % set up second figure, containing the RP and RQA-Quantity
        figure('Color',[1 1 1],...
        'Position',[150 400 600 760],...
        'Name','RP & spec. RQA-quantity',...
        'NumberTitle','off');

        % RP-Axis 
        k(1) = axes('XLim',[x(1,1) x(length(x),1)],...
        'XDir','normal',...
        'YLim',[x(1,1) x(length(x),1)],...
        'YDir','normal',...
        'Box','on',...
        'XGrid','on',...
        'Layer','top',...
        'Units','Centimeters',...
        'Position',[2.5 9 16 16],...
        'LineWidth',1.2,...
        'FontName','Helvetica',...
        'FontSize',12); hold on

        imagesc(x,x,RR); colormap([1 1 1;0 0 0]), hold on

        strx = ['global RR = ',num2str(round(Z2(1),3)*100),'%'];
        text(round(0.05*max(x)),max(x)-0.03*max(x),strx,...
           'BackgroundColor',[1 1 1],...
           'EdgeColor',[1 1 1],...
           'VerticalAlignment','top',...
           'FontSize',15,...
           'FontName','Monaco')
        if strcmp(norm,'max')
            stry = 'Maximum norm';
            text(round(0.05*max(x)),max(x)-0.09*max(x),stry,...
           'BackgroundColor',[1 1 1],...
           'EdgeColor',[1 1 1],...
           'VerticalAlignment','top',...
           'FontSize',15,...
           'FontName','Monaco')
        elseif strcmp(norm,'euc')
            stry = 'Euclidean norm';
            text(round(0.05*max(x)),max(x)-0.09*max(x),stry,...
           'BackgroundColor',[1 1 1],...
           'EdgeColor',[1 1 1],...
           'VerticalAlignment','top',...
           'FontSize',15,...
           'FontName','Monaco')
        end        
        
        set(get(k(1),'Title'),...
        'String',horzcat('        Recurrence Plot \newline \tau = ',num2str(tau),...
        ', m = ',num2str(m), ', \epsilon = ',num2str(epsilon),...
        ' (', num2str(thres),')'),...
        'FontName','Helvetica',...
        'FontSize',12), hold on

        % RQA-measurement-axis.
        k(2) = axes('XLim',[x(1,1) x(length(x),1)],...
        'XDir','normal',...
        'YDir','normal',...
        'Box','on',...
        'XGrid','on',...
        'Units','Centimeters',...
        'Position',[2.5 1.5 16 6],...
        'LineWidth',1.2,...
        'FontName','Helvetica',...
        'FontSize',12); hold on

        set(get(k(2),'Title'),...
        'String',str(RQA),...
        'FontName','Helvetica',...
        'FontSize',12), hold on
        set(get(k(2),'XLabel'),...
        'String','TIME',...
        'FontName','Helvetica',...
        'FontSize',12), hold on

        line(x,r_win_e(RQA,:),...
                'LineWidth',1)
        linkaxes(k,'x')
        set(k(1),'tag','rp')
        set(gcf,'WindowButtonMotionFcn',...
           'h=findobj(gcf,''tag'',''rp'');x = get(h,''xlim''); set(h,''ylim'',x)')
       
       
       
        % set up third figure, containing the RP and time-series
        figure('Color',[1 1 1],...
        'Position',[150 400 600 760],...
        'Name','TS & RP',...
        'NumberTitle','off');

        % RP-Axis 
        j(1) = axes('XLim',[x(1,1) x(length(x),1)],...
        'XDir','normal',...
        'YLim',[x(1,1) x(length(x),1)],...
        'YDir','normal',...
        'Box','on',...
        'Layer','top',...
        'XGrid','on',...
        'XTickMode','auto',...
        'Units','Centimeters',...
        'Position',[2.5 9 16 16],...
        'LineWidth',1.2,...
        'FontName','Helvetica',...
        'FontSize',12); hold on

        imagesc(x,x,RR); colormap(j(1),[1 1 1;0 0 0]), hold on

        strx = ['global RR = ',num2str(round(Z2(1),3)*100),'%'];
        text(round(0.05*max(x)),max(x)-0.03*max(x),strx,...
           'BackgroundColor',[1 1 1],...
           'EdgeColor',[1 1 1],...
           'VerticalAlignment','top',...
           'FontSize',15,...
           'FontName','Monaco')
        if strcmp(norm,'max')
            stry = 'Maximum norm';
            text(round(0.05*max(x)),max(x)-0.09*max(x),stry,...
           'BackgroundColor',[1 1 1],...
           'EdgeColor',[1 1 1],...
           'VerticalAlignment','top',...
           'FontSize',15,...
           'FontName','Monaco')
        elseif strcmp(norm,'euc')
            stry = 'Euclidean norm';
            text(round(0.05*max(x)),max(x)-0.09*max(x),stry,...
           'BackgroundColor',[1 1 1],...
           'EdgeColor',[1 1 1],...
           'VerticalAlignment','top',...
           'FontSize',15,...
           'FontName','Monaco')
        end          
        
        set(get(j(1),'Title'),...
        'String',horzcat('       Recurrence Plot \newline \tau = ',num2str(tau),...
        ', m = ',num2str(m), ', \epsilon = ',num2str(epsilon),...
        ' (', num2str(thres),')'),...
        'FontName','Helvetica',...
        'FontSize',12), hold on

        % TS-axis.
        j(2) = axes('XLim',[x(1,1) x(length(x),1)],...
        'XDir','normal',...
        'YDir','normal',...
        'Box','on',...
        'XGrid','on',...
        'Units','Centimeters',...
        'Position',[2.5 1.5 16 6],...
        'LineWidth',1.2,...
        'FontName','Helvetica',...
        'FontSize',12); hold on

        set(get(j(2),'Title'),...
        'String','Original Time-Series',...
        'FontName','Helvetica',...
        'FontSize',12), hold on
        set(get(j(2),'XLabel'),...
        'String','TIME',...
        'FontName','Helvetica',...
        'FontSize',12), hold on

        line(x,y)
        linkaxes(j,'x')
        set(j(1),'tag','rp')
        set(gcf,'WindowButtonMotionFcn',...
           'j=findobj(gcf,''tag'',''rp'');x = get(j,''xlim''); set(j,''ylim'',x)')
    end
end
  


if windowshape==1 || windowshape==2
    % set up first figure, containing the time-series and RQA-Quantity
    figure('Color',[1 1 1],...
    'Position',[150 400 800 600],...
    'Name','TS & spec. RQA-quantity based on new window-shape',...
    'NumberTitle','off');

    % time-series-axis
    h(1) = axes('XLim',[x(1,1) x(length(x),1)],...
    'XDir','normal',...
    'YDir','normal',...
    'Box','on',...
    'XGrid','on',...
    'Units','Centimeters',...
    'Position',[2.5 11.5 23 8],...
    'LineWidth',1.2,...
    'FontName','Helvetica',...
    'FontSize',12); hold on

    set(get(h(1),'Title'),...
    'String','Original Time-Series',...
    'FontName','Helvetica',...
    'FontSize',12), hold on

    line(h(1),x,y,'LineWidth',1.5)

    % RQA-measurement-axis.
    h(2) = axes('XLim',[x(1,1) x(length(x),1)],...
    'XDir','normal',...
    'YDir','normal',...
    'Box','on',...
    'XGrid','on',...
    'Units','Centimeters',...
    'Position',[2.5 1.5 23 8],...
    'LineWidth',1.2,...
    'FontName','Helvetica',...
    'FontSize',12); hold on

    set(get(h(2),'Title'),...
    'String',horzcat(str(RQA),' new window-shape'),...
    'FontName','Helvetica',...
    'FontSize',12), hold on
    set(get(h(2),'XLabel'),...
    'String','TIME',...
    'FontName','Helvetica',...
    'FontSize',12), hold on

    line(h(2),x,r_win_e2(RQA,:),...
            'LineWidth',2)
    linkaxes(h,'x')
    
    if runningwindow == 0
        % set up second figure, containing the RP and RQA-Quantity
        figure('Color',[1 1 1],...
        'Position',[150 400 600 760],...
        'Name','RP & spec. RQA-quantity based on new window-shape',...
        'NumberTitle','off');

        % RP-Axis 
        h(1) = axes('XLim',[x(1,1) x(length(x),1)],...
        'XDir','normal',...
        'YLim',[x(1,1) x(length(x),1)],...
        'YDir','normal',...
        'Box','on',...
        'XGrid','on',...
        'Layer','top',...
        'Units','Centimeters',...
        'Position',[2.5 9 16 16],...
        'LineWidth',1.2,...
        'FontName','Helvetica',...
        'FontSize',12); hold on

        imagesc(x,x,RR); colormap(h(1),[1 1 1;0 0 0]), hold on

        set(get(h(1),'Title'),...
        'String',horzcat('        Recurrence Plot \newline \tau = ',num2str(tau),...
        ', m = ',num2str(m), ', \epsilon = ',num2str(epsilon),...
        ' (', num2str(thres),')'),...
        'FontName','Helvetica',...
        'FontSize',12), hold on

        % RQA-measurement-axis.
        h(2) = axes('XLim',[x(1,1) x(length(x),1)],...
        'XDir','normal',...
        'YDir','normal',...
        'Box','on',...
        'XGrid','on',...
        'Units','Centimeters',...
        'Position',[2.5 1.5 16 6],...
        'LineWidth',1.2,...
        'FontName','Helvetica',...
        'FontSize',12); hold on

        set(get(h(2),'Title'),...
        'String',horzcat(str(RQA),' new window-shape'),...
        'FontName','Helvetica',...
        'FontSize',12), hold on
        set(get(h(2),'XLabel'),...
        'String','TIME',...
        'FontName','Helvetica',...
        'FontSize',12), hold on

        line(h(2),x,r_win_e2(RQA,:),...
                'LineWidth',1)
        linkaxes(h,'x')
        set(h(1),'tag','rp')
        set(gcf,'WindowButtonMotionFcn',...
           'j=findobj(gcf,''tag'',''rp'');x = get(j,''xlim''); set(j,''ylim'',x)')

        if windowshape==1
            % set up third figure, containing the RP and time-series
            figure('Color',[1 1 1],...
            'Position',[150 400 600 760],...
            'Name','RP & TS',...
            'NumberTitle','off');

            % RP-Axis 
            j(1) = axes('XLim',[x(1,1) x(length(x),1)],...
            'XDir','normal',...
            'YLim',[x(1,1) x(length(x),1)],...
            'YDir','normal',...
            'Box','on',...
            'XGrid','on',...
            'Layer','top',...
            'XTickMode','auto',...
            'Units','Centimeters',...
            'Position',[2.5 9 16 16],...
            'LineWidth',1.2,...
            'FontName','Helvetica',...
            'FontSize',12); hold on

            imagesc(x,x,RR); colormap(j(1),[1 1 1;0 0 0]), hold on

            set(get(j(1),'Title'),...
            'String',horzcat('       Recurrence Plot \newline \tau = ',num2str(tau),...
            ', m = ',num2str(m), ', \epsilon = ',num2str(epsilon),...
            ' (', num2str(thres),')'),...
            'FontName','Helvetica',...
            'FontSize',12), hold on

            % TS-axis.
            j(2) = axes('XLim',[x(1,1) x(length(x),1)],...
            'XDir','normal',...
            'YDir','normal',...
            'Box','on',...
            'XGrid','on',...
            'Units','Centimeters',...
            'Position',[2.5 1.5 16 6],...
            'LineWidth',1.2,...
            'FontName','Helvetica',...
            'FontSize',12); hold on

            set(get(j(2),'Title'),...
            'String','Original Time-Series',...
            'FontName','Helvetica',...
            'FontSize',12), hold on
            set(get(j(2),'XLabel'),...
            'String','TIME',...
            'FontName','Helvetica',...
            'FontSize',12), hold on

            line(x,y)
            linkaxes(j,'x')
            set(j(1),'tag','rp')
            set(gcf,'WindowButtonMotionFcn',...
           'j=findobj(gcf,''tag'',''rp'');x = get(j,''xlim''); set(j,''ylim'',x)')
        end
    end
end


% Additional figures containing all RQA-Measurements

if ShowAll == 1
    
    if windowshape == 0 
        
            figure('Position',[250 110 600 760],'Color',[1 1 1],...
                        'Name','RQA-Quantities 1-6',...
                        'NumberTitle','off');
            % allocate vector storing axis-objects:
            k = zeros(1,6);
            for i = 1 : 6
                if i == 1 || i == 2 || i == 5 || i == 6
                    k(i) = axes('Units','Centimeters',...
                          'Position',[1.5 22-(i-1)*4.2 ...
                          17 3],'XGrid','on',...
                          'XLim',[x(1,1) x(length(x),1)]); hold on
                else
                    k(i) = axes('Units','Centimeters',...
                          'Position',[1.5 22-(i-1)*4.2 ...
                          17 3],'XGrid','on',...
                          'XLim',[x(1,1) x(length(x),1)]); hold on   
                end
                line(x,r_win_e(i,:),...
                    'LineWidth',1)
                title(str(i))
            end
            linkaxes(k,'x')

            figure('Position',[250 110 600 760],'Color',[1 1 1],...
                        'Name','RQA-Quantities 7-12',...
                        'NumberTitle','off');
            index=1;
            % allocate vector storing axis-objects:
            p = zeros(1,6);
            for i = 7:12

                p(index) = axes('Units','Centimeters',...
                      'Position',[1.5 22-(index-1)*4.2 ...
                      17 3],'XGrid','on',...
                      'XLim',[x(1,1) x(length(x),1)]); hold on   
            
                line(x,r_win_e(i,:),...
                    'LineWidth',1)
                title(str(i))
                index = index + 1;
            end
            linkaxes(p,'x')

            figure('Position',[250 110 600 760],'Color',[1 1 1],...
                        'Name','RQA-Quantities 13-18',...
                        'NumberTitle','off');
            index=1;
            % allocate vector storing axis-objects:
            o = zeros(1,6);
            for i = 13:18
                o(index) = axes('Units','Centimeters',...
                      'Position',[1.5 22-(index-1)*4.2 ...
                      17 3],'XGrid','on',...
                      'XLim',[x(1,1) x(length(x),1)]); hold on
                line(x,r_win_e(i,:),...
                    'LineWidth',1)
                title(str(i))
                index = index + 1;
            end
            linkaxes(o,'x')
            
            if runningwindow == 0
                figure('Position',[250 110 600 760],'Color',[1 1 1],...
                            'Name','RQA-Quantity 19',...
                            'NumberTitle','off');
                % allocate vector storing axis-objects:
                l = zeros(1,1);
                index=1;

                for i = 19 : 20
                    l(index) = axes('Units','Centimeters',...
                          'Position',[1.5 22-(index-1)*4.2 ...
                          17 3],'XGrid','on',...
                          'XLim',[x(1,1) x(length(x),1)]); hold on
                    line(x,r_win_e(i,:),...
                        'LineWidth',1)
                    title(str(i))
                    index = index + 1;
                end
                linkaxes(l,'x')
            end

    end

    
    if windowshape == 1 
        
            figure('Position',[250 110 600 760],'Color',[1 1 1],...
                        'Name','RQA-Quantities 1-6, new window-shape',...
                        'NumberTitle','off');
            % allocate vector storing axis-objects:
            k = zeros(1,6);
            for i = 1 : 6
                k(i) = axes('Units','Centimeters',...
                      'Position',[1.5 22-(i-1)*4.2 ...
                      17 3],'XGrid','on',...
                      'XLim',[x(1,1) x(length(x),1)]); hold on
                line(x,r_win_e2(i,:),...
                    'LineWidth',1)
                title(str(i))
            end
            linkaxes(k,'x')
            
            
            figure('Position',[250 110 600 760],'Color',[1 1 1],...
                        'Name','RQA-Quantities 7-12, new window-shape',...
                        'NumberTitle','off');
            index=1;
            % allocate vector storing axis-objects:
            p = zeros(1,6);
            for i = 7:12
                p(index) = axes('Units','Centimeters',...
                      'Position',[1.5 22-(index-1)*4.2 ...
                      17 3],'XGrid','on',...
                      'XLim',[x(1,1) x(length(x),1)]); hold on
                line(x,r_win_e2(i,:),...
                    'LineWidth',1)
                title(str(i))
                index = index + 1;
            end
            linkaxes(p,'x')

            figure('Position',[250 110 600 760],'Color',[1 1 1],...
                        'Name','RQA-Quantities 13-18, new window-shape',...
                        'NumberTitle','off');
            index=1;
            % allocate vector storing axis-objects:
            o = zeros(1,6);
            for i = 13:18
                o(index) = axes('Units','Centimeters',...
                      'Position',[1.5 22-(index-1)*4.2 ...
                      17 3],'XGrid','on',...
                      'XLim',[x(1,1) x(length(x),1)]); hold on
                line(x,r_win_e2(i,:),...
                    'LineWidth',1)
                title(str(i))
                index = index + 1;
            end
            linkaxes(o,'x')
    end    
    
 
end

end



% Helper functions
function y = rqa(varargin)
% RQA performs recurrence quantification analysis
%
%    Y=RQA(X,L,T,line_correct) calculates measures of recurrence 
%    quantification analysis from a recurrence plot X using minimal line 
%    length L and a Theiler window T. Optional input 'line_correct' can
%    be used in order to correct for border effects while counting the
%    diagonal line length histograms (Kraemer & Marwan 2019).
%    If 'line_correct = 1', the "kelo"-correction is applied,
%    if 'line_correct = 0', the conventional method is applied (Default)
%
%    Output:
%      Y(1) = RR     (recurrence rate)
%      Y(2) = DET    (determinism)
%      Y(3) = <L>    (mean diagonal line length)
%      Y(4) = Lmax   (maximal diagonal line length)
%      Y(5) = ENTR   (normalized entropy of the diagonal line lengths)
%      Y(6) = LAM    (laminarity)
%      Y(7) = TT     (trapping time)
%      Y(8) = Vmax   (maximal vertical line length)
%      Y(9) = RTmax  (maximal white vertical line length)
%      Y(10) = T2     (recurrence time of 2nd type)
%      Y(11) = RTE    (recurrence time entropy, i.e., RPDE)
%      Y(12) = Clust  (clustering coefficient)
%      Y(13) = Trans  (transitivity)
% 
%
%    Example:
%         N = 300; % length of time series
%         x = .9*sin((1:N)*2*pi/70); % exemplary time series
%         xVec = embed(x,2,17);
%         R = rp(xVec,.1);
%         Y = rqa(R);
%
% Copyright (c) 2016
% Norbert Marwan, Potsdam Institute for Climate Impact Research, Germany
% http://www.pik-potsdam.de
%
% Modified by Hauke Krämer, Potsdam Institute for Climate Impact Research,
% Germany
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or any later version.

%% check input
narginchk(1,4)
nargoutchk(0,1)

try 
    line_correct = varargin{4};
    if line_correct ~= 1 && line_correct ~= 0
        warning('input "line_correct" must be set to 0 (conventional) or 1 (kelo correction). Now set to 0.')
        line_correct = 0;
    end
catch
    line_correct = 0;
end

try
    theiler_window = varargin{3};
catch
    theiler_window = 1;
end

try
    l_min = varargin{2};
catch
    l_min = 2;
end
x = varargin{1};

% size of recurrence plot
N = size(x);

%% calculation
y = zeros(13,1);

% applt Theiler window to recurrence plot
if theiler_window
   x_theiler = double(triu(x,theiler_window) + tril(x,-theiler_window));
else
   x_theiler = double(x);
end

% reduce the number of possible states by the Theiler window
N_all = N(1)*N(2);
N_all = N_all - N(1) - 2*((theiler_window-1)*N(1) - sum(1:(theiler_window-1)));

% recurrence rate
N_recpoints = sum(x_theiler(:)); % number of rec. points (in complete RP)
y(1) = N_recpoints/N_all; 

% histogram of diagonal lines (look at complete RP)
if line_correct == 0
    [~,l_hist] = dl_conventional(x_theiler);
else
    [~,l_hist] = dl_kelo(x_theiler);
end

% make proper histogram with respect to l_min
l_hist = hist(l_hist,1:N(1));

% determinism
y(2) = sum(l_hist(l_min:N(1)) .* (l_min:N(1))) / N_recpoints;

% mean diagonal line length
if isnan(sum(l_hist(l_min:N(1)) .* (l_min:N(1))) / sum(l_hist(l_min:N(1))))
    y(3) = 0;
else
    y(3) = sum(l_hist(l_min:N(1)) .* (l_min:N(1))) / sum(l_hist(l_min:N(1)));
end

% maximal line length
if any(l_hist)
   y(4) = find(l_hist,1,'last');
end


% diagonal line length entropy
l_hist = l_hist(l_min:end);
% check how many classes occupied in diagonal line histogram
l_classes = sum(l_hist~=0); 

l_prob = l_hist/sum(l_hist); % get probability distribution from histogram
ent_Sum = (l_prob .* log(l_prob));

if l_classes > 1
    y(5) = -nansum(ent_Sum)/log(length(l_min:N(1)));
else
    y(5) = -nansum(ent_Sum);
end

% histogram of vertical lines
v_hist = zeros(1,N(1)); % allocate vector
for i = 1:N(1) % walk along the columns
   cnt = 0;
   for j = 1:N(2)
      if x_theiler(j,i) % are we on a rec. point? (walk along a column)
         cnt = cnt+1; % count number of points on a vertical line
      else % line has ended
         if cnt
             v_hist(cnt) = v_hist(cnt) + 1; % store line length
         end
         cnt = 0; % set back to zero for a new line
      end
   end
   if cnt
       v_hist(cnt) = v_hist(cnt) + 1;
   end
end

% laminarity
y(6) = sum(v_hist(l_min:N(1)) .* (l_min:N(1))) / N_recpoints;

% mean vertical line length (trapping time)
if isnan(sum(v_hist(l_min:N(1)) .* (l_min:N(1))) / sum(v_hist(l_min:N(1))))
    y(7) = 0;
else
    y(7) = sum(v_hist(l_min:N(1)) .* (l_min:N(1))) / sum(v_hist(l_min:N(1)));
end

% maximal vertical length
if any(v_hist)
   y(8) = find(v_hist,1,'last');
end

% recurrence times ("white" vertical lines)
rt_hist = zeros(1,N(1)); % allocate vector
for i = 1:N(1)
   cnt = 0;
   
   % boolean variable to avoid counting white lines at the edges of RP
   first_flag = false;
   
   for j = 1:N(2)
      if ~x(j,i) % are we on a white line?
         if first_flag % line does not cross the RP's edges
            cnt = cnt + 1; % count number of points along the vertical line
         end
      else % we meet a recurrence point
         first_flag = true; % we are for sure within the RP
         if cnt
             rt_hist(cnt) = rt_hist(cnt) + 1; % store line length
         end
         cnt = 0;
      end
   end
end

% maximal white vertical line length
if any(rt_hist)
    y(9) = find(rt_hist,1,'last');
end

% recurrence time
y(10) = sum(rt_hist .* (1:N(1))) / sum(rt_hist);
if isnan(y(10))
    y(10)=0;
end

% recurrence time entropy
rt_classes = sum(rt_hist~=0); % number of occupied bins (for normalization of entropy)
rt_prob = rt_hist/sum(rt_hist); % get probability distribution from histogram
ent_Sum = (rt_prob .* log(rt_prob));
if rt_classes > 1
    y(11) = -nansum(ent_Sum)/log(N(1));
else
    y(11) = -nansum(ent_Sum);
end

% clustering
kv = sum(x_theiler,1); % degree of nodes
y(12) = nanmean(diag(x_theiler*x_theiler*x_theiler)' ./ (kv .* (kv-1)));

% transitivity
denom = sum(sum(x_theiler * x_theiler));
y(13) = trace(x_theiler*x_theiler*x_theiler)/denom;

end

function [a_out, b_out]=dl_kelo(varargin)
% DL_KELO   Mean of the diagonal line lengths and their distribution - 
% Correction for border lines: KEep LOngest diagonal line (kelo)
%    A=dl_kelo(X) computes the mean of the length of the diagonal 
%    line structures in a recurrence plot X. In this correction, just the 
%    longest border line (in each triangle) of the RP is counted. All other 
%    border lines are discarded.
%
%    A=dl_kelo(X,'semi') computes the mean of the length of the diagonal 
%    line structures in a recurrence plot X using the mentionded correction. 
%    Not only lines starting AND ending at a border of the RP, but also semi
%    border lines - lines, that start OR end at a border of the RP - are 
%    denoted as border lines. The longest of these count.
%
%    [A B]=dl_kelo(X,'semi') computes the mean A and the lengths of the
%    found diagonal lines of the recurrence plot X, stored in B, using the 
%    correction mentioned above and also accounts for semi-border diagonals.
%    In order to get the histogramme of the line lengths, simply call 
%    HIST(B,[1 MAX(B)]).
%
%    Examples (CRP toolbox needs to be installed):
%       x = sin(linspace(0,5*2*pi,1050));
%       xe = embed(x,2,50);
%       r = rp(xe,.2);
%       [l l_dist] = dl_kelo(r);
%       subplot(1,2,1)
%       imagesc(r), colormap([1 1 1;0 0 0]), axis xy square
%       title('underlying RP')
%       subplot(1,2,2)
%       histogram(l_dist,1000)
%       xlim([0 1000])
%       xlabel('diagonal line length')
%       ylabel('counts')
%       title('diagonal line length histogram - kelo correction')
%
%
% Copyright (c) 2019-
% K.Hauke Kraemer, Potsdam Institute for Climate Impact Research, Germany
% http://www.pik-potsdam.de
% Institute of Geosciences, University of Potsdam,
% Germany
% http://www.geo.uni-potsdam.de
% hkraemer@pik-potsdam.de, hkraemer@uni-potsdam.de
% Norbert Marwan, Potsdam Institute for Climate Impact Research, Germany
% http://www.pik-potsdam.de
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or any later version.

X = varargin{1};
styleLib={'normal','semi'}; % the possible borderline-style to look for
try
    type = varargin{2};
    if ~isa(type,'char') || ~ismember(type,styleLib)
        warning(['Specified RP type should be one of the following possible values:',...
           10,sprintf('''%s'' ',styleLib{:})])
    end
catch
    type = 'normal';
end

[Y,~] = size(X);
if issymmetric(X)
    [lines(1),borderlines(1)] = getLinesOnDiag(X,-Y+1,type); % init with first (scalar) diagonal
    for j=-Y+2:-1
        [ll,bl] = getLinesOnDiag(X,j,type);
        lines = horzcat(lines,ll);
        borderlines = horzcat(borderlines,bl);
    end
    % append lines for second triangle 
    lines = horzcat(lines,lines);
    % append longest border lines for second triangle (but exclude LOI)
    lines = horzcat(lines,max(borderlines),max(borderlines));
else
    [lines(1),borderlines(1)] = getLinesOnDiag(X,-Y+1,type); % init with first (scalar) diagonal
    for j=-Y+2:Y-1
        [ll,bl] = getLinesOnDiag(X,j,type);
        lines = horzcat(lines,ll);
        borderlines = horzcat(borderlines,bl);
    end
    borderlines = sort(borderlines,'descend');
    % add longest borderlines to lines
    lines = horzcat(lines,borderlines(2:3));
end

% remove lines of length zero (=no line)
zero_lines = lines(:)==0;
lines(zero_lines) = []; 

b_out= sort(lines,'descend')';
a_out = mean(b_out);
end

function [lines, borderline] = getLinesOnDiag(M,j,type)
    d = diag(M,j);
    border_line_length = length(d);
    if ~any(d)
        lines = 0;
        borderline = 0;
        return
    end
    starts = find(diff([0; d],1)==1);
    ends = find(diff([d; 0],1)==-1);

    lines = zeros(1,numel(starts));
    borderline = zeros(1,numel(starts));
    
    if strcmp(type,'normal')
        for n=1:numel(starts)
            if ends(n) - starts(n) + 1 < border_line_length
                lines(n) = ends(n) - starts(n) +1;
            elseif ends(n) - starts(n) + 1 == border_line_length
                borderline = ends(n) - starts(n) +1;
            end
        end
    elseif strcmp(type,'semi')
        for n=1:numel(starts)
            if ends(n) ~= border_line_length && starts(n) ~=1               
                lines(n) = ends(n) - starts(n) +1;
            else
                borderline(n) = ends(n) - starts(n) +1;
            end
        end    
    end
    
end

function [a_out, b_out]=dl_conventional(x)
% DL_CONVENTIONAL   Mean of the diagonal line lengths and their distribution.
%    A=dl_conventional(X) computes the mean of the length of the diagonal 
%    line structures in a recurrence plot.
%
%    [A B]=dl_conventional(X) computes the mean A and the lengths of the
%    found diagonal lines, stored in B. In order to get the 
%    histogramme of the line lengths, simply call 
%    HIST(B,[1 MAX(B)]).
%
%    Examples (CRP toolbox needs to be installed):
%       x = sin(linspace(0,5*2*pi,1050));
%       xe = embed(x,2,50);
%       r = rp(xe,.2);
%       [l l_dist] = dl_conventional(r);
%       subplot(1,2,1)
%       imagesc(r), colormap([1 1 1;0 0 0]), axis xy square
%       title('underlying RP')
%       subplot(1,2,2)
%       histogram(l_dist,1000)
%       xlabel('diagonal line length')
%       ylabel('counts')
%       title('diagonal line length histogram - conventional counting')
%
%    See also CRQA, TT.
%
% Copyright (c) 2008-
% Norbert Marwan, Potsdam Institute for Climate Impact Research, Germany
% http://www.pik-potsdam.de
%
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or any later version.


narginchk(1,1)
nargoutchk(0,2)

warning off
if any(x(:))

  if min(size(x))>100000       % this should speed up the routine; the value
                             % depends on the available memory
    x2=uint8(x);
    N=size(x2);
    x3=zeros(2*N(2)+N(1),N(2));
    x3(N(2)+1:N(2)+N(1),1:N(2))=x2;
    N3=size(x3);
    
    i2=repmat(((1:1+N(2))+N(1)+N(2))',1,N(2));
    i4=i2+repmat((2*N(2)+N(1)+1)*[0:N(2)-1],size(i2,1),1);
    i4(:,end)=[];
    i4=reshape(i4,size(i4,1)*size(i4,2),1);
    x3(i4)=[];
    x3(end)=[];
    x2=(reshape(x3,N(1)+N(2),N(2)))';
  
    x2(end+1,:)=0;
    x=reshape(x2,size(x2,1)*size(x2,2),1);
    x2=x(2:end);x(end)=[];
    z0=find(x==0&x2==1);
    z1=find(x2==0&x==1);
  
  else
  
    N=size(x);
    x1=spdiags(double(x));
    z=reshape(x1,size(x1,1)*size(x1,2),1);
    z2(2:length(z)+1)=z;z2(1)=0;z2(end+1)=0;
    z=diff(z2);
    z0=find(z==1);
    z1=find(z==-1);
  
  end
  if length(z0)>length(z1), z0(end)=[]; end
  if length(z1)>length(z0), z1(end)=[]; end
  
  if isempty(z0), z0=0; end
  if isempty(z1), z1=0; end
  
  if z0(1)>z1(1)
    z0(2:end+1)=z0(1:end);z0(1)=0; 
    if length(z0)>length(z1) 
       z0(end)=[];
    end
  end

  l=sort(z1-z0); %l(end)=[];
  l1=l(find(l-1));
  
  if nargout==2
     b_out=zeros(length(l),1);
     b_out=l';
  end
  
  if nargout>0
     a_out=mean(l1);
  else
     mean(l1)
  end
  
else

  if nargout==2
     b_out=NaN;
  end

  if nargout>0
     a_out=NaN;
  else
     NaN
  end

end

warning on
end

function [X_new,dl_new] = rp_diagonal(varargin)
% RP_DIAGONAL  RP with corrected lines ...
% 
%    [RP_new, dl_new] = rp_diagonal(RP) 
%    computes a new recurrence plot 'RP_new' by altering diagonal line structures 
%    in the input recurrence plot 'RP': Slubs, but also block structures, 
%    are deleted in favour of the longest diagonal lines (skeletonization). 
%    Whenever a diagonal line (starting with the longest lines contained in 
%    the diagonal line length histogram) encounters an adjacent diagonal 
%    line, this adjacent line and - recursively - all its consecutive 
%    adjacent lines, get deleted. 
%
%    Output:
%    You receive the corrected recurrence plot 'RP_new' INCLUDING the line 
%    of identity. Optional you also receive a matrix containing all lines 
%    of this new RP. This matrix has three lines and as many columns as 
%    there are lines in the RP. In the first line the total length of each 
%    line is stored. In the second and third line, the corresponding line 
%    and column indices are stored.
%
%    Example (CRP toolbox needs to be installed):
%      x = sin(linspace(0,5*2*pi,1000));
%      xe = embed(x,2,50);
%      r = rp(xe,.2);
%      [r2, ~] = rp_diagonal(r);
%      figure
%      subplot(1,2,1)
%      imagesc(r), colormap([1 1 1; 0 0 0]), axis xy square
%      title('input RP')
%      subplot(1,2,2)
%      imagesc(r2), colormap([1 1 1; 0 0 0]), axis xy square
%      title('diagonal RP')
%
%   
% Copyright (c) 2019-
% K.Hauke Kraemer, Potsdam Institute for Climate Impact Research, Germany
% http://www.pik-potsdam.de
% Institute of Geosciences, University of Potsdam,
% Germany
% http://www.geo.uni-potsdam.de
% hkraemer@pik-potsdam.de, hkraemer@uni-potsdam.de
% Norbert Marwan, Potsdam Institute for Climate Impact Research, Germany
% http://www.pik-potsdam.de
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or any later version.

%% check input
narginchk(1,3)
nargoutchk(1,2)

X = varargin{1};

% size of the RP
[N_org,M_org] = size(X);
if N_org~=M_org
    error('Input needs to be a squared, binary recurrence matrix')
end
if sum(sum(X>1))~=0 || sum(sum(X<0))~=0 || sum(sum(rem(X,1)))~=0
    error('Input needs to be a squared, binary recurrence matrix')
end

% check whether input RP is symmetric
if issymmetric(X)
    symm = true;
    % if yes, just take the lower triangle
    X2 = tril(X);
    
    % convert this RP into a close returns map and just use the upper half
    X_cl = convertRP(X2);
 
    % get line distributions
    [~, lines_1] = dl_h(X_cl); % black horizontal lines
        
    % make a copy of the line matrix
    lines_1_copy = lines_1;
    
    Nlines = size(lines_1,2); % number of found lines
    
    % create a close returns map with horizontal lines represented by 
    % numbers, equal to its lengths
    X_hori = zeros(size(X_cl));
    for i = 1:size(lines_1,2)
        line_ind = lines_1(2,i);
        column_ind = lines_1(3,i);
        for j = 0:lines_1(1,i)-1
            X_hori(line_ind,column_ind+j) = lines_1(1,i);
        end
    end
      
    
else
    symm = false;
    % if not, store lower triangle in X2 and upper triangle transposed in
    % X3
    X2 = tril(X);
    X3 = triu(X)';
    
    % convert these RPs into close returns maps and just use the upper half
    X_cl = convertRP(X2);
    X_cl2 = convertRP(X3);
    
    % get line distributions
    [~, lines_1] = dl_h(X_cl); % black horizontal lines
    [~, lines_2] = dl_h(X_cl2); % black horizontal lines
    
    % make a copy of the line matrices
    lines_1_copy = lines_1;
    lines_2_copy = lines_2;
    
    Nlines = size(lines_1,2); % number of found lines in lower triangle
    Nlines2 = size(lines_2,2); % number of found lines in upper triangle
    
    % create a close returns map with horizontal lines represented by 
    % numbers, equal to its lengths
    X_hori = zeros(size(X_cl));
    for i = 1:size(lines_1,2)
        line_ind = lines_1(2,i);
        column_ind = lines_1(3,i);
        for j = 0:lines_1(1,i)-1
            X_hori(line_ind,column_ind+j) = lines_1(1,i);
        end
    end
    X_hori2 = zeros(size(X_cl2));
    for i = 1:size(lines_2,2)
        line_ind = lines_2(2,i);
        column_ind = lines_2(3,i);
        for j = 0:lines_2(1,i)-1
            X_hori2(line_ind,column_ind+j) = lines_2(1,i);
        end
    end
end


% scan the lines, start with the longest one and discard all adjacent lines

% initialize final line matrix
line_matrix_final = zeros(3,1);

% go through all lines stored in the sorted line matrix
[N,M] = size(X_hori);
for l_ind = 1:Nlines

    % check if line is still in the rendered line matrix
    if ~ismember(lines_1(:,l_ind)',lines_1_copy','rows')
        continue
    end
    
    % get index pair for start of the line
    linei = lines_1(2,l_ind); 
    columni = lines_1(3,l_ind);
    
    % copy this line in the final line matrix
    line_matrix_final = horzcat(line_matrix_final,lines_1(:,l_ind));
    
    % delete this line from the RP
    X_hori = delete_line_from_RP(X_hori,lines_1(:,l_ind));
 
    % go along each point of the line and check for neighbours
    l_max = lines_1(1,l_ind);
    for l = 1:l_max
        
        % scan each line twice - above and underneth
        for index = -1:2:1
            
            % make sure not to exceed RP-boundaries
            if linei+index > N | linei+index == 0
                break
            end
                           
            % if there is a neighbouring point, call recursive scan-function
            if X_hori(linei+index,columni+l-1)

                [X_hori,lines_1_copy]=scan_lines(X_hori,lines_1_copy,...
                    linei+index,columni+l-1);

            end

        end
        
    end

end

% if not symmetric input RP, than compute for the upper triangle as well
if ~symm
    % initialize final line matrix
    line_matrix_final2 = zeros(3,1);
    for l_ind = 1:Nlines2

        % check if line is still in the rendered line matrix
        if ~ismember(lines_2(:,l_ind)',lines_2_copy','rows')
            continue
        end

        % get index pair for start of the line
        linei = lines_2(2,l_ind); 
        columni = lines_2(3,l_ind);

        % copy this line in the final line matrix
        line_matrix_final2 = horzcat(line_matrix_final2,lines_2(:,l_ind));

        % delete this line from the RP
        X_hori2 = delete_line_from_RP(X_hori2,lines_2(:,l_ind));

        % go along each point of the line and check for neighbours
        l_max = lines_2(1,l_ind);
        for l = 1:l_max

            % scan each line twice - above and underneth
            for scan = 1:2

                if scan == 1
                    index = 1;
                    % make sure not to exceed RP-boundaries
                    if linei+index > N
                        break
                    end
                else
                    index = -1;
                    % make sure not to exceed RP-boundaries
                    if linei+index == 0
                        break
                    end 
                end

                % if there is a neighbouring point, call recursive scan-function
                if X_hori2(linei+index,columni+l-1)

                    [X_hori2,lines_2_copy]=scan_lines(X_hori2,lines_2_copy,...
                        linei+index,columni+l-1);

                end

            end

        end      

    end
end

% build RP based on the histogramm of the reduced lines

X_cl_new = zeros(N,M);
if ~symm
   X_cl2_new = zeros(N,M); 
end

% fill up close returns map with lines stored in the new line matrix
for i = 1:size(line_matrix_final,2)
    l_max = line_matrix_final(1,i);
    linei = line_matrix_final(2,i);
    columni = line_matrix_final(3,i);
    for j = 1:l_max
        X_cl_new(linei,columni+j-1) = 1;
    end
end

if symm
    % revert this close returns map into a legal RP
    XX = revertRP(X_cl_new);    
    X_new = XX + (XX-eye(size(XX)))';
else 
    % fill up close returns map with lines stored in the new line matrix
    for i = 1:size(line_matrix_final2,2)
        l_max = line_matrix_final2(1,i);
        linei = line_matrix_final2(2,i);
        columni = line_matrix_final2(3,i);
        for j = 1:l_max
            X_cl2_new(linei,columni+j-1) = 1;
        end
    end
    % revert this close returns map into a legal RP
    XX = revertRP(X_cl_new);
    XXX= revertRP(X_cl2_new);
    X_new = XX + (XXX-eye(size(XXX)))';
    
end

% bind optional output

% get line distributions of new RP
if nargout > 1
    [~, lines3] = dl_e(X_new);
    dl_new = sortrows(lines3','descend')';
end

end

function [a_out, b_out] = dl_e(X)
% DL_E   Mean of the diagonal line lengths and their distribution
% additionally with the corresponding indices of the lines.
%    A=DL_E(X) computes the mean of the length of the diagonal 
%    line structures in a recurrence plot.
%
%    [A B]=DL_E(X) computes the mean A and the lengths of the
%    found diagonal lines, stored in the first line of B. B is a 3 line
%    matrix storing the found diagonal line lengths in its columns. Line 2
%    and 3 store the indices i, j of the startpoint of the diagonal line
%    stored in the same column in the first line.
%    In order to get the 
%    histogramme of the line lengths, simply call 
%    HIST(B(1,:),[1 MAX(B(1,:))]).
%
%    Examples: X = crp(rand(200,1),1,1,.3,'fan','silent');
%              [l l_dist] = dl_e(X);
%              hist(l_dist(1,:),200)
%
%    See also CRQA, TT, DL.

% Copyright (c) 2019-
% Norbert Marwan, Potsdam Institute for Climate Impact Research, Germany
% http://www.pik-potsdam.de
% K.Hauke Kraemer, Potsdam Institute for Climate Impact Research, Germany
% http://www.pik-potsdam.de
% Institute of Geosciences, University of Potsdam,
% Germany
% http://www.geo.uni-potsdam.de
% hkraemer@pik-potsdam.de, hkraemer@uni-potsdam.de
%
% $Date: 2018/09/19 $
% $Revision:  $
%
% $Log: dl.m,v $
%
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or any later version.
[Y,~] = size(X);
lines(:,1) = getLinesOnDiag(X,-Y+1); % init with first (scalar) diagonal
for j=-Y+2:Y-1
    lines = horzcat(lines,getLinesOnDiag(X,j)); 
end

% remove lines of length zero (=no line)
zero_lines = lines(1,:)==0;
lines(:,zero_lines) = []; 

b_out= sortrows(lines','descend')';
a_out = mean(b_out(1,:));
end

function tuples = get_indices(indices,position_from_LOI)
% GET_INDICES gets "true" indices from of the lines represented in the column 
% vector 'indices' and its diagonal determined by 'position_from_LOI'.
% "True" indices in this case means the indices in the corresponding RP,
% not the close returns map.
%
% tuples = get_indices(indices,position_from_LOI,sourceRP)
%
% Input:
% 'indices' is a column vector, 'position_from_LOI' a integer
%
% Output: a cell array of tuples conating the true line and column index in
% the sourceRP.
%
% Copyright (c) 2019
% K.Hauke Kraemer, Potsdam Institute for Climate Impact Research, Germany
% http://www.pik-potsdam.de
% Institute of Geosciences, University of Potsdam,
% Germany
% http://www.geo.uni-potsdam.de
% hkraemer@pik-potsdam.de, hkraemer@uni-potsdam.de
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or any later version.

%% check input
narginchk(2,2)
nargoutchk(1,1)

if size(indices,1)<size(indices,2)
    indices = indices';
end
if rem(position_from_LOI,1)~=0
    error('position_from_LOI needs to be a integer')
end



%%
tuples = cell(1,length(indices));

for i = 1:length(indices)
    
    if position_from_LOI < 0
        start_line = indices(i) + abs(position_from_LOI);
        start_column = indices(i);
    elseif position_from_LOI > 0
        start_line = indices(i);
        start_column = indices(i) + abs(position_from_LOI); 
    elseif position_from_LOI == 0
        start_line = indices(i);
        start_column = indices(i);  
    end    
   
    tuples{i}=[start_line start_column];
end


end

function [a_out, b_out]=dl_h(x)
% DL_H   Mean of the horizontal line lengths and their distribution in a
% close returns map, additionally with the corresponding indices of the lines.
%    A=DL_H(X) computes the mean of the length of the horizontal 
%    line structures in a close returns map of a recurrence plot.
%
%    [A B]=DL_H(X) computes the mean A and the lengths of the
%    found horizontal lines, stored in the first line of B. B is a 3 line
%    matrix storing the found horizontal line lengths in its columns. Line 2
%    and 3 store the indices i, j of the startpoint of the horizontal line
%    stored in the same column in the first line.
%    In order to get the 
%    histogramme of the line lengths, simply call 
%    HIST(B(1,:),[1 MAX(B(1,:))]).
%
%    Examples: X = crp(rand(200,1),1,1,.3,'fan','silent');
%              [l l_dist] = dl_h(convertRP(X));
%              hist(l_dist(1,:),200)
%
%    See also CRQA, TT, DL, convertRP, revertRP

% Copyright (c) 2018-
% K.Hauke Kraemer, Potsdam Institute for Climate Impact Research, Germany
% http://www.pik-potsdam.de
% Institute of Geosciences, University of Potsdam,
% Germany
% http://www.geo.uni-potsdam.de
% hkraemer@pik-potsdam.de, hkraemer@uni-potsdam.de
% Norbert Marwan, Potsdam Institute for Climate Impact Research, Germany
% http://www.pik-potsdam.de
%
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or any later version.


narginchk(1,1)
nargoutchk(0,2)

[N,~] = size(x);

liness = zeros(3,1);
for j = 1:N
    d = x(j,:)';
    starts = find(diff([0; d],1)==1);
    ends = find(diff([d; 0],1)==-1);
    
    if ~isempty(starts)
        lines = zeros(3,numel(starts));
        for n=1:numel(starts)
            lines(2,n) = j;
            lines(3,n) = starts(n);
            lines(1,n) = ends(n) - starts(n) +1;        
        end
    else
        lines = zeros(3,1);
    end
    liness = horzcat(liness,lines);
end

% remove lines of length zero (=no line)
zero_lines = liness(1,:)==0;
liness(:,zero_lines) = []; 

b_out= flipud(sortrows(liness',1))';
a_out = mean(b_out(1,:));
end

function Y = convertRP(X)
% CONVERTRP   Transforms the standard RP to a close returns map
%    Y = convertRP(X)
%
%    Example:
%      x = sin(linspace(0,5*2*pi,200));
%      xe = embed(x,2,5);
%      r = rp(xe,.2);
%      imagesc(r)
%      c = convertRP(r);
%      imagesc(c)

%% size of RP
N = size(X);

%% initialize new matrix
Y = zeros(2*N(1)+1,N(1));

%% fill rows of Y by the diagonals of X
% upper triangle
for i = 0:N(1)-1
   Y(N(1)+i+1,(1:(N(1)-i))) = diag(X,i);
end
   
% lower triangle
for i = 0:N(1)-1
   Y(N(1)-i+1,(1:(N(1)-i))+i) = diag(X,-i);
end

end

function X = revertRP(Y)
% REVERTRP   Transforms a close returns map to the standard RP
%    X = revertRP(Y)
%
%    Example:
%      x = sin(linspace(0,5*2*pi,200));
%      xe = embed(x,2,5);
%      r = rp(xe,.2);
%      c = convertRP(r);
%      r2 = revertRP(c);
%      imagesc(r2)
%      imagesc(r-r2) % shows the difference between original and reverted

%% size of close returns map
N = size(Y);

%% initialize new matrix
X = zeros(N(2),N(2));

%% make Y to a square matrix, fill the new part with zeros
Z = [Y zeros(N(1),N(2)+1)];
Z = flipud(Z); % flip upside down

%% fill columns of  by the diagonals of Z (but only the first N points) 
for i = 1:N(2)
    di = diag(Z,-i);
    X(:,N(2)-i+1) = di(1:N(2));
end
end

function RP = delete_line_from_RP(RP,l_vec)
% deletes a line, specified in 'l_vec' (line vector, with first line being
% the total line length, the second line the line-index of the starting point
% and the third line the column-index of the starting pint) from the 'RP'.
    
    RP(l_vec(2),l_vec(3)+(1:l_vec(1))-1) = 0;
%    X = RP;
end

function [XX,YY]= scan_lines(XX,l_vec,line,column)
    
    % for the input index tuple look for the start indices
    index = 0;
    while true        
        % check whether the input index tuple is a listed index for starting
        % points of line lengths in the line matrix 
        loc_line = find(line==l_vec(2,:));
        del_ind = loc_line(column+index==l_vec(3,loc_line));
        
        if del_ind
            break
        else
            index = index - 1;
        end
    end   
    
    % delete the line from RP
    %XX = delete_line_from_RP(XX,l_vec(:,del_ind));
    XX(l_vec(2,del_ind),l_vec(3,del_ind)+(1:l_vec(1,del_ind))-1) = 0;
    
    % bind line length, line & column starting index
    len = l_vec(1,del_ind);
    li = l_vec(2,del_ind);
    co = l_vec(3,del_ind);
    
    % delete the line from the line matix
    l_vec(:,del_ind) = [];
    
    [N,M] = size(XX);
    
    %%%%%% check for borders of the RP %%%%%%
    if li-1 < 1
        flag1 = false;
    else
        flag1 = true;
    end
    if li+1 > N
        flag2 = false;
    else
        flag2 = true;
    end
    
    for i = 1:len
        
        % check for borders of the RP
        if li-1 < 1 || co+i-2 == 0
            flag1b = false;
        else
            flag1b = true;
        end
        if li-1 < 1 || co+i > M
            flag1c = false;
        else
            flag1c = true;
        end
        if li+1 > N || co+i-2 == 0
            flag2b = false;
        else
            flag2b = true;
        end
        if li+1 > N || co+i > M
            flag2c = false;
        else
            flag2c = true;
        end
        
        % check above left for a neighbour
        if flag1b && XX(li-1,co+i-2)
            % call function itself
            [XX,l_vec] = scan_lines(XX,l_vec,li-1,co+i-2);
            
        % check above the line for a neighbour
        elseif flag1 && XX(li-1,co+i-1)
            % call function itself
            [XX,l_vec] = scan_lines(XX,l_vec,li-1,co+i-1);
            
        % check above right for a neighbour
        elseif flag1c && XX(li-1,co+i)
            % call function itself
            [XX,l_vec] = scan_lines(XX,l_vec,li-1,co+i);
            
        % check underneeth left for a neighbour    
        elseif flag2b && XX(li+1,co+i-2)
            % call function itself
            [XX,l_vec] = scan_lines(XX,l_vec,li+1,co+i-2);
            
        % check underneeth the line for a neighbour    
        elseif flag2 && XX(li+1,co+i-1)
            % call function itself
            [XX,l_vec] = scan_lines(XX,l_vec,li+1,co+i-1);
            
        % check underneeth right for a neighbour    
        elseif flag2c && XX(li+1,co+i)
            % call function itself
            [XX,l_vec] = scan_lines(XX,l_vec,li+1,co+i);
        end
    end
    
    YY = l_vec;

end

function y = embed(varargin)
%
%   version 1.0
%
%
%     Create embedding vector using time delay embedding
%     Y=EMBED(X,M,T) create the embedding vector Y from the time
%     series X using a time delay embedding with dimension M and
%     delay T. The resulting embedding vector has length N-T*(M-1),
%     where N is the length of the original time series.
%
%     Example:
%         N = 300; % length of time series
%         x = .9*sin((1:N)*2*pi/70); % exemplary time series
%         y = embed(x,2,17);
%         plot(y(:,1),y(:,2))
%
% Copyright (c) 2012
% Norbert Marwan, Potsdam Institute for Climate Impact Research, Germany
% http://www.pik-potsdam.de
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or any later version.

narginchk(1,3)
nargoutchk(0,1)

try
    t = varargin{3};
catch
        t = 1;
end
        
try
    m = varargin{2};
catch
        m = 1;
end    

x = varargin{1}(:);

% length of time series
Nx = length(x);
NX = Nx-t*(m-1);
if t*(m-1) > Nx
   warning('embedding timespan exceeding length of time series')
end

%% create embeeding vector using Matlab's vectorization
% index series considering the time delay and dimension
for mi = 1:m
    jx(1+NX*(mi-1):NX+NX*(mi-1)) = 1+t*(mi-1):NX+t*(mi-1);
end

% the final embedding vector
y = reshape(x(jx),NX,m);
end



