%This fuction is desinged to calculate a piecewise correaltion between two
%data sets x and y with the user defined window size ws.
%The function uses the Speaman correlation
%The outputs are corrt_t with the first row containing the correlation
%coefficioents and the second row p values
%  t is the time axis that has been used to evenly sample the two data sets
%  xn = evenly spaced data set according to input x 
%  yn = evenly spaced data set according to input x




function [corr_t,t,xn,yn,w] = win_corr(x,y,ws)

% We determine the shortest data set 

if max(x(:,1)) <= max(y(:,1))
    
    tmax = max(x(:,1));

    
else
    
    tmax = max(y(:,1));

end

%We determine the sampling rate

inv_x = mean(diff(x(:,1)));

inv_y = mean(diff(y(:,1)));
a = 1;


if inv_x >= inv_y
    
    inv = inv_x;
    inv_2 = inv_y;
    xs = x;
    ys = y;
    a = a +1;
    
    
else
    
    inv = inv_y;
    inv_2 = inv_x;
    xs = y;
    ys = x;

    
end

%we determine the starting age
if min(x(:,1)) >= min(y(:,1))
    start = min(x(:,1));
else 
    
    start = min(y(:,1));
end
%% we gauss filter the data set with the higher data set before we downsample it


sz = round(inv/inv_2);
if sz >  10
    
   
sigma = round(1/5*sz);% sigma
    % length of gaussFilter vector
h = linspace(-sz / 2, sz / 2, sz);
gaussFilter = exp(-h .^ 2 / (2 * sigma ^ 2));
gaussFilter = gaussFilter / sum (gaussFilter); % normalize


ys(:,2) =filtfilt(gaussFilter,1,ys(:,2));%We apply the filter with the function filtfilt

end

%%

%we store the data sets in one 
data_b{1} = xs;
data_b{2} = ys;
      
clear data_g
clear data_d

for i = 1 : 2
            % Interpolate data upon an evenly spaced time axes.
            clear data2i_2
            datai = data_b{i};

             % Remove identical values in the first column.
            clear data2b
            datai = sortrows(datai,1);
            [dataii,IA,IC] = unique(datai(:,1));
            dataii(:,2) = datai(IA,2);

            % Remove NaNs in the 2nd column
            dataii(isnan(dataii(:,2))==1,:) = [];

            data2i_2(:,1) = start : inv : tmax;
            data2i_2(:,2) = interp1(dataii(:,1),dataii(:,2),data2i_2(:,1),'pchip');

            % normalization
            data_c(:,i) = (data2i_2(:,2) - mean(data2i_2(:,2)))/std(data2i_2(:,2));
end

        % define window size w (in datapoint units)
        w = round(ws/inv);
        % compute the selected window size in years
        % tf_1 = num2str(t(w));
        %THe window size can only be a even number
        if mod(w,2) == 1

            w = w +1;

        end

        % preallocate storing matrices for corrcoefficients
        corr_spearman = zeros(2,length(w/2+1:length(data_c)-w/2));
    %%  Data padding    
        data_cp = zeros(length(data_c)+ w,2);
        
        data_cp(1:w/2,:) = rand(w/2,2);
        
        data_cp(length(data_c)+w/2 +1:end,:) = rand(w/2,2);
        
        data_cp(w/2+1 : length(data_c) + w/2,:) = data_c;
        


        %%

        % perform running window corrcoeff approach
        for i = w/2+1:length(data_cp)-w/2
            % compute rank correlation
            [p1,h1] = corr(data_cp(i-w/2:i+w/2,1),data_cp(i-w/2:i+w/2,2),'Type','Spearman');
            % store results
            corr_spearman(1,i) = p1;
            corr_spearman(2,i) = h1;


        end

        %%

        % spearman correlation - phase correction accounting for sliding window
        corr_spearman(:,1:w/2) = corr_spearman(:,w/2+1) * ones(1,w/2);
        corr_spearman(:,length(data_cp)-w/2+1:length(data_cp)) = ...
           corr_spearman(:,length(data_cp)-w/2) * ones(1,w/2);
%%

        corr_t = corr_spearman(:,w/2+1 : length(data_c) + w/2)';

        t = data2i_2(:,1);

        if a >1

                xn = data_c(:,1);
                yn = data_c(:,2);

        else

                xn = data_c(:,2);
                yn = data_c(:,1);
        end




end