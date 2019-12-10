
function [xn,yn,t] = even_space(x,y)
%%
%%This function deterimies the sampling rate of two data sets to eenly
%%space both time sereis of the same time vector.. 

%Input variables : x and y need to have time in the first and y values in
%                   the second column

%Output variables : xn is the correspondinx y values to the input x and y
%                   are the corresponding y-values to the input y. t ios the time vector that
%                   is equal for both output xn and yn

% We determine the shortest data set 

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

%%

if sz >= 10 
sigma = round(1/5*sz);% sigma
    % length of gaussFilter vector
h = linspace(-sz / 2, sz / 2, sz);
gaussFilter = exp(-h .^ 2 / (2 * sigma ^ 2));
gaussFilter = gaussFilter / sum (gaussFilter); % normalize





ys(:,2) =filtfilt(gaussFilter,1,ys(:,2));%We apply the filter with the function filtfilt
end
%we store the data sets in one 
data_b{1} = xs;
data_b{2} = ys;

      
clear data_g
clear data_d

%%

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


%we store the data sets in one 



t = data2i_2(:,1);

if a >1

        xn = data_c(:,1);
        yn = data_c(:,2);

else
    
        xn = data_c(:,2);
        yn = data_c(:,1);
end

end



