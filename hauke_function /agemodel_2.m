function [data_tt,inv] = agemodel_2(ages,data,tmax)
 clear agemodel
 clear datai

%We calculate a new data set with a the maximum resolution of the input
%data to display a wavelt transform and the prox data over time

% We calculate the median sampling rate of the input data set
%We calculate a sampling rate that produces data length- similar to
%original spatial data set

% Remove identical values in the first column.
clear data2b
datai = sortrows(data,1);
[dataii,IA,IC] = unique(datai(:,1));
dataii(:,2) = data(IA,2);
%%
% Remove NaNs in the 2nd column
dataii(isnan(dataii(:,2))==1,:) = [];

data = dataii;
c = max(data(:,1))/max(ages(:,2));
inv = median(diff(data(:,1)))/c;
%%
clear agemodel
clear datai

agemodel(:,1) = 0 : inv : tmax;
agemodel(:,2) = interp1(ages(:,1),ages(:,2),agemodel(:,1),'pchip');
datai(:,1) = interp1(agemodel(:,1),agemodel(:,2),data(:,1),'pchip');
datai(:,2) = data(:,2);
datai(:,3) = data(:,1);


clear data2i_2
data2i_2(:,1) = min(datai(:,1)) : inv : max(datai(:,1));
data2i_2(:,2) = interp1(datai(:,1),datai(:,2),data2i_2(:,1),'pchip');
data2i_2(:,3) = interp1(datai(:,1),datai(:,3),data2i_2(:,1),'pchip');

     
        
count = 1;
                
%   if the interpolated data series is larger then the original data this while loops
%   increases the sampling frequency until the length is lower       
%%        
while length(data2i_2) >= length(data)

clear agemodel
clear datai
    
    
agemodel(:,1) = 0 : inv : tmax;
agemodel(:,2) = interp1(ages(:,1),ages(:,2),agemodel(:,1),'pchip');
datai(:,1) = interp1(agemodel(:,1),agemodel(:,2),data(:,1),'pchip');
datai(:,2) = data(:,2);
datai(:,3) = data(:,1);


clear data2i_2
data2i_2(:,1) = min(datai(:,1)) : inv : max(datai(:,1));
data2i_2(:,2) = interp1(datai(:,1),datai(:,2),data2i_2(:,1),'pchip');
data2i_2(:,3) = interp1(datai(:,1),datai(:,3),data2i_2(:,1),'pchip');
inv = inv*1.001; 


count = count +1;
end

if count > 2
inv = inv*1/1.001
data_tt = data2i_2;
end
end