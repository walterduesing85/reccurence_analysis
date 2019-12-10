
clear 

clc 
%%


dataxrf = load('data_raw_xrf_11_11_2019.txt');
% % % 
% % % dataxrd = load('data_xrd.txt');
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

%remove outlier
[data_1(:,2)] = filloutliers(data_1(:,2),'pchip','gesd');

[wt,f] = cwt(data_1(:,2),1/inv,'amor');
%%
clear data_12
[data_12(:,2),data_12(:,1)] = resample(data_1(:,2),data_1(:,1),10);



Order = 7;

 Wn = [median(cf_n_2(:,1)), median(cf_n_2(:,2))];
 
                x12 = data_i_2(:,2);
      
                [b12,a12] = butter(Order,Wn/(0.5*(1/inv_d)),'bandpass');     %Filter Generator 
                [h_2,w] = freqz(b12,a12,2^nextpow2(length(x12))); %Filter frequenz response       
                f_2 = (1/inv_d)*w/(2*pi);                  
                xf12 = filtfilt(b12,a12,x12);
                
   



