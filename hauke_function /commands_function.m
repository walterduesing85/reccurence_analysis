clear 
clc

close all 

%This script uses the reccurence function of Kraemer et al. 2019

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
  data1 = [dataxrf(:,1) log(dataxrf(:,Vn1)./dataxrf(:,13))];

%apply age model
[data_1,inv] = agemodel_2(ages,data1,650);

%remove outlier
[data_1(:,2)] = filloutliers(data_1(:,2),'pchip','gesd');



%%
m = 2;
tau = 1;
epsilon = 1;
T = 1;
L = 2;
w = 500;
ws = 50;
RQA = 1;


 [Z1,Z2,DM] = rqaplot(data_1(:,1),data_1(:,2),m,tau,epsilon,T,L,w,ws,RQA);
 
 
 
 