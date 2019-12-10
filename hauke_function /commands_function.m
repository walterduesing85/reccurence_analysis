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
% zeit vektor:
 t = data_1(:,1);
 % data vektor:
 y = data_1(:,2);
 % embedding dimension
 m = 2;
 % time delay
 tau = 1;
 % threshold value
 threshold = 0.05;
 % Theiler window
 T = 1;
 % minimum line length
 l_min = 2;
 % window size
 w = 100;
 % window step
 ws = 10;
 % RQA quantifier you'd like to enlarge in a separate window
 RQA = 2; % 2 = DET
 % norm
 norm = 'euc';
 % threshold selection method
 threshold_meth = 'var';
 % window shape of running window
 window_shape = 0;
 % diagonal RP
 diagonal_RP = 0;
 % border line correction
 line_correct = 0;
 % running window over a global RP
 running_window = 0;
 % Show all RQA quantifiers
 ShowAll = 1;


 [Z1,Z2,DM] = rqaplot(t,y,m,tau,threshold,T,l_min,w,ws,RQA,norm,...
                      threshold_meth,window_shape,diagonal_RP,...
                      line_correct,running_window,ShowAll);
 
 
 
 