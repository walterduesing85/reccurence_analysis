%% Coherence wavelet script
clc
clear 
close all 

% In this script we use the coherence 
% 
%and 

% 1)    load the data (different time series that resulted from preprocessing script)
% 2)   run wavelet coherence on data sets and plots 
  

%%
%1)

dataxrf = load('data_raw_xrf_11_11_2019.txt');
% % % 
% % % dataxrd = load('data_xrd.txt');
% % % 
  ages = load('data_ages_mubawa_11_11_2019.txt');
Vn1 = 7;
%   data1 = [dataxrf(:,1) dataxrf(:,Vn1)];
% % 
% data1 = load('data_magsus.txt');
% % 
  data1 = [dataxrf(:,1) log(dataxrf(:,Vn1)./dataxrf(:,14))];

%apply age model
[data_1,inv] = agemodel_2(ages,data1,650);

data_2 = load('data_laskar_equator_wet_dry_ratio.txt');

%% we even space the two data sets 
[xn,yn,t] = even_space(data_1,data_2);
%%
ws = 100;
[corr_t,t2,xn2,yn2,w2] = win_corr(data_1,data_2,ws);

%%

%2)

figure

subplot(2,1,1) 

yyaxis right 

plot(t2,xn2)

yyaxis left 

plot(t2,xyn2)

subplot(2,1,2) 




figure

subplot(2,1,1)

wcoherence(xn,yn,years(1));
gca 

title('wavelet coherence [equator Inso grad. CB log(K/Zr)')
xlabel('kyrs BP');
ylabel('period')




subplot(2,1,2)


yyaxis right
laskar = load('data_eccentricity.txt');

plot(laskar(:,1),laskar(:,2))
ylabel('eccentricity [°]')

yyaxis left
plot(t2,corr_t(:,1))

ylabel('spearman corr. coeff.')

axis = gca;
axis.YGrid = 'on';



xlim([0 630])


%%
figure

wcoherence(data_c(:,2),data_c(:,4),years(0.5));
gca 

title('wavelet coherence between magsus CB and and Grant wet dry')
xlabel('kyrs BP');
ylabel('period')

yyaxis right

laskar = load('data_eccentricity.txt');

plot(laskar(:,1),laskar(:,2))


yyaxis left

corr_2 = corr_spearman'

xlim([0 630])

%%



% define window size w (in datapoint units)
w = 48;
% compute the selected window size in years
tf_1 = num2str(t(w));


% preallocate storing matrices for corrcoefficients
corr_corr = zeros(length(w/2+1:length(t)-w/2),2);


% perform running window corrcoeff approach
for i = w/2+1:length(corr_2)-w/2
    % compute rank correlation
    [p1,h1] = corr(corr_2(i-w/2:i+w/2,1),data_ie(i-w/2:i+w/2,2),'Type','Spearman');
    % compute linear correlation
    [p2,h2] = corr(corr_2(i-w/2:i+w/2,1),data_ie(i-w/2:i+w/2,2));
    % store results
    corr_corr(i,1) = p1;
    corr_corr(i,2) = h1;

    
end


corr_corr = corr_corr';

corr_corr(:,1:w/2) = corr_corr(:,w/2+1) * ones(1,w/2);
corr_corr(:,length(data_c)-w/2+1:length(data_c)) = ...
   corr_corr(:,length(data_c)-w/2) * ones(1,w/2);

corr_corr = corr_corr';
%%


figure('Units','centimeters','Position',[40 40 30 30])

subplot(2,1,1)

wcoherence(corr_2(:,1),data_ie(:,2),years(0.5));


title('wavelet coherence correlation and CO_2')
xlabel('kyrs BP');
ylabel('period')


subplot(2,1,2)

yyaxis right

plot(t,corr_2(:,1))

ylabel('corr coeff')

yyaxis left 


plot(t,data_ie(:,2))

ylabel('C02')



xlabel('kyrs BP')
%%

figure('Units','centimeters','Position',[40 40 30 30])

subplot(2,1,1)

wcoherence(corr_2(:,1),data_ie(:,3),years(0.5));


title('wavelet coherence correlation and Ice Volume')
xlabel('kyrs BP');
ylabel('period')


subplot(2,1,2)

yyaxis right

plot(t,corr_2(:,1))

ylabel('corr coeff')

yyaxis left 


plot(t,data_ie(:,3))

ylabel('Ice Volume')



xlabel('kyrs BP')

figure 

yyaxis right

line(t,corr_corr(:,1))

yyaxis left

line(t,corr_corr(:,2))




