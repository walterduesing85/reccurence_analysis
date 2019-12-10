
clear

close all 

clc



dataxrf = load('data_raw_xrf_11_11_2019.txt');
% % % 
% % % dataxrd = load('data_xrd.txt');
% % % 
%%
figure 

subplot(4,2,1)

histogram(dataxrf(:,4))
title('Ca')
ylabel('counts Ca')
xlabel('counts/s')

xlim([0 1000000])

subplot(4,2,2) 

histogram(dataxrf(:,7))


title('K')
ylabel('counts K')
xlabel('counts/s')


subplot(4,2,3) 

histogram(dataxrf(:,13))

title('Ti')

ylabel('counts Ti')
xlabel('counts/s')

subplot(4,2,4) 

histogram(dataxrf(:,14))

title('Zr')

ylabel('counts Zr')
xlabel('counts/s')

subplot(4,2,5) 

histogram(log(dataxrf(:,4)./dataxrf(:,13)))

title('log(Ca/Ti)')

ylabel('counts log(Ca/Ti)')
xlabel('counts/s')


subplot(4,2,6) 

histogram(log(dataxrf(:,7)./dataxrf(:,14)))

title('log(K/Zr)')

ylabel('counts log(K/Zr)')
xlabel('counts/s')


subplot(4,2,7) 

scatter(dataxrf(:,4),dataxrf(:,13),0.5)

title('Ca vs Ti')

ylabel('Ti')
xlabel('Ca')

subplot(4,2,8) 

scatter(dataxrf(:,7),dataxrf(:,14),0.5)

title('K vs Zr')

ylabel('Zr')
xlabel('K')
%%

figure

sz = 20

x = normalize(dataxrf(:,4));


[xn] = gaussfilter(x,sz)
y = normalize(dataxrf(:,13));

[yn] = gaussfilter(y,sz)
scatter(x,y,0.5)

title('Ca vs Ti')

ylabel('Ti')
xlabel('Ca')
%%
figure

a = log(dataxrf(:,7)./dataxrf(:,14));
b = log(dataxrf(:,4)./dataxrf(:,13));



scatter(a,b,0.5)

%%
figure

subplot(2,1,1)

d = dataxrf(:,14);
c = load('data_magsus.txt');

histogram(c)

% c = normalize(c)
subplot(2,1,2)

scatter(d,c(1:length(d),2),0.5)
































