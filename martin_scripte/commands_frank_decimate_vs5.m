%% Script fuer Frank's Daten
% Erstmal alles loeschen
clear, clc, close all

% Print option, 0=no, 1=yes.
printoption1 = 0;
printoption2 = 1;

% Save option, 0=no, 1=yes.
saveoption = 1;

% Scaling option, magsus and colref will be log10 transformed, xrf data are
% Aitchison log10-ratio transformed with ratios by division by Al, 0=no,
% 1=yes.
scaleoption = 1;

% Loading splice depth.
splicedepth = load('data_frank_vs3.txt');

% Loading magnetic susceptibility point measurements.
load datamagsusexport.mat

% Loading XRF data.
% 1	Al	Chew Bahir XRF		
% 2	Ba	Chew Bahir XRF
% 3	Ca	Chew Bahir XRF
% 4	Cl	Chew Bahir XRF
% 5	Fe	Chew Bahir XRF
% 6	K	Chew Bahir XRF
% 7	Mn	Chew Bahir XRF
% 8	Rb	Chew Bahir XRF
% 9	S	Chew Bahir XRF
% 10	Si	Chew Bahir XRF
% 11	Sr	Chew Bahir XRF
% 12	Ti	Chew Bahir XRF
% 13	Zr	Chew Bahir XRF
load dataxrfexport.txt

% Loading color reflectance values [nm].
% 1     2   3   4   5   6   7   8   9   10  11  12  13  14  15  16  17
% Depth 360 370 380 390 400 410 420 430 440 450 460 470 480 490 500 510
%
% 18  19  20  21  22  23  24  25  26  27  28  29  30  31  32  33  34  35
% 520 530 540 550 560 570 580 590 600 610 620 630 640 650 660 670 680 690
%
% 36  37  38  49  40
% 700 710 720 730 740
load datacolrefexport.mat

%%
% Remove duplicate and negative values.
for i = 2 : length(datamagsusexport);
    if datamagsusexport(i,1) == datamagsusexport(i-1,1)
        datamagsusexport(i,:) = NaN;
    end
    
end
datamagsusexport(isnan(datamagsusexport(:,1))==1 | ...
   isnan(datamagsusexport(:,2))==1 ,:) = [];

for i = 2 : length(datacolrefexport);
    if datacolrefexport(i,1) == datacolrefexport(i-1,1)
        datacolrefexport(i,:) = NaN;
    end
end
datacolrefexport(isnan(datacolrefexport(:,1))==1 | ...
   isnan(datacolrefexport(:,2))==1 ,:) = [];

for i = 2 : length(dataxrfexport);
    if dataxrfexport(i,1) == dataxrfexport(i-1,1)
        dataxrfexport(i,:) = NaN;
    end
end
dataxrfexport(isnan(dataxrfexport(:,1))==1,:) = [];
for i = 1 : 15
    dataxrfexport(isnan(dataxrfexport(:,i))==1,i) = 0;
end

%%
% The splicedepth contains duplicate values. We add 0.001 m to the second
% value.
% for i = 1 : length(splicedepth)-1
%     if splicedepth(i) == splicedepth(i+1)
%         i
%         splicedepth(i+1);
%         splicedepth(i+1) = splicedepth(i+1)+0.01;
%     end
% end

%%
% Interpolate
datamagsus(:,1) = splicedepth;
datamagsus(:,2) = interp1(datamagsusexport(:,1),datamagsusexport(:,2),...
               splicedepth,'pchip');
datamagsus(1,2) = 0; % For some reasons the first value is extreme

for i = 1 : 39
    datacolref(:,1) = splicedepth;
    datacolref(:,1+i) = interp1(datacolrefexport(:,1),...
        datacolrefexport(:,1+i),splicedepth,'pchip');
end

datacolref(1,:) = 0; % For some reasons the first value is extreme(1,2) = 0;

for i = 1 : 14
    dataxrf(:,1) = splicedepth;
    dataxrf(:,1+i) = interp1(dataxrfexport(:,1),...
        dataxrfexport(:,1+i),splicedepth,'pchip');
end

dataxrf(1,:) = 0; % For some reasons the first value is extreme(1,2) = 0;

%%
% Calculate blue/red 450 nm/650 nm color reflectance
databluered(:,1) = datacolref(:,1);
databluered(:,2) = datacolref(:,11)./datacolref(:,31);

%%
% Optional display mag sus data.
if printoption1 == 1
figure('Position',[200 600 1800 600],...
       'Color',[1 1 1])
line(datamagsusexport(:,1),datamagsusexport(:,2),...
    'LineWidth',1)
line(datamagsus(:,1),datamagsus(:,2),...
    'LineWidth',1,...
    'Marker','o',...
    'Color',[0.8 0.3 0.1],...
    'MarkerEdgeColor',[0.8 0.3 0.1])
title('Magnetic Susceptibility')
end

% Optional display 360 nm color reflectance data.
if printoption1 == 1
figure('Position',[200 600 1800 600],...
       'Color',[1 1 1])
line(datacolrefexport(:,1),datacolrefexport(:,2),...
    'LineWidth',1)
line(datacolref(:,1),datacolref(:,2),...
    'LineWidth',1,...
    'Marker','o',...
    'Color',[0.8 0.3 0.1],...
    'MarkerEdgeColor',[0.8 0.3 0.1])
title('Color Reflectance Values, 360 nm')
end

% Optional display blue/red 450 nm/650 nm color reflectance
if printoption2 == 1
figure('Position',[200 600 1800 600],...
       'Color',[1 1 1])
line(datacolrefexport(:,1),datacolrefexport(:,11)./datacolrefexport(:,31),...
    'LineWidth',1)
line(databluered(:,1),databluered(:,2),...
    'LineWidth',1,...
    'Marker','o',...
    'Color',[0.8 0.3 0.1],...
    'MarkerEdgeColor',[0.8 0.3 0.1])
title('Color Reflectance Values, 450 nm/650 nm Ratio')
end

% Optional display potassium data.
if printoption1 == 1
figure('Position',[200 600 1800 600],...
       'Color',[1 1 1])
line(dataxrfexport(:,1),dataxrfexport(:,17),...
    'LineWidth',1)
line(dataxrf(:,1),dataxrf(:,17),...
    'LineWidth',1,...
    'Marker','o',...
    'Color',[0.8 0.3 0.1],...
    'MarkerEdgeColor',[0.8 0.3 0.1])
title('XRF Data, Potassium Concentration')
end

%%
% Scaling, magsus and colref will be log10 transformed, xrf data are
% Aitchison log10-ratio transformed with ratios by division by Al, 0=no,
% 1=yes
if scaleoption == 1
    datamagsus(:,2) = log10(datamagsus(:,2));
    datamagsus(isinf(datamagsus)==1) = 0;
    dataxrf(:,2:end) = log10(dataxrf(:,2:end)./dataxrf(:,2));
    dataxrf(isinf(dataxrf)==1) = 0;
    datacolref(:,2:end) = log10(datacolref(:,2:end));
    datacolref(isinf(datacolref)==1) = 0;
    databluered(:,2) = log10(databluered(:,2));
    databluered(isinf(databluered)==1) = 0;
end

%%
% Optional display mag sus data.
if printoption2 == 1
figure('Position',[200 600 1800 600],...
       'Color',[1 1 1])
line(datamagsus(:,1),datamagsus(:,2),...
    'LineWidth',1,...
    'Marker','o',...
    'Color',[0.8 0.3 0.1],...
    'MarkerEdgeColor',[0.8 0.3 0.1])
title('Magnetic Susceptibility')
end

% Optional display 360 nm color reflectance data.
if printoption2 == 1
figure('Position',[200 600 1800 600],...
       'Color',[1 1 1])
line(datacolref(:,1),datacolref(:,2),...
    'LineWidth',1,...
    'Marker','o',...
    'Color',[0.8 0.3 0.1],...
    'MarkerEdgeColor',[0.8 0.3 0.1])
title('Color Reflectance Values, 360 nm')
end

% Optional display blue/red 450 nm/650 nm color reflectance
if printoption2 == 1
figure('Position',[200 600 1800 600],...
       'Color',[1 1 1])
line(databluered(:,1),databluered(:,2),...
    'LineWidth',1,...
    'Marker','o',...
    'Color',[0.8 0.3 0.1],...
    'MarkerEdgeColor',[0.8 0.3 0.1])
title('Color Reflectance Values, 450 nm/650 nm Ratio')
end

% Optional display potassium data.
if printoption2 == 1
figure('Position',[200 600 1800 600],...
       'Color',[1 1 1])
line(dataxrf(:,1),dataxrf(:,7),...
    'LineWidth',1,...
    'Marker','o',...
    'Color',[0.8 0.3 0.1],...
    'MarkerEdgeColor',[0.8 0.3 0.1])
title('XRF Data, Potassium Concentration')
end

%%
% Saving data
if saveoption == 1
    save datacolref.txt datacolref -ascii
    save datamagsus.txt datamagsus -ascii
    save dataxrf.txt dataxrf -ascii
    save databluered.txt databluered -ascii
end





