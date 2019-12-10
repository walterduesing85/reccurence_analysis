
clear 

close all

clc 

%we load the xrf data set
dataxrf = load('data_raw_xrf_11_11_2019.txt');
% % % 
% we define the variable names
 VariableNames = {'SpliceDepth','Al','Ba','Ca','Cl',...
           'Fe','K','Mn','Rb','S','Si','Sr','Ti','Zr','qualityflag'};
       
 %%      
%%----SECTION 1 -------------
%%------ We prepare the data for the first Proxy log(Ca/Ti) (data_1)-------


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





% sz = round(1/15*dim(1));%Filter length
sz = 40;
sigma = round(1/5*sz);% sigma
    % length of gaussFilter vector
x = linspace(-sz / 2, sz / 2, sz);
gaussFilter = exp(-x .^ 2 / (2 * sigma ^ 2));
gaussFilter = gaussFilter / sum (gaussFilter); % normalize

y = data_1(:,2);

yfilt =filtfilt(gaussFilter,1,y);%We apply the filter with the function filtfilt

% we calculate the wavelet transforamtion

[wt,f,coi] = cwt(data_1(:,2),1/inv,'amor');

%%
%%----SECTION 2 -------------
%%------ We prepare the data for the second Proxy MS (data_2)-------


% % % 
  ages = load('data_ages_mubawa_11_11_2019.txt');
Vn1 = 7;

data2 = load('data_magsus.txt');
% % 


%apply age model
[data_2,inv_2] = agemodel_2(ages,data2,650);

%remove outlier
[data_2(:,2)] = filloutliers(data_2(:,2),'pchip','gesd');


% sz = round(1/15*dim(1));%Filter length
sz = 40;
sigma = round(1/5*sz);% sigma
    % length of gaussFilter vector
x = linspace(-sz / 2, sz / 2, sz);
gaussFilter = exp(-x .^ 2 / (2 * sigma ^ 2));
gaussFilter = gaussFilter / sum (gaussFilter); % normalize

y = data_2(:,2);

yfilt2 =filtfilt(gaussFilter,1,y);%We apply the filter with the function filtfilt


% we calculate the wavelet transforamtion


[wt2,f2,coi2] = cwt(data_2(:,2),1/inv_2,'amor');

%%
%%----SECTION 3 -------------
%%------ We prepare the data for the second Proxy MS (data_2)-------


% % % 
  ages = load('data_ages_mubawa_11_11_2019.txt');
Vn1 = 7;

data3 = [dataxrf(:,1) log(dataxrf(:,Vn1)./dataxrf(:,14))];


%apply age model
[data_3,inv_3] = agemodel_2(ages,data3,650);

%remove outlier
[data_3(:,2)] = filloutliers(data_3(:,2),'pchip','gesd');


% sz = round(1/15*dim(1));%Filter length
sz = 40;
sigma = round(1/5*sz);% sigma
    % length of gaussFilter vector
x = linspace(-sz / 2, sz / 2, sz);
gaussFilter = exp(-x .^ 2 / (2 * sigma ^ 2));
gaussFilter = gaussFilter / sum (gaussFilter); % normalize

y = data_3(:,2);

yfilt3 =filtfilt(gaussFilter,1,y);%We apply the filter with the function filtfilt


% we calculate the wavelet transformation


[wt3,f3,coi3] = cwt(data_3(:,2),1/inv_3,'amor');


%%


figure(...
    'Units','Centimeters',...
     'Position',[40 1 18 27.7],...
    'color',[1 1 1])



%------Section 1 log(Ca/Ti) --------

axes1 = axes(...
    'Box','off',...
    'Units','Centimeters',...
    'Position',[2 23 14 3.5],...
    'LineWidth',1,...
    'FontName','Helvetica');

 line(data_1(:,1),yfilt,'Color','k','LineWidth',1)
% %line(data_1(:,1),xrec,'Color','k','LineWidth',2,'Color','k')


ylabel('log(Ca/Ti)')

       
axis = gca;
axis.YColor = 'k';
axis.XGrid = 'on';
axis.YGrid = 'on';
axis.GridAlpha = 1;
axis.XTick = 50:50:650;
axis.YTick = -2:1:4;
axis.FontSize = 8;
axis.GridLineStyle = '--';
axis.XColor = 'k';
axis.Color = 'none';
axis.FontSize = 10;
axis.Layer = 'top';
axis.Box = 'off'
axis.YDir = 'reverse';
 xlim([0 620])


%          
%         axis.YColor = 'none';
        axis.XColor = 'none';
 

  
  


axes2 = axes(...
    'Box','off',...
    'Units','Centimeters',...
    'Position',[2 19.5 14 3.5],...
    'LineWidth',1,...
    'FontName','Helvetica');

    pcolor(axes2,data_1(:,1),1./f,abs(wt))
    colormap(jet)
      
    caxis([0 0.7]) % magsus 50

        line(axes2,data_1(:,1),1./coi,'LineWidth',2,'Color','k')
         


         ylabel('log(Ca/Ti) [kyrs period]')
         xlabel('age [kyrs BP]')

         ylim([0 120])
         
         
axis = gca;
axis.YColor = 'k';
axis.XGrid = 'on';
axis.YGrid = 'on';
axis.GridAlpha = 1;
axis.XTick = 50:50:650;
axis.YTick = 10:10:110;
axis.FontSize = 16;
axis.GridLineStyle = '--';
axis.XColor = 'k';
axis.Color = 'none';
axis.FontSize = 8;
axis.Layer = 'top';
axis.Box = 'off'

% axis.YColor = 'none';
% axis.XColor = 'none';

   shading interp
         
 xlim([0 620])

%-----Section 2 MS record--------

axes3 = axes(...
    'Box','off',...
    'Units','Centimeters',...
    'Position',[2 14.5 14 3.5],...
    'LineWidth',1,...
    'FontName','Helvetica');

 line(data_2(:,1),yfilt2,'Color','k','LineWidth',1)


ylabel('MS [?]')

axis = gca;
axis.YColor = 'k';
axis.XGrid = 'on';
axis.YGrid = 'on';
axis.GridAlpha = 1;
axis.XTick = 50:50:650;
axis.YTick = 60:60:450;
axis.FontSize = 8;
axis.GridLineStyle = '--';
axis.XColor = 'k';
axis.Color = 'none';
axis.FontSize = 10;
axis.Layer = 'top';
axis.Box = 'off'
 
xlim([0 620])



%          
%         axis.YColor = 'none';
        axis.XColor = 'none';
 
axes4 = axes(...
    'Box','off',...
    'Units','Centimeters',...
    'Position',[2 11 14 3.5],...
    'LineWidth',1,...
    'FontName','Helvetica');

    pcolor(axes4,data_2(:,1),1./f2,abs(wt2))
    colormap(jet)
      
    caxis([0 50]) % magsus 50

        line(axes4,data_2(:,1),1./coi2,'LineWidth',2,'Color','k')
         


         ylabel('log(Ca/Ti) [kyrs period]')
         xlabel('age [kyrs BP]')


         ylim([0 120])
         
         
axis = gca;
axis.YColor = 'k';
axis.XGrid = 'on';
axis.YGrid = 'on';
axis.GridAlpha = 1;
axis.XTick = 50:50:650;
axis.YTick = 10:10:110;
axis.FontSize = 16;
axis.GridLineStyle = '--';
axis.XColor = 'k';
axis.Color = 'none';
axis.FontSize = 8;
axis.Layer = 'top';
axis.Box = 'off'
% axis.YColor = 'none';
% axis.XColor = 'none';

   shading interp
         
 xlim([0 620])

  
 
%-----Section 3 log(K/Zr) record--------

axes5 = axes(...
    'Box','off',...
    'Units','Centimeters',...
    'Position',[2 6 14 3.5],...
    'LineWidth',1,...
    'FontName','Helvetica');

 line(data_3(:,1),yfilt3,'Color','k','LineWidth',1)


ylabel('log(K/Zr)')

axis = gca;
axis.YColor = 'k';
axis.XGrid = 'on';
axis.YGrid = 'on';
axis.GridAlpha = 1;
axis.XTick = 50:50:650;
axis.YTick = 0.5:0.5:4;
axis.FontSize = 8;
axis.GridLineStyle = '--';
axis.XColor = 'k';
axis.Color = 'none';
axis.FontSize = 10;
axis.Layer = 'top';
axis.Box = 'off'
axis.YDir = 'reverse';
 xlim([0 620])


%          
%         axis.YColor = 'none';
        axis.XColor = 'none';
 
axes6 = axes(...
    'Box','off',...
    'Units','Centimeters',...
    'Position',[2 2.5 14 3.5],...
    'LineWidth',1,...
    'FontName','Helvetica');

    pcolor(axes6,data_3(:,1),1./f3,abs(wt3))
    colormap(jet)
      
    caxis([0 0.5]) % magsus 50

        line(axes6,data_3(:,1),1./coi3,'LineWidth',2,'Color','k')
         


         ylabel('log(Ca/Ti) [kyrs period]')
         xlabel('age [kyrs BP]')


         ylim([0 120])
         
         
axis = gca;
axis.YColor = 'k';
axis.XGrid = 'on';
axis.YGrid = 'on';
axis.GridAlpha = 1;
axis.XTick = 50:50:650;
axis.YTick = 10:10:110;
axis.FontSize = 16;
axis.GridLineStyle = '--';
axis.XColor = 'k';
axis.Color = 'none';
axis.FontSize = 8;
axis.Layer = 'top';
axis.Box = 'off'
% axis.YColor = 'none';
% axis.XColor = 'none';

   shading interp
         
 xlim([0 620])

 

