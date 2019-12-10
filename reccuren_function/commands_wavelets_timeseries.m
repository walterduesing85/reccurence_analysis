
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

xrec = icwt(wt,f,[1/15 1/9],'amor');


% sz = round(1/15*dim(1));%Filter length
sz = 40;
sigma = round(1/5*sz);% sigma
    % length of gaussFilter vector
x = linspace(-sz / 2, sz / 2, sz);
gaussFilter = exp(-x .^ 2 / (2 * sigma ^ 2));
gaussFilter = gaussFilter / sum (gaussFilter); % normalize

y = data_1(:,2);

yfilt =filtfilt(gaussFilter,1,y);%We apply the filter with the function filtfilt




% 
 VariableNames = {'SpliceDepth','Al','Ba','Ca','Cl',...
           'Fe','K','Mn','Rb','S','Si','Sr','Ti','Zr','qualityflag'};

%VariableNames = {'SpliceDepth' 'Smectite Abundance Group' 'Illite Abundance Group' 'Analcime Group' 'Analcime Intensity' 'Smectite Inverted Group'};




inv = mean(diff(data_1(:,1)));
%%

data_1(:,2) = detrend(data_1(:,2),2);

[wt,f_t,coi] = cwt(data_1(:,2),1/inv); 


%%
figure(...
    'Units','Centimeters',...
     'Position',[40 1 40 20],...
    'color',[1 1 1])


axes1 = axes(...
    'Box','off',...
    'Position',[0.1 0.55 0.8 0.4],...
    'LineWidth',1,...
    'FontName','Helvetica');
% 
        pcolor(axes1,data_1(:,1),1./f_t,abs(wt))
    colormap(jet)
%       
          caxis([0 0.5]) %K 17000 Cl 1800 %K 17000 ca 140000 magsus 50
%          
%       
%          
         line(axes1,data_1(:,1),1./coi,'LineWidth',2,'Color','k')
         
%       line(axes1,IG(:,1),IG(:,3)*140,'color','k','LineWidth',2);
%    

% 
 %     ylabel({VariableNames{Vn1} 'period [kyrs]'})
%ylabel({'Ca' 'period [kyrs]'})
% 
%   ylabel('log(k/Ti) [kyrs period]')
% 
 ylabel('Ca [period kyrs]')
 
 colorbar

         ylim([0 120])
         
         
axis = gca;
axis.YColor = 'k';
axis.XGrid = 'on';
axis.YGrid = 'on';
axis.GridAlpha = 1;
axis.XTick = 50:50:650;
axis.YTick = 10:10:120;
axis.FontSize = 16;
axis.GridLineStyle = '--';
axis.XColor = 'k';
axis.Color = 'none';
axis.FontSize = 10;
axis.Layer = 'top';
axis.Box = 'off'
% axis.YColor = 'none';
% axis.XColor = 'none';

   shading interp
         
 xlim([0 620])

axes2 = axes(...
    'Box','off',...
    'Position',[0.1 0.1 0.8 0.4],...
    'LineWidth',1,...
    'FontName','Helvetica');
% 
 yyaxis left
  line(data_1(:,1),yfilt,'Color','k','LineWidth',1.5)
% %line(data_1(:,1),xrec,'Color','k','LineWidth',2,'Color','k')

% axis = gca;
% axis.YDir = 'reverse';
% %ylabel({VariableNames{Vn1} '[counts/s]'})
ylabel('log(Ca/Ti)')
% % ylabel('MagSus')

 %ylabel({'filtered' 'log(Cl(Ti)'})
%  xlabel('age [kyrs BP]')



 yyaxis  right
% % 
% Inso_m = load('data_laskar_inso_0_march.txt');
% 
% Inso_s = load('data_laskar_inso_0_sept.txt');
% 
% l1 = line(Inso_m(:,1),Inso_m(:,2),'LineWidth',1);
% 
% hold on 
% % 
% l2 = line(Inso_s(:,1),Inso_s(:,2),'LineWidth',1,'LineStyle','--');
% 
%legend([l1 l2],'March','September')
% 
% 
% 
%  ylabel({'Equator' 'Inso M/S' '[W/m^2]'})
% % % 


data_11 = load('data_laskar_equator_wet_dry_ratio.txt');
data_22 = load('data_laskar_inso_grad_JD_0.txt');
   prec = [data_11(:,1) data_11(:,2)+data_22(:,2)*0.5];
%   prec_2 = load('data_laskar_inso_0_sept.txt');

 hold on
   line(prec(:,1),prec(:,2),'LineWidth',2)
%    line(prec_2(:,1),prec_2(:,2),'LineWidth',2,'LineStyle','--')
% 

axis = gca;
axis.YDir = 'reverse';

 ylabel('SITIG')


       
axis = gca;
axis.YColor = 'k';
% axis.XGrid = 'on';
% axis.YGrid = 'on';
axis.GridAlpha = 1;
axis.XTick = 50:50:650;

axis.FontSize = 16;
axis.GridLineStyle = '--';
axis.XColor = 'k';
axis.Color = 'none';
axis.FontSize = 10;
axis.Layer = 'top';
axis.Box = 'off'
 xlim([0 620])


%          
%         axis.YColor = 'none';
%         axis.XColor = 'none';
 

  
  
  
  
 

