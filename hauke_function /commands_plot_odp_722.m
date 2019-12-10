%%
clear
% close all

data_11 = load('data_laskar_equator_wet_dry_ratio.txt');
data_22 = load('data_laskar_inso_grad_JD_0.txt');


% figure
% 
% yyaxis right 
% 
% plot(data_22(:,1),data_22(:,2))
% 
% yyaxis left 
% hold on 

%     
  data_1 = [data_11(:,1) data_22(:,2)./data_11(:,2)];

line(data_1(:,1),data_1(:,2))


%%

% 
% inv = mean(diff(data_11(:,1)));
% 
% [data_1(:,2),data_1(:,1)] = resample(data_11(:,2),data_11(:,1),inv);

[wt,f_t,coi] = cwt(data_1(:,2),1); 







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
          caxis([0 1]) %K 17000 Cl 1800 %K 17000 ca 140000 magsus 50
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
 ylabel('prec. - evap. 0° period kyrs]')
 
 colorbar

        
         
         
axis = gca;
axis.YColor = 'k';
axis.XGrid = 'on';
axis.YGrid = 'on';
axis.GridAlpha = 1;
axis.XTick = 50:50:650;
axis.YTick = 10:10:100;
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
         
 xlim([0 600])
ylim([0 120])
axes2 = axes(...
    'Box','off',...
    'Position',[0.1 0.1 0.8 0.4],...
    'LineWidth',1,...
    'FontName','Helvetica');
% 
%  yyaxis left
  line(data_1(:,1),data_1(:,2),'Color','k','LineWidth',1.5)
% %line(data_1(:,1),xrec,'Color','k','LineWidth',2,'Color','k')

% axis = gca;
% axis.YDir = 'reverse';
% %ylabel({VariableNames{Vn1} '[counts/s]'})
 ylabel('prec. - evap. 0°')

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
 xlim([0 600])
% ylim([20 30])

%          
%         axis.YColor = 'none';
%         axis.XColor = 'none';
 

  
  
  
  
 

