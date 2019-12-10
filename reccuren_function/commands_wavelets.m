
clear 

clc 
%%

% 
% dataxrf = load('data_raw_xrf_11_11_2019.txt');

dataxrd = load('data_xrd.txt');

ages = load('data_ages_mubawa_11_11_2019.txt');
Vn1 = 2;
data1 = [dataxrd(:,1) dataxrd(:,Vn1)];

[data_1,inv] = agemodel_2(ages,data1,650);

   [data_1(:,2)] = filloutliers(data_1(:,2),'pchip','gesd');


% VariableNames = {'SpliceDepth','Al','Ba','Ca','Cl',...
%           'Fe','K','Mn','Rb','S','Si','Sr','Ti','Zr','qualityflag'};


VariableNames = {'SpliceDepth' 'Smectite Abundance Group'
                 'Illite Abundance Group'
                  'Analcime Group'
                 'Analcime Intensity'                   
                 'Smectite Inverted Group'};

inv = mean(diff(data_1(:,1)));

[wt,f_t,coi] = cwt(data_1(:,2),1/inv); 
%%


inv_2 = mean(diff(data_2(:,1)));


[wt_2,f_t_2,coi_2] = cwt(data_2(:,2),1/inv_2); 

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

        pcolor(axes1,data_1(:,1),1./f_t,abs(wt))
        colormap(jet)
      
         caxis([0 1800]) %K 17000
         
      
         
         line(axes1,data_1(:,1),1./coi,'LineWidth',2,'Color','k')
         
%       line(axes1,IG(:,1),IG(:,3)*140,'color','k','LineWidth',2);
%    


     ylabel({VariableNames{Vn1} 'period [kyrs]'})
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

   shading interp
         
 xlim([0 620])

axes2 = axes(...
    'Box','off',...
    'Position',[0.1 0.1 0.8 0.4],...
    'LineWidth',1,...
    'FontName','Helvetica');


        pcolor(axes2,data_2(:,1),1./f_t_2,abs(wt_2))
        colormap(jet)
        
         
        shading interp
          caxis([0 0.35]) %K 17000 ca 140000 magsus 50
        hold on
        line(axes2,data_2(:,1),1./coi_2,'LineWidth',2,'Color','k')
%       line(axes2,IG(:,1),IG(:,3)*140,'color','k','LineWidth',2);
       

         xlabel('age [kyrs BP]') 
         ylabel({'log(K/Zr)' 'period [kyrs]'})
%          ylabel({VariableNames{Vn2} 'period [kyrs]'})
         ylim([0 120])
%      
    xlim([0 620])      
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


 
 
 
 
 
 

  
  
  
  
 

