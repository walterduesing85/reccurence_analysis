% Script to display the wavelet power spectrum of Chew Bahir and other data
%
% 4 Sep 2019 - Trauth
f(9) = figure('Position',[100 100 800 500],...
    'Color',[1 1 1]);

a(9) = axes('Position',[0.1 0.55 0.8 0.3],...
    'Visible','Off');
b(9) = axes('Position',[0.1 0.55 0.8 0.3],...
    'XTickLabels',[],...
    'XGrid','On',...
    'Color','None');
l(9) = line(newdata_a9(:,1),newdata_a9(:,2));

b(10) = axes('Position',[0.1 0.15 0.8 0.3],...
    'XGrid','On',...
    'Color','None');
if wavedisp == 1
    w(9) = pcolor(newdata_a9(:,1),1./fr,abs(wt));
    shading interp, colormap jet, hold on
elseif wavedisp == 2
    contourf(newdata_a9(:,1),1./fr,abs(wt)); 
end
caxis(waveaxis)
set(b(10),...
    'XLim',wavexlim,...
    'YLim',waveylim)
xl(10) = xlabel('Depth (m)');
yl(10) = ylabel('Wavelength (m)');
l(10) = line(newdata_a9(:,1),coi,'Color','w',...
    'LineStyle','--',...
    'LineWidth',2);
a(10) = axes('Position',[0.1 0.15 0.8 0.3],...
    'Visible','Off');
t(10) = text(0.99,0.95,datawavestring,...
    'HorizontalAlignment','Right',...
    'VerticalAlignment','Top');











