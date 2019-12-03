%% Script fuer Frank's Daten
% 11 Sep 2019 Trauth, latest age model
clear, clc, close all

% Loading data, previously known as data_frank_vs3.txt.
splicedepth = load('data_frank_grainsize_depths.txt');

%%
%    1 = mubawa_best_1 (Duesing's best tuned age model)
%    2 = oxcal_2 (Roberts+Ramsey R+R May 2019 age model)
%    3 = martin_550_1 (Trauth's 550 kyr tuned age model)
%    4 = martin_560_1 (Trauth's 560 kyr tuned age model)
%    5 = martin_570_1 (Trauth's 570 kyr tuned age model)
%    6 = martin_merge_CBR+HR_MHT_500 (composite of R+R and Trauth's 500)
%    7 = martin_merge_CBR+HR_MHT_550 (composite of R+R and Trauth's 550)
%    8 = martin_merge_CBR+HR_MHT_570 (composite of R+R and Trauth's 570)
for agemodeloption = 1 : 8
if agemodeloption == 1
    agemodeltiepoints = load('agemodeltiepoints_mubawa_best_1.txt');
elseif agemodeloption == 2
    agemodeltiepoints = load('agemodeltiepoints_oxcal_2.txt');
    agemodeltiepoints(:,2) = agemodeltiepoints(:,2)/1000;
elseif agemodeloption == 3
    agemodeltiepoints = load('agemodeltiepoints_martin_550_1.txt');
elseif agemodeloption == 4
    agemodeltiepoints = load('agemodeltiepoints_martin_560_1.txt');
elseif agemodeloption == 5
    agemodeltiepoints = load('agemodeltiepoints_martin_570_1.txt');
elseif agemodeloption == 6
    agemodeltiepoints = load('agemodeltiepoints_martin_merge_CBR+HR_MHT_500.txt'); 
elseif agemodeloption == 7
    agemodeltiepoints = load('agemodeltiepoints_martin_merge_CBR+HR_MHT_550.txt');   
elseif agemodeloption == 8
    agemodeltiepoints = load('agemodeltiepoints_martin_merge_CBR+HR_MHT_570.txt');   
end
ageinterp(:,agemodeloption) = interp1(agemodeltiepoints(:,1),...
    agemodeltiepoints(:,2),...
    splicedepth,'pchip');
end

save data_frank_agemodels.txt ageinterp -ascii






