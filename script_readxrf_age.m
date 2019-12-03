% Script to read the Chew Bahir XRF Data
%
% 1. Reads revised (Aug 2019) XRF data from LacCore File.
% 2. Removes bad points according to Verena's quality flags 0 and 1. The
%    column with quality flags is replaced by a column of ones after the 
%    flags were used to remove bad points. The columns of ones is used to
%    display elements instead of element ratios.
% 3. Removes duplicate data points.
% 4. Replaces NaNs by zeros.
% 5. Interpolates data to evenly spaced time axis to create a new double
%    array |dataxrfage|:
%
%     1 = age (kyrs BP)
%     2-47 = chemical elements
%     48 = ones
%
% 6. Creates string array |dataxrflabels| with element names.
%
% 26 Aug 2019 - Trauth

% Reading data from LacCore File.
dataall = datastore('data_HSPDP_CHB14_XRF_Aug2019_20190827.txt');
dataall.SelectedVariableNames = {'SpliceDepth','Al','Ba','Ca','Cl',...
          'Fe','K','Mn','Rb','S','Si','Sr','Ti','Zr','qualityflag'};
dataxrftable = readall(dataall);

% Creating labels.
for i = 1 : 14
   lbs(:,i) = ["CHB ";string(dataall.SelectedVariableNames(:,i))];
end
lbs = lbs';
dataxrfstring = strcat(lbs(:,1),lbs(:,2));
dataxrfstring(1,:) = [];

% Converting table |dataxrftable| to double |dataxrf|.
dataxrf = table2array(dataxrftable);

% Removing bad points according to Verena's quality flags 0 and 1.
dataxrf(dataxrf(:,15) == 0 | dataxrf(:,15) ==1,:) = [];

% Removing duplicate data points.
dataxrf2 = dataxrf;
for i = 1:size(dataxrf2,1)-1
    if dataxrf2(i,1) == dataxrf2(i+1,1)
        dataxrf2(i,1) = NaN;
    end
end
dataxrf2(isnan(dataxrf2(:,1))==1,:) = [];
dataxrf = dataxrf2;
clear dataxrf2

% Removing gaps.
dataxrf(isnan(dataxrf(:,2))==1,:) = [];

% Surpress warning concerning the ignorance of NaNs during interpolation.
warning('off','MATLAB:interp1:NaNstrip')

% Interpolating data to age model.
dataxrf2 = dataxrf;
dataxrf2(:,1) = interp1(agemodeltiepoints(:,1),agemodeltiepoints(:,2),...
    dataxrf(:,1),inttype);
dataxrf = dataxrf2;
clear dataxrf2

% FOR FRANK: Export data interpolated to age model but not evenly spaced
frankfilename = strcat("data_frank_xrf_",...
    agemodelstring(agemodeloption),".txt");
save(frankfilename,'dataxrf','-ascii')

% Interpolating data to evenly spaced time axis.
dataxrfage(:,1) = agemodelmin : agemodelres : agemodelmax;
for i = 2:14
    dataxrfage(:,i) = interp1(dataxrf(:,1),...
        dataxrf(:,i),dataxrfage(:,1),inttype);
end
dataxrfage(:,15) = 1;

clear i dataall dataxrftable dataxrf lbs




