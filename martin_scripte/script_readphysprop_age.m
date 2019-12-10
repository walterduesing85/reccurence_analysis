% Script to read the Chew Bahir Physical Properties Data
%
% 1. Reads data from LacCore File.
% 2. Removes bad points according to Verena's quality flags 0 and 1. The
%    column with quality flags is replaced by a column of ones after the 
%    flags were used to remove bad points.
% 4. Removes duplicate data points.
% 5. Replaces NaNs by zeros.
% 6. Interpolates data to evenly spaced time axis to create a new double
%    array |dataphyspropage|:
%
%     1 = age (kyrs BP)
%     2 = magnetic susceptibility
%
% 7. Creates string array |dataphyspropstring| with element names.
%
% 21 Aug 2019 - Trauth

% Reading data from LacCore File.
fidA = fopen('data_HSPDP_CHB14_XYZ_Oct2017_20180913.txt');
formstrg = [repmat('%s ',1,3),repmat('%u %s ',1,2),...
    '%u %s ',repmat('%f ',1,55)];
A2 = textscan(fidA,formstrg,'Delimiter','\t','Headerlines',3);
fclose(fidA);
clear fidA

clear fidB ans formstrg i

% Extract magnetic susceptibility data from cell array A2, where A2{13}=
% splice depth, A2{17}=magsus and A2{64}=Verena's quality flag with 0=
% outlier; 1=compromised data; 2= good data.
dataphysprop(:,1) = A2{13};  % Splice depth
dataphysprop(:,2) = A2{17};  % Magnetic susceptibility
dataphysprop(:,3) = A2{64};  % Quality flags
dataphysprop = sortrows(dataphysprop,1);

% Remove bad values according to quality flag.
dataphysprop2 = dataphysprop;
dataphysprop2(dataphysprop2(:,3)==0 | dataphysprop2(:,3)==1,:) = [];
dataphysprop2(:,3) = [];

dataphysprop = dataphysprop2;
clear dataphysprop2

% Removing duplicate data points.
dataphysprop2 = dataphysprop;
for i = 1:size(dataphysprop2,1)-1
    if dataphysprop2(i,1) == dataphysprop2(i+1,1)
        dataphysprop2(i,1) = NaN;
    end
end
dataphysprop2(isnan(dataphysprop2(:,1))==1,:) = [];
dataphysprop = dataphysprop2;
clear dataphysprop2

% Replacing NaNs by zeros.
dataphysprop(isnan(dataphysprop)==1) = 0;

% Interpolating data to age model.
dataphysprop2 = dataphysprop;
dataphysprop2(:,1) = interp1(agemodeltiepoints(:,1),agemodeltiepoints(:,2),...
    dataphysprop(:,1),inttype);
dataphysprop = dataphysprop2;
clear dataphysprop2

% Interpolating data to evenly spaced time axis.
dataphyspropage(:,1) = agemodelmin : agemodelres : agemodelmax;
for i = 2:2
    dataphyspropage(:,i) = interp1(dataphysprop(:,1),...
        dataphysprop(:,i),dataphyspropage(:,1),inttype);
end
dataphyspropage(:,48) = 1;

% Create string for graphics.
dataphyspropstring = "CHB Mag Sus";

clear i dataphysprop A2 labels






