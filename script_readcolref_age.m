% Script to read the Chew Bahir color reflectance data
%
% 1. Reads data from LacCore File.
% 2. Remove bad values according to quality flag.
% 3. Run PCA on color reflectance values.
% 4. Create double and char array with wavelengths.
% 6. Interpolates data to evenly spaced time axis to create a new double
%    array |datacolrefage|. We only use the first five PCs of the color
%    reflectance values.
%
%     1 = age (kyrs BP)
%     2-6 = PC1 to PC5 of col ref
%
% 7. Creates string array |datacolrefwavelengthlabel|.
% 8. Creates double array |datacolrefwavelength|.
%
% 21 Aug 2019 - Trauth

% Reading data from LacCore File.
fidA = fopen('data_HSPDP_CHB14_XYZ_Oct2017_20180913.txt');
formstrg = [repmat('%s ',1,3),repmat('%u %s ',1,2),...
    '%u %s ',repmat('%f ',1,55)];
A3 = textscan(fidA,formstrg,'Delimiter','\t','Headerlines',3);
fclose(fidA);

fidB = ...
    fopen('data_HSPDP_CHB14_XYZ_Oct2017_20180913_LABELS.txt');
labels = textscan(fidB,'%s');

% Extract color reflectance data from cell array A3, where A3{13}=
% splice depth, A3{25}-A2{63}=color reflectance 360 nm-740 nm and A3{64}=
% Verena's quality flag with 0=outlier; 1=compromised data; 2= good data.
colref(:,1) = A3{13};
colref = sortrows(colref,1);
for i = 1 : 39
    colref(:,1+i) = A3{24+i};
end
colref(:,41) = A3{64};

% Remove bad values according to quality flag.
colref(colref(:,41)==0 | colref(:,41)==1,:) = [];
colref(:,41) = [];

% Run PCA on color reflectance values.
colrefpca(:,1) = colref(:,1);
[coeff,score,~,~,explained] = pca(colref(:,2:40));
colrefpca(:,2:40) = score;
explained(1:5,:);

% Create double and char array with wavelengths.
% datacolrefwavelength = 360 : 10 : 740;
% for i = 1 : length(datacolrefwavelength)
%     datacolrefwavelengthlabel(i,:) = num2str(datacolrefwavelength(:,i));
% end

% Removing duplicate data points.
colrefpca2 = colrefpca;
for i = 1:size(colrefpca2,1)-1
    if colrefpca2(i,1) == colrefpca2(i+1,1)
        colrefpca2(i,1) = NaN;
    end
end
colrefpca2(isnan(colrefpca2(:,1))==1,:) = [];
colrefpca = colrefpca2;
clear colrefpca2

% Replacing NaNs by zeros.
colrefpca(isnan(colrefpca)==1) = 0;

% Interpolating data to age model.
colrefpca2 = colrefpca;
colrefpca2(:,1) = interp1(agemodeltiepoints(:,1),agemodeltiepoints(:,2),...
    colrefpca(:,1),inttype);
colrefpca = colrefpca2;
clear colrefpca2

% Interpolating data to evenly spaced time axis.
datacolrefpcaage(:,1) = agemodelmin : agemodelres : agemodelmax;
for i = 2:40
    datacolrefpcaage(:,i) = interp1(colrefpca(:,1),...
        colrefpca(:,i),datacolrefpcaage(:,1),inttype);
end

% Let us use only the first five PCs of the color reflectance values.
datacolrefpcaage = datacolrefpcaage(:,1:6);

% Create string for graphics.
datacolrefpcastring = ["CHB PC1 Col Ref";
                    "CHB PC2 Col Ref";
                    "CHB PC3 Col Ref";
                    "CHB PC4 Col Ref";                    
                    "CHB PC5 Col Ref"];       

clear colref* A3 explained fid* formstrg
clear i labels score ans coeff