% Script to read the ICDP Malawi data (Johnson et al. 2016)
%
% Johnson, T.C.; Werne, J.P.; Brown, E.T.; Abbott, A.; Berke, M.; Steinman,
% B.A.; Halbur, J.; Contreras, S.; Grosshuesch, S.; Deino, A.L.; Lyons,
% R.P.; Scholz, C.A.; Schouten, S.; Sinninghe Damste, J.S. (2016): A
% Progressively Wetter Climate In Southern East Africa Over The Past 1.3
% Million Years. doi:10.1038/nature19065 Nature
%
% 1. Reads data from file.
% 6. Interpolates data to evenly spaced time axis to create a new double
%    array |dataicdpmalawiage|:
%
%     1 = age (kyrs BP)
%     2 = calcite (%)
%     3 = d18C (permille)
%     4 = lake level (m)
%
% 7. Creates string array |dataodp967string| with the type of data.
%
% 30 Aug 2019 - Trauth

% Read data from file.
malawi_ca   = load('data_malawi_ca.txt');
malawi_ca(:,1) = malawi_ca(:,1)/1000;
malawi_d13c = load('data_malawi_d13C.txt');
malawi_d13c(:,1) = malawi_d13c(:,1)/1000;
malawi_lakelevel = load('data_malawi_lakelevel.txt');
%malawi_lakelevel(:,1) = malawi_lakelevel(:,1)/1000;

% Create string array for graphics.
dataicdpmalawistring = ["ICDP Malawi Calcite";
                        "ICDP Malawi d13C";
                        "ICDP Malawi Lake Level"];
                    
% Remove rows without data
malawi_ca(isnan(malawi_ca(:,1))==1,:) = [];
malawi_d13c(isnan(malawi_d13c(:,1))==1,:) = [];
malawi_lakelevel(isnan(malawi_lakelevel(:,1))==1,:) = [];

% Interpolating data to evenly spaced time axis.
dataidcpmalawicaage(:,1) = agemodelmin : agemodelres : agemodelmax;
dataidcpmalawicaage(:,2) = interp1(malawi_ca(:,1),...
        malawi_ca(:,2),dataidcpmalawicaage(:,1),inttype);
    
dataidcpmalawid13cage(:,1) = agemodelmin : agemodelres : agemodelmax;
dataidcpmalawid13cage(:,2) = interp1(malawi_d13c(:,1),...
        malawi_d13c(:,2),dataidcpmalawid13cage(:,1),inttype);
    
dataidcpmalawillage(:,1) = agemodelmin : agemodelres : agemodelmax;
dataidcpmalawillage(:,2) = interp1(malawi_lakelevel(:,1),...
        malawi_lakelevel(:,2),dataidcpmalawillage(:,1),inttype);

% Merge
dataicdpmalawiage(:,1) = dataidcpmalawicaage(:,1);
dataicdpmalawiage(:,2) = dataidcpmalawicaage(:,2);
dataicdpmalawiage(:,3) = dataidcpmalawid13cage(:,2);
dataicdpmalawiage(:,4) = dataidcpmalawillage(:,2);

clear malawi*
