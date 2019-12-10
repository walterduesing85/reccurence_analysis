% Script to read the Chew Bahir stratigraphic column from Verena
%
% 1. Reads data from Verena's file.
% 2. Creates image array |stratcol|.
% 3. Creates x and y for the image.
%
% 27 Aug 2019 - Trauth

stratcol = imread('data_datastratcolumn.png');
stratcolx = 1/size(stratcol,2) : 1/size(stratcol,2) : 1;
stratcoly = 0.28 : ...
    (292.873-0.28)/(size(stratcol,1)-1) : 292.873;
