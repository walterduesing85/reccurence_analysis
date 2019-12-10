% Script to managed the Chew Bahir data for display.
%
%           x       y
%     1   = empty
%  2- 6   = 1,1 and 2,2-6    % PC color reflectance
%  7-20   = 2,1 and 2,2-15   % XRF
% 21-25   = 3,1 and 3,2-6    % XRD
%    26   = 4,1 and 4,2      % Mag sus
% 27-31   = 5,1 and 5,2-6    % Grainsize EM
% 32      = 6,1 and 6,2      % Grainsize mean
% 33-35   = 7,1 and 7,2-4    % TC, TIC, TOC
%
% 5 Sep 2019 - Trauth

if any(varselect>= 35)==1
    msg = 'Please choose varselect=<35';
    error(msg)
end

for i = 1 : 9
if varselect(i)>=2 & varselect(i)<=6
    varselectcol(i,:) = [1,varselect(i)];
elseif varselect(i)>= 7 & varselect(i)<=20
    varselectcol(i,:) = [2,varselect(i)-5];
elseif varselect(i)>=21 & varselect(i)<=25
    varselectcol(i,:) = [3,varselect(i)-19];
elseif varselect(i)==26
    varselectcol(i,:) = [4,varselect(i)-24];
elseif varselect(i)>=27 & varselect(i)<=31
    varselectcol(i,:) = [5,varselect(i)-25];
elseif varselect(i)==32
    varselectcol(i,:) = [6,varselect(i)-30]; 
elseif varselect(i)>=33 & varselect(i)<=35
    varselectcol(i,:) = [7,varselect(i)-31];
end
end

newdata_a1_1 = data{varselectcol(1,1)};
newdata_a1(:,1) = newdata_a1_1(:,1);
newdata_a1(:,2) = newdata_a1_1(:,varselectcol(1,2));
titlestr_a1_1 = datastr{varselectcol(1,1)};

newdata_a2_1 = data{varselectcol(2,1)};
newdata_a2(:,1) = newdata_a2_1(:,1);
newdata_a2(:,2) = newdata_a2_1(:,varselectcol(2,2));

newdata_a3_1 = data{varselectcol(3,1)};
newdata_a3(:,1) = newdata_a3_1(:,1);
newdata_a3(:,2) = newdata_a3_1(:,varselectcol(3,2));

newdata_a4_1 = data{varselectcol(4,1)};
newdata_a4(:,1) = newdata_a4_1(:,1);
newdata_a4(:,2) = newdata_a4_1(:,varselectcol(4,2));

newdata_a5_1 = data{varselectcol(5,1)};
newdata_a5(:,1) = newdata_a5_1(:,1);
newdata_a5(:,2) = newdata_a5_1(:,varselectcol(5,2));

newdata_a6_1 = data{varselectcol(6,1)};
newdata_a6(:,1) = newdata_a6_1(:,1);
newdata_a6(:,2) = newdata_a6_1(:,varselectcol(6,2));

newdata_a7_1 = data{varselectcol(7,1)};
newdata_a7(:,1) = newdata_a7_1(:,1);
newdata_a7(:,2) = newdata_a7_1(:,varselectcol(7,2));

newdata_a8_1 = data{varselectcol(8,1)};
newdata_a8(:,1) = newdata_a8_1(:,1);
newdata_a8(:,2) = newdata_a8_1(:,varselectcol(8,2));

newdata_a9_1 = data{varselectcol(9,1)}; % For wavelet power spectrum
newdata_a9(:,1) = newdata_a9_1(:,1);
newdata_a9(:,2) = newdata_a9_1(:,varselectcol(9,2));

titlestr_a1_1 = datastr{varselectcol(1,1)};
titlestr_a1 = titlestr_a1_1(varselectcol(1,2)-1,:);

titlestr_a2_1 = datastr{varselectcol(2,1)};
titlestr_a2 = titlestr_a2_1(varselectcol(2,2)-1,:);

titlestr_a3_1 = datastr{varselectcol(3,1)};
titlestr_a3 = titlestr_a3_1(varselectcol(3,2)-1,:);

titlestr_a4_1 = datastr{varselectcol(4,1)};
titlestr_a4 = titlestr_a4_1(varselectcol(4,2)-1,:);

titlestr_a5_1 = datastr{varselectcol(5,1)};
titlestr_a5 = titlestr_a5_1(varselectcol(5,2)-1,:);

titlestr_a6_1 = datastr{varselectcol(6,1)};
titlestr_a6 = titlestr_a6_1(varselectcol(6,2)-1,:);

titlestr_a7_1 = datastr{varselectcol(7,1)};
titlestr_a7 = titlestr_a7_1(varselectcol(7,2)-1,:);

titlestr_a8_1 = datastr{varselectcol(8,1)};
titlestr_a8 = titlestr_a8_1(varselectcol(8,2)-1,:);

titlestr_a9_1 = datastr{varselectcol(9,1)}; % For wavelet power spectrum
titlestr_a9 = titlestr_a9_1(varselectcol(9,2)-1,:);

clear i *_1 varselectcol

