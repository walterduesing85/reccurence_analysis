% Script to mark gaps in potassium record
%
% 1. Gaps were detected visually in the potassium record
% 2. Replace x values by NaNs
%
% 12 Sep 2019 - Trauth
  
if gapsoption == 1
    
    gaps = [
     -555.4 -554.3
     -545.7 -542.3
     -343.9 -340.2
     -288.7 -286.3
     -204.1 -202.0
     -185.3 -184.2
     -180.8 -179.8
     -171.0 -167.2
     -121.8 -116.5
     -115.3 -108.1
     -102.3 -101.4
      -91.6  -90.0];
  
    for i = 1 : length(gaps)
        x(t>gaps(i,1) & t<gaps(i,2)) = NaN;
        r_win_e(RQA_Select(1),t>gaps(i,1) & t<gaps(i,2)) = NaN;
        RR(t>gaps(i,1) & t<gaps(i,2),:) = NaN;
        RR(:,t>gaps(i,1) & t<gaps(i,2)) = NaN;
    end
end

