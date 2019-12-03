% Script to compute the recurrence plot
%
% 1. Create embedding vector.
% 2. Compute recurrence plot.
% 3. Phase correction
%
% 12 Sep 2019 - Trauth/Kraemer
% 17 Sep 2019 - Kraemer
% 24 Oct 2019 - Trauth


% Starting value e for threshold calculation, examples are e = 1 for fix,
% e = 0.06 for var and e = 0.08 for fan.
if strcmp(threshold_calculation,'var')
   e = 0.08;
elseif strcmp(threshold_calculation,'fan')
   e = 0.08;
elseif strcmp(threshold_calculation,'fix')
   e = 0.1*range(x);
end

% Create embedding vector.
xVec = embed(x,m,tau);

% Compute normal RP.
[RP,~,eps] = rp(xVec,e,threshold_calculation,norm);

% Compute diagonal RP if needed.
if diagonal == 1
    R = rp_diagonal(RP);
else
    R = RP;  
end

% Phase correction of RP regarding embedding parameters.
RR = NaN(size(x,2),size(x,2));
RR(1+round(timespan_diff/2):size(x,2)-floor(timespan_diff/2),...
   1+round(timespan_diff/2):size(x,2)-floor(timespan_diff/2)) = R;
