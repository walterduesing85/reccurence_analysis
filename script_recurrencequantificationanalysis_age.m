% Script to perform the recurrence quantification analysis
%
% 1. Perform RQA
%
% 12 Sep 2019 - Trauth/Kraemer
% 24 Oct 2019 - Trauth
   
% Calculate RQA measures in moving windows.
r_win=zeros(13,ceil((length(R)-w)/ws));  % Preallocate memory
cnt = 1; % Counter
for i = 1:ws:(length(R)-w)
     r_win(:,cnt) = rqa(R(i:(i+w),i:(i+w)),l_min,theiler,line_correct);
     cnt = cnt+1;
end

% Phase correction for windowed measures added by M.H. Trauth by
% shifting the measures by half the window size and sampled at the
% resolution ws. To draw a continuous line the phase corrected measures
% are interpolated linearly using fillmissing.
r_win_w(1:size(r_win,1),1:size(R,1)) = NaN;
r_win_w(:,1+w/2:ws:size(R,1)-w/2) = r_win;
for i = 1 : 13
     r_win_w(i,1+w/2:size(R,1)-w/2) = ...
     fillmissing(r_win_w(i,1+w/2:size(R,1)-w/2),'linear');
end

% Phase correction for embedding delay added by M.H. Trauth by shifting 
% the measures by half the embedding delay. To draw a continuous line
% the phase corrected measures are interpolated linearly using
% fillmissing.
r_win_e(1:size(r_win,1),1:size(x,2)) = NaN;
r_win_e(:,1+round(timespan_diff/2): ...
 size(x,2)-floor(timespan_diff/2)) = r_win_w;

% Add gaps in the RQA measures where the data x have gaps.
for i = 1 : 13
   r_win_e(i,isnan(x)==1) = NaN;
end

% Flag the NaNs by their relative location: The first NaN of an NaN block
% is assigned the value of 1, the last NaN of the block gets the value 2,
% and if the block consists of only one single NaN, then marked it by 3.
index_vector_first = [];
index_vector_last = [];
nan_vector = isnan(x);
flag = true;
for i = 1 : size(r_win_e,2)
   if i == 1
     if nan_vector(i) == 0
        flag = false;
     end
   end
   if nan_vector(i) == 0 && flag == true && i == 1
     flag = false;
   elseif nan_vector(i) == 1 && flag == false 
     index_vector_first = [index_vector_first i];
     flag = true;     
   elseif nan_vector(i) == 0 && flag == true 
     index_vector_last = [index_vector_last i-1];
     flag = false;
   end
end

% Correct the RQA values for gaps. NaN-padding before NaN block.
for i = 1 : length(index_vector_first)
     j = index_vector_first(i);
     span = j - (w) - 2;
     gap = (w/2) + mod(span,ws);
     begin = j - gap;    
     for k = 1 : 13
       r_win_e(k,begin:j) = NaN;
     end
end

for i = 1 : length(index_vector_last)
     j = index_vector_last(i);
     span = j + (w);
     gap = (w/2) + mod(span,ws);
     ende = j + gap;    
     for k = 1 : 13
        r_win_e(k,j:ende) = NaN;
     end
end

  