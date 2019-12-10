function  [yn] = gaussfilter(y,sz)


    
   
sigma = round(1/5*sz);% sigma
    % length of gaussFilter vector
h = linspace(-sz / 2, sz / 2, sz);
gaussFilter = exp(-h .^ 2 / (2 * sigma ^ 2));
gaussFilter = gaussFilter / sum (gaussFilter); % normalize


yn =filtfilt(gaussFilter,1,y);%We apply the filter with the function filtfilt

end


