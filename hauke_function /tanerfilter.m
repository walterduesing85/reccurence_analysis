
% tanerfilter.m
function [tanerbandx,filtout,f]=tanerfilter(x,dt,fc,fl,fh,c)
%            Produces a bandpassed version of an input real-
%            valued time series by FFT, multiplication by filter
%            designed in frequency domain, and inverse FFT.
%            Has zero phase response.
%
%    INPUTS:
%           x = 1 column vector of data to filter
%           dt = sample rate of x 
%           fc = center frequency of passband
%           f1 = lower cut-off frequency 
%           fh = upper cut-off frequency
%           c = roll-off/octave (default=12) (>12 is steeper)
%           tanerbandx = output filtered series
%
%    OUTPUTS:
%     tanerbandx = 1 column vector with filtered data
%     filtout = zero and positive frequencies of filter gain
%     f = frequency-scale for filter
%
%    PROCEDURE:
%  1. Prepare data length to next power of 2 + npad*zero-padding 
%     (default npad=10)
%  2. Set up filter parameters 
%  4. Set up filter across Fourier domain
%  5. Apply filter to FFT of zero-padded data
%  6. Inverse FFT filtered data and remove zero-padding
%       
% EXAMPLE:
% dt=5.;
% len=2000;
% t=(1:dt:len*dt)';
% x=randn(len,1)+cos(2*pi*t/20);
% fc=1/20;
% fl=fc-1/50;
% fh=fc+1/50;
% [tx,tfilter,tf]=tanerfilter(x,dt,fc,fl,fh);
% figure;plot(tf,tfilter);xlim([0,0.1]);
% figure;plot(t,x,t,tx+4);xlim([0 500]);legend('x','tx');
%
%      Reference:
%        Taner, M.T. (2000), Attributes revisited, Technical Publication, 
%         Rock Solid Images, Inc., Houston, Texas, URL:
%         http://www.rocksolidimages.com/pdf/attrib_revisited.htm 
%
%   Brief editing log:
%   L. Hinnov, October, 2016 - Updated with internal zero-padding 
%                              procedure to avoid problems with odd
%                              vs.even FFTs. (Suggested by S. Meyers.)
%   L. Hinnov, December, 2006 - Created from Taner (2000).
%
filter = [ ];
xftfor = [ ];
xftinv = [ ];
tanerbandx = [ ];
twopi=2.*pi;
% default c=10^12;
c=10^c;
npts=length(x);
%%%%%%%%%%%%% apply zero-padding; default npad=10x
npad=10; % default zero-padding length
nlen=length(x);
fnyq=1/(2*dt);
% force into even number(power of 2) plus npad x npts zero-padding
pow=log10(nlen)/log10(2);
npow=round(pow);
dpow=pow-npow;
if dpow < 0
    npow=npow+1;
end
nptspow=2^npow;
newnlen=npad*nptspow;
df=1./(newnlen*dt);
%%%%%%%%%%%%% end zero-padding procedure
npts=newnlen;
xnpts = npts;
wnyq = twopi/(2.*dt);
dw = twopi/(xnpts*dt);
%.....Set up filter
wl = twopi*fl;
wc = twopi*fc;
wh = twopi*fh;
bw = wh - wl;
amp = 1./sqrt(2.);
arg1 = 1. - (c*log(10.)) / (20.*log(amp));
arg1 = log(arg1);
arg2 = ((bw+2)/bw)^2;
arg2 = log(arg2);
twod = 2.*arg1/arg2;
ncut = fix(npts/2+1);
%.....Filter for zero and positive frequencies; f, filtout for output
f=[ ];
for n=1:ncut;
w = (n-1)*dw;
f(n)=(n-1)*df; 
arg = (2.*abs(w-wc)/bw)^twod;
darg = -1.0*arg;
filter(n) = (amp * exp(darg));
end
filtout=filter*(sqrt(2)); % normalized for plotting 
ncut1 = ncut+1;
%.....Filter for negative frequencies
for n=ncut1:npts;
w = (n-npts-1)*dw;
aw = abs(w);
arg = (2.*abs(aw-wc)/bw)^twod;
darg=-1.0*arg;
filter(n) = (amp * exp(darg));
end
%.....Forward FFT data
xftfor = fft(x,npts);
%.....Apply the filter over all frequencies
for n=1:npts;
xftfor(n) = xftfor(n) * filter(n);
end
%.....Inverse FFT; normalize by sqrt(2)
xftinv = ifft(xftfor,npts);
tanerbandx = sqrt(2)*real(xftinv);
%....Remove zero-padding
tanerbandx=tanerbandx(1:nlen);
end


