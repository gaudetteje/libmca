function y = mca_iffilt(x,IF,varargin)
% MCA_IFFILT  applies a time-varying bandpass filter to a signal
%
% Y = mca_iffilt(X,IF) filters the signal, X, around the desired
%     instantaneous frequency, IF.  Fs is 1 by default.
% Y = mca_iffilt(X,IF,Fs) specify a sampling rate, Fs, that matches IF
% Y = mca_iffilt(X,IF,Fs,BW) sets a different bandwidth for the filter
% Y = mca_iffilt(X,IF,Fs,B,A) applies a lowpass filter with coefficients (B,A).
%
% The following set of options will plot the procedure at each step
% Y = mca_iffilt(X,IF,true)
% Y = mca_iffilt(X,IF,Fs,true)
% Y = mca_iffilt(X,IF,Fs,BW,true)
% Y = mca_iffilt(X,IF,Fs,B,A,true)
%
% Method:  Based on the desired IF, we calculate the phase law, PL, and
% warp the signal along this line.  Effectively, the signal is basebanded
% with a time dependent modulation.  An LTI low-pass filter is applied to the
% warped complex signal to remove any components outside of this band.  The
% forward-backward filter (help filtfilt) is used to retain any phase
% information.  The complex signal is then demodulated back to its original
% form and returned.
%
% Notes:
%
% A Chebychev II lowpass filter is used by default.  The 3dB cutoff is at 0.052*Fs
%
% Estimation of the phase law from instantaneous frequency may require
% upsampling the data if results are not satisfactory
%
% The input signal, X, is converted to the analytic form and must be real.
% The output, Y, will always be complex.
%
% Example - LFM pulse:
%   fs = 1e3;
%   f0 = 400;  f1 = 100;  T = 0.5;
%   t = (0:1/fs:T);
%   phi = f0*t + 0.5*(f1-f0)/T*t.^2;
%   x = 2*cos(2*pi*phi) + randn(size(t));
%   IF = f0*ones(size(t)) + (f1-f0)/T * t;
%   mca_iffilt(x,IF,fs);
%
% see also MCA_IFESTIMATE, MCA_EXTRACT

% set default parameters
PLOTFLAG = false;
fs = 1;
B = .1;     % 10% bandwidth

% force into column vectors
if size(x,1) == 1
    x = x(:);
end
if size(IF,1) == 1
    IF = IF(:);
end

% input data and IF must be column vectors
if size(x,2) > size(x,1)
    warning('X has more channels than samples, verify data is in column vectors')
end
if size(IF,2) ~= 1 && size(IF,2) ~= size(x,2)
    warning('IF must be a single column vector or have the same number of column vectors as X')
end
assert(size(x,1) == size(IF,1), 'X and IF must be equal length vectors')

% define filter
switch nargin
    case 2
    case 3
        if islogical(varargin{1})
            PLOTFLAG = varargin{1};
        else
            fs = varargin{1};
            B = B * fs;
        end
    case 4
        fs = varargin{1};
        if islogical(varargin{2})
            PLOTFLAG = varargin{2};
            B = B * fs;
        else
            B = varargin{2};
        end
    case 5
        fs = varargin{1};
        if islogical(varargin{3})
            B = varargin{2};
            PLOTFLAG = varargin{3};
        else
            b = varargin{2};
            a = varargin{3};
            N = max(numel(b),numel(a));
        end
    case 6
        fs = varargin{1};
        b = varargin{2};
        a = varargin{3};
        N = max(numel(b),numel(a));
        PLOTFLAG = varargin{4};
    otherwise
        error('Incorrect number of parameters entered')
end

% design FIR filter with Parks-McClellan equiripple algorithm
if ~exist('b','var')
    N = 10;
    B = B/.66;      % for 10th order cheby2 filter, 1dB point is at .66 x B. (3dB @ .705 X B)
    [b,a] = cheby2(N,80,B/fs);
    %from original Process_via_EMD.m:  [b,a] = butter(20,2*(xloc(LOCS(1))-20000)/fs);
end

% meet required data length for filtfilt (3*N)
nPad = 0;
if length(x) < 3*N
    nPad = 3*(N+1)-length(x);
    z = zeros(nPad,1);
    x = [x; z];
    IF = [IF; z];
    warning('Data must be 3x length of filter for filtfilt.m to work - zero padding data by %d samples.',nPad);
end
t = ((0:length(x)-1)/fs)';

% convert to real
if ~isreal(x)
    x = real(x);  %we could potentially save complex data by shifting everything up
end

% make signal analytic
xhat = hilbert(x);

% integrate over the IF to get phase, i.e. int d phi/dt = phi
PL = 2*pi*cumtrapz(t',IF')';     % estimate phase law from IF

% demodulate the signal - this will warp the freqency domain
demod = xhat.*exp(-1i*PL);

% filter the demodulated signal
resdemod = filtfilt(b,a,demod);

% remodulate the filtered result
y = resdemod.*exp(1i*PL);
%y = y(1:end-length(z));

% remove extra samples, if added originally
if nPad
    y = y(1:end-nPad);
end

% return results
if PLOTFLAG
    %fs = 250e3;         % HARDCODED UNTIL BUG FIX
    
    nfft = 256;
    winlen = 64;
    %wind = gausswin(nfft)';
    
    %figure; freqz(b,a,1024);
    %title('Filter Response')
    
    fh1 = figure;
    plotSpecgram(xhat,winlen,nfft,fs);
    cMax = get(gca,'clim');
    hold on; plot3(t/fs*1e3,IF*fs*1e-3,zeros(size(IF)),'g','linewidth',2)
    %title('Original Signal (real)'); colorbar
    
    fh2 = figure;
    plotSpecgram(demod,winlen,nfft,fs);
    %title('Demodulated Signal (complex)'); colorbar
    
    fh3 = figure;
    plotSpecgram(resdemod,winlen,nfft,fs);
    %title('Filtered Demodulated Signal (complex)'); colorbar
    
    fh4 = figure;
    plotSpecgram(y,winlen,nfft,fs);
    %title('Real part of filtered signal (real)'); colorbar
    
end


function plotSpecgram(x,winlen,nfft,fs)
    dBrange = 40;

    [~,F,T,P] = spectrogram(x,hamming(winlen),winlen-2,nfft,fs,'yaxis');
    if ~isreal(x)
        P = fftshift(P,1);      % shift imaginary part below DC
        F = [-flipud(F(1:end/2+1)); F(2:end/2)];
    end
    P = P./max(max(P));     % normalize
    imagesc(T*1e3,F*1e-3,10*log10(P));
    colormap(flipud(hot))
    set(gca,'ydir','normal');
    %set(gca,'xLim',[0 3.5]);

    cLim = get(gca,'cLim');
    cLim(1) = cLim(2) - dBrange;
    set(gca,'cLim',cLim)
    colorbar

    xlabel('Time (ms)')
    ylabel('Frequency (kHz)')
    %title(sprintf('f_i(t) = %g t + %g',mu0,f0))
