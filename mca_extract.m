function [IMF,IA,IF,IP] = mca_extract(X,fs,varargin)
% MCA_EXTRACT  calculates the IA and IF for each IMF extracted from a monocomponent signal
%
% IMF = mca_extract(X,fs) calculates up to 3 intrinsic mode functions on signal
%       X sampled at rate fs.
% [IMF,IA,IF] = mca_extract(X,fs) returns the IMFs along with instantaneous
%       frequency, IF, and instantaneous amplitude, IA, for each IMF.
% [IMF,IA,IF,PL] = mca_extract(X,fs) also returns the phase law, PL, used to
%       estimate IF
% [...] = mca_extract(X,fs,nModes) calculates M-1 IMF functions along with
%       the residual signal, IMF(:,end).  if M=0, bypass the EMD process and
%       only perform Hilbert Spectral Analysis on the input signal, X.
% [...] = mca_extract(X,fs,nModes,true) plots the resulting HSA for each IMF
% [...] = mca_extract(X,fs,nModes,pOrder,pLength) adjusts the parameters of
%       the Savitsky-Golay smoothing filter.
% [...] = mca_extract(X,fs,nModes,pOrderIF,pLengthIF,pOrderIA,pLengthIA)
%       uses different filter options for IF and IA.
% [...] = mca_extract(X,fs,nModes,pOrder,pLength,true) also plots the HSAs
%
% Note:  By default, IA and IF are smoothed using a least-square optimal
%       Savitsky-Golay filter.  This can be turned off using pOrder = 0,
%       and pLength = 1.
%
% see also MCA_IFFILT, SGOLAYFILT


% set default and optional user parameters
modes = 3;                          % Number of IMFs to extract
PLOTFLAG = false;
pOrderIF = 3;                           % Savitsky-Golay polymomial order
pLengthIF = ceil(0.2*length(X));         % Savitzky-Golay filter length
pOrderIA = pOrderIF;
pLengthIA = pLengthIF;
switch nargin
    case 3
        modes = varargin{1};
    case 4
        modes = varargin{1};
        PLOTFLAG = varargin{2};
    case 5
        modes = varargin{1};
        pOrderIF = varargin{2};
        pLengthIF = varargin{3};
        pOrderIA = pOrderIF;
        pLengthIA = pLengthIF;
    case 6
        modes = varargin{1};
        pOrderIF = varargin{2};
        pLengthIF = varargin{3};
        PLOTFLAG = varargin{4};
        pOrderIA = pOrderIF;
        pLengthIA = pLengthIF;
    case 7
        modes = varargin{1};
        pOrderIF = varargin{2};
        pLengthIF = varargin{3};
        pOrderIA = varargin{4};
        pLengthIA = varargin{5};
    case 8
        modes = varargin{1};
        pOrderIF = varargin{2};
        pLengthIF = varargin{3};
        pOrderIA = varargin{4};
        pLengthIA = varargin{5};
        PLOTFLAG = varargin{6};
end

% ensure filter length less than sequence length
pLengthIF = min(pLengthIF,numel(X)-2);
pLengthIA = min(pLengthIA,numel(X)-2);

% force SG filter length to be odd (required by algorithm)
pLengthIF = pLengthIF + mod(pLengthIF+1,2);
pLengthIA = pLengthIA + mod(pLengthIA+1,2);

% extract IMFs from signal (when modes > 0)
if modes
    IMF = emd(real(X),'maxmodes',modes)';
else
    IMF = real(X);
end

% extract IF/IA from the Hilbert spectrum using the analytic approach
HT = hilbert(IMF);

% filter with Savitzky-Golay filter - LS optimal
IA = sgolayfilt(abs(HT),pOrderIA,pLengthIA);

% calculate phase law by taking the central finite difference
HT = [zeros(1,size(HT,2)); HT; zeros(1,size(HT,2))];   % pad start/end with extra sample so length(IA) == length(IF)
idx = 2:size(HT,1)-1;
IP = -HT(idx+1,:) .* conj(HT(idx-1,:));

%filter with Savitzky-Golay filter - LS optimal
IF = sgolayfilt(fs*(angle(IP)+pi)./(4*pi),pOrderIF,pLengthIF);


% plot results of first IMF
if PLOTFLAG
    t = 1e3.*(0:length(IF)-1)./fs;

    nfft = 256;
    winlen = 64;
    dBrange = 40;
    
    % plot spectrogram of each IMF
    for n=1:size(IF,2)
        figure;
        [~,F,T,P] = spectrogram(IMF(:,n),hamming(winlen),winlen-2,nfft,fs,'yaxis');
        P = P./max(max(P));     % normalize
        imagesc(T,F,10*log10(P));
        colormap(flipud(hot))
        set(gca,'ydir','normal');
        %set(gca,'ylim',[0 125])
        %set(gca,'xlim',[0 t(end)])
        
        cLim = get(gca,'cLim');
        cLim(1) = cLim(2) - dBrange;
        set(gca,'cLim',cLim)
        colorbar

        xlabel('Time (ms)')
        ylabel('Frequency (kHz)')

%        spectrogram(IMF(:,n),winlen,winlen-2,nfft,ts.fs*1e-3,'yaxis')
%        plotcol(t,1e-3.*IF(:,n).',IA(:,n).',IA(:,n).');
%        hold on
%        colormap(flipud(hot))
    end
    
    % plot lines for each IMF
    figure;
    
    subplot(2,1,1)
    plot(t, IF,'linewidth',2,'color',[0 0 .6])
    ylabel('Freq. (kHz)')
    xlabel('Time (ms)')
    %axis([0 t(end) 0 125])
    grid on
    
    subplot(2,1,2)
    plot(t, db(IA),'linewidth',2,'color',[0 0 .6])
    ylabel('Amp. (dB)')
    %axis([0 t(end) -100 -20])
    grid on
    %legend('IMF1','IMF2','IMF3','Res.','Location','SouthOutside','Orientation','Horizontal')
    
    % plot HSA for each IMF
    for n=1 %:size(IF,2)
        figure;
        c = db(IA(:,n));
        cmap = flipud(hot);
        ccplot(t,IF(:,n).',c.',cmap)
%        plotcol(t,IF(:,n).',IA(:,n).',db(IA(:,n)).');
        hold on
        colormap(cmap)
        set(gca,'clim',[min(c) max(c)])
        colorbar
        view(2)
        %set(gca,'ylim',[0 125])
        %set(gca,'xlim',[0 t(end)])
        xlabel('Time (ms)')
        ylabel('Frequency (kHz)')
    end
end
