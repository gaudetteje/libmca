function [IF,p] = mca_ifestimate(x,varargin)
% MCA_IFESTIMATE  iteratively estimates the IF of a modulated signal
%
% IF = mca_ifestimate(X) computes the FrFT over a default set of parameters
%      and returns the estimated IF of the dominant component.
% IF = mca_ifestimate(X,Fs,U,A) overrides the default parameters.
% IF = mca_ifestimate(X,Fs,U,A,THRESH) sets the threshold
% IF = mca_ifestimate(X,Fs,U,A,THRESH,true) will turn plotting on/off
%
% Input parameters:
%   FS = 1       sampling rate of analytic signal X
%   UMAX = 0.2   window range in alpha-domain to search for next peak along ridge
%   AMAX = 0.2   window range in u-domain to search for next peak along ridge
%   THRESH = 0.45  normalized threshold to locate peaks relative to global maximum
%
%   NFFT0 = 1024;
%   ADEL0 = 0.05;
%   NFFT1 = 8192;
%   ADEL1 = 0.01;
%
% Notes:
%   The signal X should first be made analytic by the Hilbert transform
%   since the FrFT will be symmetric about U and cause ambiguity.
%
% see also MCA_FRFT, MCA_IFFILT, MCA_EXTRACT

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% default user parameters

% first search stage
nfft0 = 1024;                   % total number of FrFT points in u domain - 1st cut
aDel0 = 0.05;                   % resolution of alpha domain - 1st cut

% 2nd search stage
nfft1 = 4096;                   % total number of FrFT points in u domain - 2nd cut
aDel1 = 0.01;                   % resolution of alpha domain - 2nd cut

% default parameters
aRange = 0.2;                   % window range in alpha-domain to search for ridge
uRange = 0.2;                   % window range in u-domain to search for ridge
thresh = 0.45;                  % normalized threshold for finding points along ridge 
aBound = [0.1 .975];                % minimum and maximum initial search bounds on alpha


% IF polynomial fitting parameters
FITMODE = 'poly'; %'spline'; %
polyOrd = 3;

% parse optional parameters
PLOTFLAG = 0;
GENAVI = 0;
switch nargin
    case 6
        Fs = varargin{1};
        uRange = varargin{2};
        aRange = varargin{3};
        thresh = varargin{4};
        PLOTFLAG = varargin{5};
    case 5
        Fs = varargin{1};
        uRange = varargin{2};
        aRange = varargin{3};
        thresh = varargin{4};
    case 4
        Fs = varargin{1};
        uRange = varargin{2};
        aRange = varargin{3};
    case 2
        Fs = varargin{1};
    case 1
        Fs = 1;
    otherwise
        error('Bad number of input parameters')
end


% extend frft points in u-domain if necessary
N = numel(x);       % number of samples in data
nfft0 = max(nfft0,N);
nfft1 = max(nfft1,N);

% force x to be column vector
x = x(:);
t = (0:N-1)./Fs;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Use coarse FrFT to find approximate (alpha,u) peak
alpha0 = (aBound(1):aDel0:aBound(2));
u0 = linspace(0,1,nfft0);
Z0 = mca_frft(x,alpha0,u0,PLOTFLAG);
[aMax0,uMax0,zMax0] = find_peaks(Z0,alpha0,u0);
if PLOTFLAG
    fh0 = figure(gcf); title(sprintf('(%g, %g, %g)', aMax0, uMax0, zMax0));
    hold on; plot(aMax0,uMax0,'g+','markersize',5);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Recompute localized FrFT for precise (alpha,u) peak estimate

% THIS 2ND STAGE SHOULD BE ITERATIVE IN BLOCKS, RECOMPUTING NEXT ALHPA BLOCK IF THRESHOLD
% NOT MET IN EITHER DIRECTION

% remap (aMax0,uMax0) peak to (aMax1,uMax1) if nfft1 != nfft0
[~,mu0,f0] = lfm_est(aMax0,uMax0,Fs,nfft0);
[aMax1,uMax1] = frft_est(mu0,f0,Fs,nfft1);

% constrain search for alpha within bounds
alpha1 = (max(aMax1-aRange/2, aBound(1)) : aDel1 : min(aMax1+aRange/2, aBound(2)));
u1 = linspace(0,1,nfft1);
Z1 = mca_frft(x,alpha1,u1,PLOTFLAG);

% constrain search for u within bounds
uLower = find(u1 >= uMax1-uRange/2, 1, 'first');
uUpper = find(u1 <= uMax1+uRange/2, 1, 'last');
uIdx = uLower:uUpper;
u1 = u1(uIdx);
Z1 = Z1(uIdx,:);
[aMax1,uMax1,zMax1] = find_peaks(Z1,alpha1,u1);

if PLOTFLAG
    fh1 = figure(gcf); title(sprintf('(%g, %g, %g)', aMax1, uMax1, zMax1));
    set(gca,'yLim',[u1(1) u1(end)]);
    hold on; plot(aMax1,uMax1,'g+','markersize',5);
    
    % show localized bounding box
    if nfft0 == nfft1
        figure(fh0);
        plot([alpha1(1) alpha1(end) alpha1(end) alpha1(1) alpha1(1)], ...
            [u1(1) u1(1) u1(end) u1(end) u1(1)],'c','linewidth',2)
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Iteratively trace along ridge until below threshold

% TBD - need better search algorithm for finding ridge (look for peak
% within uRange of *previous* peak and avoid using large windows

% now find the ridge within the specified threshold
Z1 = Z1./max(max(Z1));     % normalize Z
[aPeaks,uPeaks,zPeaks] = find_ridge(Z1,alpha1,u1,thresh);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% compute LFMs for each point along ridge
M = numel(aPeaks);      % number of lines found
[lfm, mu0, f0] = lfm_est(aPeaks,uPeaks,Fs,nfft1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% find intersection points of adjacent lines

% TBD - interpolate additional points between boundaries and endpoints
xInt(M+1) = N;
yInt(M+1) = mu0(end) * N + f0(end);
xInt(1) = 1;
yInt(1) = f0(1);
for k = 2:M
    xInt(k) = (f0(k)-f0(k-1)) / (mu0(k-1)-mu0(k));
    yInt(k) = mu0(k) * xInt(k) + f0(k);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% interpolate points and return final IF estimate
[xInt,sortIdx] = sort(xInt);
yInt = yInt(sortIdx);
switch (FITMODE)
    case 'poly'
        p = polyfit(xInt,yInt,polyOrd);
        IF = polyval(p,t);
    case 'spline'
        IF = spline(xInt,yInt,t);
    otherwise
        error('Unknown IF fitting mode')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot results
if PLOTFLAG

    Fs = 250e3;     % HARD CODED FOR QUICK PLOTS - NEED TO USE ACTUAL FS IN SCRIPT
    winlen = 64;
    nfft = 256;
    dBrange = 40;
    
    %% plot spectrogram of signal
    fh2 = figure;
    %subplot(2,1,1)
    [S,F,T,P] = spectrogram(real(x),hamming(winlen),winlen-2,nfft,Fs,'yaxis');
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
    
    
    % set visible range
%    zLim = get(gca,'zlim');
%    set(gca,'zlim',[zLim(1) 1]);
    hold on
    
    drawnow
    
    
    % initialize video
    if GENAVI
        aviobj1 = avifile(fullfile('vids','frft_ridge.avi'));
        aviobj1.Quality = 100;
        aviobj1.Fps = 2;
        F1 = getframe(fh1);
        aviobj1 = addframe(aviobj1,F1);

        aviobj2 = avifile(fullfile('vids','frft_lines.avi'));
        aviobj2.Quality = 100;
        aviobj2.Fps = 2;
        F2 = getframe(fh2);
        aviobj2 = addframe(aviobj2,F2);
    end
    
    
    % iteratively show progression over FrFT & STFT
    figure(fh1)
    plot(aPeaks(1),uPeaks(1),'g+','markersize',6);
    %figure(fh2)
    %lh0 = plot(1e3*t./Fs, lfm(1:N,1).*Fs*1e-3, 'c', 'linewidth', 1);
    %axis([0 1e3*N/Fs 0 (1e-3*Fs)/2])
    %lh0 = plot((1:N)./Fs, Fs*lfm(:,1), 'k');
    
    
    if GENAVI
        F1 = getframe(fh1);
        aviobj1 = addframe(aviobj1,F1);

        F2 = getframe(fh2);
        aviobj2 = addframe(aviobj2,F2);
    end
    
    %tilefigs([fh1 fh2], 2,2)
    for k=2:M
        % mark next point on FrFT
        figure(fh1)
        plot(aPeaks(k),uPeaks(k),'g+','linewidth',2,'markersize',2);
        
%         % draw next line & its intersection point on STFT
%         figure(fh2)
%         lh1 = plot(1e3*(1:N)./Fs, lfm(1:N,k).*Fs*1e-3, 'c', 'linewidth', 1);
%         %plot(1e3*xInt(k)./Fs,1e-3*Fs*yInt(k),'w+','linewidth',2,'markersize',2)
%         
%         % capture next frames
%         drawnow
%         if GENAVI
%             F1 = getframe(fh1);
%             aviobj1 = addframe(aviobj1,F1);
% 
%             F2 = getframe(fh2);
%             aviobj2 = addframe(aviobj2,F2);
%         end
%         
%         % remove lines after delay
%         %pause(0.1);
%         delete(lh0)
%         lh0 = lh1;
        
    end
%     if exist('lh1','var')
%         delete(lh1)
%     end
    
    % add delay to 2nd last frame
    drawnow
    if GENAVI
        for cnt=1:3
            F1 = getframe(fh1);
            aviobj1 = addframe(aviobj1,F1);

            F2 = getframe(fh2);
            aviobj2 = addframe(aviobj2,F2);
        end
    end
    
    % show resulting IF over spectrogram
    figure(fh2)
    plot(1e3*t./Fs, IF*Fs*1e-3, 'g', 'linewidth',2);
   
    % get last frame and save avi
    drawnow
    if GENAVI
        F1 = getframe(fh1);
        aviobj1 = addframe(aviobj1,F1);

        F2 = getframe(fh2);
        aviobj2 = addframe(aviobj2,F2);
        
        aviobj1 = close(aviobj1);
        aviobj2 = close(aviobj2);
    end
    
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%     Subfunctions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x,y,z] = find_peaks(Z,x,y,varargin)
% FIND_PEAKS  find peak locations in 2D surface
%
% [xMax,yMax,zMax] = FIND_PEAKS(Z,x,y) returns the highest peak coordinates
%     and value given the surface, Z, and index vectors, x and y.
% [x,y,z] = FIND_PEAKS(Z,x,y,N) returns column vectors of the N highest peaks

N = 1;
if nargin > 3
    N = varargin{1};
end

% sort surface as column and save index
Zcol = Z(:);
[val,idx] = sort(Zcol,'descend');

% extract N highest peaks
[ypk,xpk] = ind2sub(size(Z),idx(1:N));
x = x(xpk);
y = y(ypk);
z = val(1:N);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [lfm,mu0,f0] = lfm_est(a0,u0,Fs,T)
% LFM_EST  estimate LFM parameters from FrFT rotation angle and peak location
mu0 = -(Fs^2 / T) .* cot(a0.*pi/2);
%fc = N*(u0 - 0.5) * (Fs/N) * csc(a0*pi/2); % for u = (1:N)
fc = (u0 - 0.5) .* Fs .* csc(a0.*pi/2);        % for u = (0:1)

% convert parameters
BW = mu0.*T;
f0 = fc - BW./2;

% construct IF curve
t = (0:T-1)';
lfm = ones(size(t))*f0' + t*mu0';
lfm = Fs.*lfm;    % correct for sampling rate


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [a0,u0] = frft_est(mu0,f0,Fs,T)
% FRFT_EST  estimate peak (alpha,u) position of LFM with specified parameters
a0 = -2/pi * atan(Fs^2/T / mu0);
BW = mu0.*T;
fc = f0 + BW./2;
u0 = (fc/Fs) .* sin(a0*pi/2) + 0.5;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% deprecated subfunctions

function [xpk,ypk,zpk] = find_ridge(Z,x,y,varargin)
% FIND_RIDGE  find points along the ridge of a 2D surface starting at (x0,y0)
%
% [x,y,z] = FIND_RIDGE(Z,x,y,x0,y0,GAMMA) returns column vectors of the peaks
%     within a normalized linear threshold GAMMA

% set threshold
gamma = .75;        % default threshold
if nargin > 3
    gamma = varargin{1};
end

% set y-bound for each iteration
yDel = round(0.1 * numel(y));  % hardcode to 10% of range

N = numel(x);
xpk = NaN(N,1);
ypk = NaN(N,1);
zpk = NaN(N,1);

% find local maximum
[yLoc0,xLoc0,zMax0] = find(Z == max(max(Z)));

% % search in positive direction
% for i = 1:N
%     [val,idx] = sort(Z(:,i),'descend');
%     if val(1) >= gamma*zMax
%         xpk(i) = x(i);
%         ypk(i) = y(idx(1));
%         zpk(i) = val(1);
%     end
% end

%%%% attempt to prevent jumping to next ridge within block

% search in positive direction
yLoc = yLoc0;
for i = (xLoc0:N)
    yRng = (yLoc-yDel:yLoc+yDel);
    [zVal,yIdx] = sort(Z(yRng,i),'descend');
    if zVal(1) >= gamma*zMax0
        yLoc = yRng(yIdx(1));
        xpk(i) = x(i);
        ypk(i) = y(yLoc);
        zpk(i) = zVal(1);
    else
        break
    end
end

% search in negative direction
yLoc = yLoc0;
for i = (xLoc0-1:-1:1)
    yRng = (yLoc-yDel:yLoc+yDel);
    [zVal,yIdx] = sort(Z(yRng,i),'descend');
    if zVal(1) >= gamma*zMax0
        yLoc = yRng(yIdx(1));
        xpk(i) = x(i);
        ypk(i) = y(yLoc);
        zpk(i) = zVal(1);
    else
        break
    end
end

% remove null values
xpk(isnan(xpk)) = [];
ypk(isnan(ypk)) = [];
zpk(isnan(zpk)) = [];
