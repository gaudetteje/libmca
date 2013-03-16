function [t,y,z,t0,y0,z0] = stfrft(x,varargin)
% STFrFT  computes the short-time fractional Fourier transform
%
% y = stfrft(x) uses the default values
%
% wind, noverlap, nfft, nalpha, 

% TODO:
%
% findpeaks subfunction needs to search for "local" maximums on surface
% around some area - there must be some existing algorithms to do this!
%
% Once we find a peak, search for its harmonic or subharmonics by
% constraining window

DEBUG = true;

% frft search parameters
nfft = 2^8;
nalpha = 1e3;
gamma = .001;                   % normalized threshold for valid points in TF plane

% window parameters
wind = @hamming;
noverlap = nfft * 0;

% define signal assumptions
nComp = 3;                      % maximum number of components to retain

% curve fitting parameters
FITMODE = 'poly'; %'spline';
polyOrd = 1;                    % polynomial fit order

% assign optional parameters
if nargin > 5, error('Incorrect number of parameters entered'); end
params = {wind, noverlap, nfft, nalpha};
params(1:numel(varargin)) = varargin;
[wind, noverlap, nfft, nalpha] = params{:};

% input data formatting
if size(x,1) == 1
    x = x(:);       % force into column vector if single channel
end
t = 1:numel(x);     % assign index if no sampling data entered

% ensure signal is analytic
if isreal(x)
    x = hilbert(x);
end

% split time series into multiple overlapping blocks
L = size(x,1);      % total number of samples
M = size(x,2);      % total number of channels

% define window
if ischar(wind)
    eval(sprintf('w = %s(nfft)',wind));
elseif isa(wind,'function_handle')
    w = wind(nfft);
else
    error('Unknown window entered - must be string or function handle')
end
w = w(:);           % force window into column vector

% find starting index for each block
block = floor(1:nfft-noverlap:L);

% FrFT indices
alpha = linspace(0,1.99,nalpha);
mu = linspace(0,1,nfft);

% init values
mu0 = zeros(numel(block),1);
y0 = zeros(numel(block),nComp);
z0 = zeros(numel(block),nComp);


% iterate over each block
for bNum = 1:numel(block)
    bIdx = (block(bNum):block(bNum)+nfft-1);
    
    % zero pad if last block too short
    if bIdx(end) > L
        x(L+1:bIdx(end)) = 0;
    end
    
    % extract and window data block
    b = x(bIdx) .* repmat(w,1,M);
    
    for cNum = 1:nComp
        % compute FrFT across entire rotation fraction plane
        Z = mca_frft(b,alpha,mu,DEBUG);
% DO LOW RES SEARCH HERE
        % find dominant peak and map corresponding LFM
        [aMax0,uMax0,z0(bNum,cNum)] = find_peaks(Z,alpha,mu);
        [~,mu0(bNum),y0(bNum,cNum)] = lfm_est(aMax0,uMax0,1,nfft);
        
        % find harmonic and subharmonic locations in frac-rot plane
        harm = [1./(nComp:-1:2) 1:nComp];      % possible harmonic factors
        [a0,u0] = frft_est(mu0(bNum)*harm,y0(bNum,cNum)*harm,1,nfft);
        
        % select strongest peaks from each possible location
        for h = 1:numel(harm)
            % select alpha search region
            aLower = find(alpha >= a0(h)-aRange/2, 1, 'first');
            aUpper = find(alpha <= a0(h)+aRange/2, 1, 'last');
            aIdx = aLower:aUpper;
            a1 = alpha(aIdx);
            
            uLower = find(mu >= u0(h)-uRange/2, 1, 'first');
            uUpper = find(mu <= u0(h)+uRange/2, 1, 'last');
            uIdx = uLower:uUpper;
            u1 = mu(uIdx);
            
            Z1 = Z1(uIdx,aIdx);
            [aMax1(h),uMax1(h),zMax1(h)] = find_peaks(Z1,a1,u1,1);
        end
        %[~,idx1] = sort(zMax);
        %idx1 = idx1(1:nComp);
        
        % save points to matrix
        z0(bNum,:) = zMax1(idx1);
        %aMax1 = aMax1(idx1);
        %uMax1 = uMax1(idx1);
        
        [~,mu0(bNum,:),y0(bNum,:)] = lfm_est(aMax1,uMax1,1,nfft);   % convert to LFM
    end
    
end


%% split points into multiple harmonics
%%% assume that points are always monotonically decreasing!


% fit polynomial curve to points
t0 = (block + nfft/2)';
switch (FITMODE)
    case 'poly'
        ypp = polyfit(t0,y0,polyOrd);
        y = polyval(ypp,t);

        zpp = polyfit(t0,z0,polyOrd);
        z = polyval(zpp,t);
    case 'spline'
        y = spline(t0,y0,t);
        z = spline(t0,z0,t);
    otherwise
        error('Unknown fit mode')
end

%% plot results
if DEBUG
    plot3(t,y,z);
    hold on;
    grid on;
    plot3(t0,y0,z0,'o')
    xlabel('Time')
    ylabel('Frequency')
    zlabel('Magnitude')
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

% extract N highest values
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
a0 = -2/pi .* atan(Fs^2/T ./ mu0);
BW = mu0.*T;
fc = f0 + BW./2;
u0 = (fc./Fs) .* sin(a0.*pi/2) + 0.5;
