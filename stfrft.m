function [t,y,z,t0,f0,z0] = stfrft(x,varargin)
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

% plotting constants
DEBUG = true;
dBrange = 80;

% frft search parameters
nfft = 2^8;                     % size of FrFT kernel
nalpha = 1e3;                   % number of rotations to compute in 1st pass
gamma = .001;                   % normalized threshold for valid points in TF plane
aRange = 0.05;                  % harmonic peak search range in rotation
uRange = 0.025;                  % harmonic peak search range in time/frequency fraction

% window parameters
wind = @hamming;
noverlap = nfft * 0.8;

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
    x = x(:);                   % force into column vector if single channel
end
t = 1:numel(x);                 % assign index if no sampling data entered

% ensure signal is analytic
if isreal(x)
    x = hilbert(x);
end

% split time series into multiple overlapping blocks
L = size(x,1);                  % total number of samples
M = size(x,2);                  % total number of channels

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

% init cell array to hold estimates for each block
[mu0(1:numel(block))] = deal({0});
[f0(1:numel(block))] = deal({0});
[z0(1:numel(block))] = deal({0});

% iterate over each block
for bNum = 1:numel(block)
    bIdx = (block(bNum):block(bNum)+nfft-1);
    
    % zero pad if last block too short
    if bIdx(end) > L
        x(L+1:bIdx(end)) = 0;
    end
    
    % extract and window data block
    b = x(bIdx) .* repmat(w,1,M);
    
    % compute FrFT across entire rotation fraction plane
    [Z,A,U] = mca_frft(b,alpha,mu);     % NEED TO DO LOW RES SEARCH FIRST
    if DEBUG
        fh = figure;
        hold on
        pcolor(A,U,db(Z))
        shading flat
        colorbar
        cLim = get(gca,'clim');
        set(gca,'clim',[cLim(2)-dBrange cLim(2)]);
        colormap(flipud(hot));                        
    end

    % find dominant peak and map corresponding LFM
    [aMax0,uMax0] = find_peaks(Z,alpha,mu);
    [~,muMax0,fMax0] = lfm_est(aMax0,uMax0,1,nfft);
    
    % find harmonic and subharmonic locations in frac-rot plane
    harm = zeros(1,nComp^2);
    for n = 1:nComp
        harm((1:nComp)+(n-1)*nComp) = (1:nComp)./n;                    % possible harmonic factors
    end
    harm = unique(harm);
    [a0,u0] = frft_est(muMax0*harm,fMax0*harm,1,nfft);
    a0 = mod(a0,2);                             % wrap alpha points if necessary
    
    % remove harmonic locations if out of bounds
    idx = u0 < min(mu) | u0 > max(mu);
    harm(idx) = [];
    a0(idx) = [];
    u0(idx) = [];

    % select peaks from each possible location
    aMax1 = zeros(1,numel(harm));               % initialize peak vectors
    uMax1 = zeros(1,numel(harm));
    zMax1 = zeros(1,numel(harm));
    for h = 1:numel(harm)

        % select alpha search region
        aLower = find(alpha >= a0(h)-aRange/2, 1, 'first');
        aUpper = find(alpha <= a0(h)+aRange/2, 1, 'last');
        aIdx = aLower:aUpper;
        a1 = alpha(aIdx);

        % select fractional search region
        uLower = find(mu >= u0(h)-uRange/2, 1, 'first');
        uUpper = find(mu <= u0(h)+uRange/2, 1, 'last');
        uIdx = uLower:uUpper;
        u1 = mu(uIdx);

        % reduce surface to search region and find local peaks
        Z1 = Z(uIdx,aIdx);
        [aMax1(h),uMax1(h),zMax1(h)] = find_peaks(Z1,a1,u1,1);
% NEED TO FILTER OUT DETECTIONS AT EDGE TO MINIMIZE INTERFERENCE
        if DEBUG
            figure(fh)
            plot([a1(1) a1(1) a1(end) a1(end) a1(1)],[u1(1) u1(end) u1(end) u1(1) u1(1)],'k')      % enclose search region with box
            plot(aMax1(h), uMax1(h),'.g')
            drawnow
        end
    end

    % filter for uniqueness and sort by energy
    [~,idx] = unique([zMax1; aMax1; uMax1]','rows');
    idx = flipud(idx);
    
    % convert points to LFM estimate
    [~,muMax1,fMax1] = lfm_est(aMax1,uMax1,1,nfft);

    % save points to cell array
    z0(bNum) = {zMax1(idx)};
    mu0(bNum) = {muMax1(idx)};
    f0(bNum) = {fMax1(idx)};
end


%% assign time vector for each point
t0 = (block + nfft/2)';

% group points into harmonic groups

% grouping algorithm:
% - start at stongest point
% - find nearest neighbor in next point (in f0 and z0 dimension)
% - repeat for each harmonic up to nComp

%if DEBUG
    fh = figure;
    hold on
    grid on
%end

% initialize components
[f_i{1:nComp}] = deal(zeros(1,numel(t0)));
[a_i{1:nComp}] = deal(zeros(1,numel(t0)));

for bNum = 1:numel(block)-1
    f = f0{bNum};
    z = db(z0{bNum});
    
    if bNum == 1
        for n = 1:nComp
            f_i{n} = f(n);
            a_i{n} = z(n);
            continue
        end
    end
    f_next = f0{bNum+1};
    z_next = db(z0{bNum+1});
    
    plot3(t0(bNum)*ones(1,numel(f)),f,db(z),'.')
    
    % search for nearest point in frequency and amplitude
    for p = 1:nComp
        D = dist([f(p); z(p)],[f_next;z_next]);
        idx = find(min(D)); % find nearest neighbor

        % append point to component
        f_i{nComp} = f(idx);
        a_i{nComp} = z(idx);

        % remove item from list
        f_next(idx) = [];
        z_next(idx) = [];
        
    end
    
    
%     dt = (z0{bNum}.^2 ./ (mu0{bNum} + 1));
%     df = mu0{bNum} .* dt;
%     quiver(t0(bNum)*ones(1,numel(f0{bNum})),f0{bNum},dt,df)
    
end


% fit polynomial curve to points
switch (FITMODE)
    case 'poly'
        ypp = polyfit(t0,f0,polyOrd);
        y = polyval(ypp,t);

        zpp = polyfit(t0,z0,polyOrd);
        z = polyval(zpp,t);
    case 'spline'
        y = spline(t0,f0,t);
        z = spline(t0,z0,t);
    otherwise
        error('Unknown fit mode')
end

%% plot results
if DEBUG
    plot3(t,y,z);
    hold on;
    grid on;
    plot3(t0,f0,z0,'o')
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

a0 = a0(:)';        % force into row vector
u0 = u0(:)';

mu0 = -(Fs^2 / T) .* cot(a0.*pi/2);
%fc = N*(u0 - 0.5) * (Fs/N) * csc(a0*pi/2); % for u = (1:N)
fc = (u0 - 0.5) .* Fs .* csc(a0.*pi/2);        % for u = (0:1)

% convert parameters
BW = mu0.*T;
f0 = fc - BW./2;

% construct IF curve
t = (0:T-1)';
lfm = ones(size(t))*f0 + t*mu0;
lfm = Fs.*lfm;    % correct for sampling rate


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [a0,u0] = frft_est(mu0,f0,Fs,T)
% FRFT_EST  estimate peak (alpha,u) position of LFM with specified parameters
a0 = -2/pi .* atan(Fs^2/T ./ mu0);
BW = mu0.*T;
fc = f0 + BW./2;
u0 = (fc./Fs) .* sin(a0.*pi/2) + 0.5;
