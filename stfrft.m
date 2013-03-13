function y = stfrft(x,varargin)
% STFrFT  computes the short-time fractional Fourier transform
%
% y = stfrft(x) uses the default values
%
% wind, noverlap, nfft, nalpha, 

DEBUG = false;

% default values
wind = @hamming;
noverlap = 0;
nfft = 100;
nalpha = 200;

% parameter checking
if nargin > 5
    error('Incorrect number of parameters entered')
end

% assign optional parameters
params = {wind, noverlap, nfft, nalpha};
params(1:numel(varargin)) = varargin;
[wind, noverlap, nfft, nalpha] = params{:};

% input data formatting
if size(x,1) == 1
    x = x(:);       % force into column vector if single channel
end

% ensure signal is analytic signal
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
w = w(:);           % force into column vector

% find starting index for each block
block = (1:nfft-noverlap:L);

% FrFT indices
alpha = linspace(-1,1,nalpha);
mu = linspace(0,1,nfft);

% init values
mu0 = zeros(numel(block),1);
f0 = zeros(numel(block),1);

% iterate over each block
for bNum = 1:numel(block)
    bIdx = (block(bNum):block(bNum)+nfft-1);
    
    % zero pad if last block too short
    if bIdx(end) > L
        x(L+1:bIdx(end)) = 0;
    end
    
    % extract and window data block
    b = x(bIdx) .* repmat(w,1,M);
    
    % compute FrFT
    [Z,A,U] = mca_frft(b,alpha,mu,DEBUG);
    
    % estimate dominant LFM
    [aMax,uMax,zMax] = find_peaks(Z,A,U,1);
    [~,mu0(bNum),f0(bNum)] = lfm_est(aMax,uMax,1,nfft);
    
    %pause
end

y = x;



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
