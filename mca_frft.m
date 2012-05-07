function [Z,A,U] = mca_frft(x,alpha,u,varargin)
% MCA_FRFT  computes the Fractional Fourier Transform plane for each alpha
%
% see also MCA_IFESTIMATE, MCA_EXTRACT

% optional parameters
PLOTMODE = false;
if nargin > 3
    PLOTMODE = varargin{1};
    if PLOTMODE
        tic;        % start timer
    end
end

% adjust data length to length of u by zero padding
N = numel(u);
nPad = N-numel(x);
if nPad > 0;
    x = [x(:); zeros(nPad,1)];
end

% perform fractional fourier transforms over all alpha values
y = zeros(N,numel(alpha));       % init matrix
for k=1:numel(alpha)
    y(:,k) = frft(x,alpha(k));
end

% take squared modulus over each FrFT (produces the Radon Wigner Transform)
Z = abs(y).^2;

% calculate grid coordinates
if (nargout > 1) || PLOTMODE
    [A,U] = meshgrid(alpha,u);
end




%% Plot FrFT
if PLOTMODE
    toc                 % display timer result
    
    dbrange = 60;
    
    % plot FrFT
    figure;
    pcolor(A,U,db(Z))
    %surf(A,U,db(Z))
    shading flat
    colorbar
    cLim = get(gca,'clim');
    set(gca,'clim',[cLim(2)-dbrange cLim(2)]);
    colormap(flipud(hot))

end

