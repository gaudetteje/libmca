% close all
% clear
% clc
% 
% load X;
% x = X-mean(X);
% x = hilbert(x);
% x = x(301:1000);
% 
% x = [zeros(500,1); x; zeros(500,1)];
% 
% IF1 = mca_ifestimate(x);
% x1 = mca_iffilt(x,IF1,1);
% 
% IF2 = mca_ifestimate(x1-x);
% x2 = mca_iffilt(x-x1,IF2,1);
% 
% tilefigs(2,3)

%%
close all
clear
clc

% synthesize 1st LFM
N = 200;
t = (0:N-1)';
f0 = 0.05;
f1 = 0.1;
mu0 = (f1-f0)/N;
phiref = f0*t + mu0/2 *t.^2;
x1 = .05*exp(1i*2*pi*phiref);

% synthesize 2nd LFM
f0 = .1;
f1 = .2;
mu0 = (f1-f0)/N;
phiref = f0*t + mu0/2 *t.^2;
x2 = .06*exp(1i*2*pi*phiref);

% zero pad
x = [zeros(100,1); x1; zeros(100,1)];

% process
IF1 = mca_ifestimate(x);
x1 = mca_iffilt(x,IF1(:),1);
%figure; tfrwv(x1);
figure; spectrogram(real(x1),64,62,256,1,'yaxis')

IF2 = mca_ifestimate(x-real(x1));
x2 = mca_iffilt(x-real(x1),IF2,1);
%figure; tfrwv(x2);
figure; spectrogram(x2,64,62,256,1,'yaxis')

tilefigs(2,4)


% %%
% 
% close all
% clear
% clc
% 
% % synthesize 1st HFM
% T = 500;
% t = (0:T-1);
% 
% f0 = .2;
% f1 = .05;
% B = abs(f1-f0);
% 
% a = T*(f0*f1)/B;
% b = T*f1/B;
% phiref = a*log(t+b);
% 
% x1 = .05*exp(1i*2*pi*phiref);
% 
% f0 = .4;
% f1 = .1;
% B = abs(f1-f0);
% 
% a = T*(f0*f1)/B;
% b = T*f1/B;
% phiref = a*log(t+b);
% 
% x2 = .05*exp(1i*2*pi*phiref);
% 
% 
% % zero pad
% x = [zeros(1,200) x1+x2 zeros(1,200)];
% 
% % process
% IF = mca_ifestimate(x);
% tilefigs(2,2)
% 
