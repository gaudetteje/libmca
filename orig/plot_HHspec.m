clear all
close all

batsignal

data                = s;

imf                 = emd(data);
imfmast             = imf';
fs                  = 142857;  %sample_rate;

range               = (0:length(data)-1)/fs;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Prune the data for noise reduction     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[rows,cols] = size(imf);
for i = 1:rows
    for j = 1:cols
        if abs(imf(i,j))<0.001
            imf(i,j) = 0; 
        else
            imf(i,j) = imf(i,j);
        end
    end
end

imfhold             = imf;
[p,k]               = size(imfhold);
hil_imf             = hilbert(imfhold');
[M,N]               = size(hil_imf);
mag                 = abs(hil_imf);
sig_conj            = conj(hil_imf);
pha_dif             = angle(sig_conj(1:(M-1),:).*hil_imf(2:M,:));
pha_dif(M,:)        = angle(hil_imf(M,:));
pha_dif             = pha_dif*(fs/(2*pi));
zeroind             = (pha_dif(:,:))<0;
pha_dif(zeroind)    = 0;

%weighting scheme
% for n = 1:M
%     w(n) = 3*M/2/(M^2-1)*(1 - ((n-(M/2-1))/(M/2))^2);
% end


[mm,nn] = size(pha_dif);

%If you only want to look at the Hilbert Spectrum for one particular IMF,
%just enter the number you want to look at as nn.

%nn = 2;
figure
for kk = 1:nn
    %if you want each Hilbert spec in its own window uncomment subplot
    %subplot(nn,1,kk)
    hold on
    %set median filtering to some "reasonable" value - depends on how much
    %data you have
    h = plotcol(range,medfilt1(pha_dif(:,kk)',21),mag(:,kk)',mag(:,kk)');

    view(2)
    set(h,'MarkerSize',10)
end

figure
for kk = 1:p
    subplot(p,1,kk),grid
    plot(range,imf(kk,:)','k')
end


figure
norm_filt = (imf(end,:));
%norm_filt = (imf(end,:)+imf(end-1,:)+imf(end-2,:));
plot(range,(imf(1,:).*(norm_filt).*((imf(1,:).*norm_filt)>0)))
figure
plot(range,(imf(1,:)./(norm_filt).*((imf(1,:)./norm_filt)>0)))






