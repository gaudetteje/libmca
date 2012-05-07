function OUT=mca_iterate_estimate(DATA,IF_est0,filt_BW,Fs)

% perform IMD analysis
% create a rough estimate of first harmonic IF
% define the fundamental IF used for the HFM pulse
PLOTTEST=1;
OUT.Fs=Fs;


LEEDIN=0;
LEEDOUT=0;
if IF_est0<0
    T = 2.5e-3;
    t = (0:1/OUT.Fs:T)';
    f0 = 50e3;  f1 = 25e3;
    B = abs(f1-f0);
    a = T*(f0*f1)/B;
    b = T*f1/B;
    IF_est0 = a./(t+b);
    IF_est0 = [zeros(round(OUT.Fs*LEEDIN),1); IF_est0; zeros(round(OUT.Fs*LEEDOUT),1)];
    IF_est0 = IF_est0(1:length(DATA));
end


if PLOTTEST
    close all;
    FIG(1)=figure('color','w')
    plot(IF_est0);
end

OUT.DUR=T;


IF_est=IF_est0;
count=0;



while count<100
    count=count+1;
    %analize pulse
    %mca_iffilt(PULSE,IF_est,OUT.Fs);
    res = mca_iffilt(DATA,IF_est,OUT.Fs,filt_BW); % filter out component from signal
    %mca_extract(res,OUT.Fs,0);
    [IMF, IAest ,IF_est] = mca_extract(res,OUT.Fs,0); % extract IF and IA for each IMF
    
    if PLOTTEST
        hold off
        plot(IF_est0,'r','markersize',1);
        hold on
        plot(IF_est,'k')
        pause(.001)
       
    end
    
end

FIG(2)=figure('color','w')
subplot(2,1,1)
spectrogram(DATA,256,250,256,Fs,'yaxis')
subplot(2,1,2)
OUT.Tif=(1:length(IF_est))/OUT.Fs;
plot(OUT.Tif,IF_est,'k'); hold on;
plot(OUT.Tif,IF_est0,'.r','markersize',.5)


OUT.IF_est=IF_est;
OUT.IF_est0=IF_est0;


