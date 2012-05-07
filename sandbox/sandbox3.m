clear
clc
close all

% calibration data file
patt = 'calmic_1.bin';
idx1 = 90:1400;
idx2 = 1490:2800;
t0 = 100;
N = 1;


files = findfiles('.',patt);

for fnum = 1:length(files)
    fname = files{fnum};
    fname = fname(1:end-4);
    
    % read data files
    res = niLoadBin([fname '.bin'], '.', 2);
    param = niLoadPrm([fname '.prm'], '.');
    
    % extract input/output into separate time series
    ts1.fs = 5e5;
    ts1.time = ts1.fs .* (1:length(idx1));
    ts1.data = res(idx1,1);
    ts1.data = ts1.data./max(ts1.data(:));                           % normalize data
    ts1.data = ts1.data - ones(length(idx1),1)*mean(ts1.data);       % subtract DC offset
    
    ts2.fs = 5e5;
    ts2.time = ts2.fs .* (1:length(idx2));
    ts2.data = res(idx2,2);
    ts2.data = ts2.data./max(ts2.data(:));                           % normalize data
    ts2.data = ts2.data - ones(length(idx2),1)*mean(ts2.data);       % subtract DC offset
    
    
    % Examine time series and spectrogram of data for a quick-look
    % plot time series
    figure(1)
    subplot(2,1,1)
    plot(ts1.time,ts1.data)
    grid on;
    title(sprintf('%s - Channel 1',fname),'interpreter','none')
    
    subplot(2,1,2)
    plot(ts2.time,ts2.data)
    grid on;
    title(sprintf('%s - Channel 2',fname),'interpreter','none')
    
    % plot spectrograms
    figure(2)
    spectrogram(ts1.data,hann(256),250,256,ts1.fs,'yaxis');
    title(sprintf('%s - Channel 1',fname),'interpreter','none')
    set(gca,'clim',[-120 -40])
    colorbar
    colormap jet
    
    figure(3)
    spectrogram(ts2.data,hann(256),250,256,ts2.fs,'yaxis');
    title(sprintf('%s - Channel 2',fname),'interpreter','none')
    set(gca,'clim',[-120 -40])
    colorbar
    colormap jet


    %% calculate system impulse response
    
    % Take Fourier transform of input and output signals
    nfft = 8192;
    X = fft(ts1.data,nfft);
    Y = fft(ts2.data,nfft);
    F = (0:(nfft-1))*(ts1.fs/nfft);
    
    % Compute the frequency response at all frequencies
    H = Y./X;
    
    % Find which fft samples are in the chirp freq range
    f1 = 10e3; % chirp low freq
    f2 = 100e3; % chirp high freq
    k1 = ceil((f1/ts1.fs)*nfft)+1;  % which fft sample corresponds to f1
    k2 = floor((f2/ts1.fs)*nfft)+1; % which fft sample corresponds to f2
    
    % Set H to unity below f1 and -20dB above f2 (lowpass filter effect)
    H(1:(k1-1)) = ones(k1-1,1);
    H((k2+1):(nfft/2)) = 0.1 .* ones((nfft/2)-k2,1);
    H((nfft/2)+1) = 0.1;
    H(((nfft/2)+2):nfft) = conj(flipud(H(2:(nfft/2))));
    
    % Smooth output signal's magnitude response, phase is effectively 0
    % Not sure why there are lots of deviations:
    %   Could be due to reverb or other reflections interfering with
    %   recorded microphone signal)
    Hdmag = abs(H);
    Hdmag = sgolayfilt(Hdmag,2,1001);       % use SG filter
    Hdmag = filtfilt(ones(50,1)./50,1,Hdmag);        % use MA filter
    [Hdr,Hdi] = pol2cart(zeros(size(Hdmag)),Hdmag);
    Hd = Hdr + 1i.*Hdi;
    
    
%     [Hphs,Hmag] = cart2pol(real(H),imag(H));
%     Hmag = sgolayfilt(Hmag,3,1001);
%     [Hr,Hi] = pol2cart(Hphs,Hmag);
%     Hd = Hr + 1i.*Hi;
    
    % calculate impulse response
%    hd = real(ifft(Hd));
    
    % flip signal to get correct time origin
%    hd = fftshift(hd);
    
    % prune down impulse response to more managable time length
    %M = 150;
    %idx = (-M/2:M/2)+4100;
    %hd = hd(idx);
    
    % plot input output response
    figure(4)
    subplot(2,1,1)
    plot(F,db(abs(X)),F,db(abs(Y)),'r',F,db(abs(H)),'k',F,db(abs(Hd)),'--k')
    grid on
    legend('Input Signal','Output Signal')
    ylabel('Magnitued [dB]')
    title('Magnitude response of Input/Output')
    
    subplot(2,1,2)
    plot(F,unwrap(angle(X))*pi/180,F,unwrap(angle(Y))*pi/180,'r',F,unwrap(angle(H))*pi/180,'k')
    grid on
    title('Phase response of Input/Output')
    ylabel('Phase [Degrees]')
    xlabel('Frequency [Hz]')
    
    % plot resulting impulse response
%     figure(5)
%     plot(hd);
%     grid on;
%     title('System Impulse Response')
%     xlabel('Time (samples)')
%     ylabel('Amplitude')
    
    

    %% estimate poles/zeros from impulse response, h(t)
    Nb = 10;                % number of zeros
    Na = 10;                % number of poles
    maxIter = 10;           % maximum iterations with Steiglitz-McBride algorithm
    
    Hdmag = abs(Hd(1:64:end));      % downsample TF and extract magnitude only
    L = length(Hdmag);
    W = dftmtx(L); Wb = W; Wa = W;
    Wb(:,Nb+2:L) = []; Wa(:,Na+2:L) = [];
    
    % generate the autocorrelation function
    r = ifft(Hdmag.^2);
    
    % construct an initial system model by Levinson-Durbin (AR), follow with Prony (ARMA)
    aL = levinson(r,floor(L/2));
    hL = impz(1,aL,Nb+2*Na+2);
    [b,a] = prony(hL,Nb,Na);
    
    % iteratively refine pole/zero positions with frequency domain Steiglitz-McBride (ARMA)
    for i = 1:maxIter,
        [Hi,w] = freqz(b,a,L,'whole');
        Hai = freqz(1,a,L,'whole');
        Pi = exp(1i*angle(Hi));
        HdPi = Hdmag.*Pi;
        b = (diag(Hai)*Wb)\HdPi; B = fft(b,L);
        a = (diag(HdPi.*Hai)*Wa)\(diag(Hai)*Wb*b);
    end
    
    % force real filter coefficients (Hd should be symmetric)
    if (sum(imag(b)) + sum(imag(a)) > 1e-10)
        warning('Poles and/or zeros not entirely real.  Possibly throwing away significant imaginary parts.')
        fprintf('Note:  The desired magnitude response, Hd, should be kept symmetric to ensure real coefficients!\n')
    end
    b = real(b.');
    a = real(a.');
    
    % scale all coefficients evenly, forcing a0=1
    b = b/a(1);
    a = a/a(1);
    
    
    figure(6)
    [Hz1,FF] = freqz(b,a,1024);
    Hz2 = freqz(a,b,1024);
    subplot(2,1,1)
    plot(FF,db(abs(Hz1)),FF,db(abs(Hz2)),'g')
    grid on
    title('Magnitude Response of Estimated and Inverted Filter')
    ylabel('Magnitude [dB]')
    
    subplot(2,1,2)
    plot(FF,unwrap(angle(Hz1))*pi/180,FF,unwrap(angle(Hz2))*pi/180,'g')
    grid on
    title('Phase Respanse of Estimated and Inverted Filter')
    ylabel('Phase [Degrees]')
    xlabel('Frequency [Hz]')
    
    figure(7)
    zplane(a,b)
    
    %% test filter by running output through filter - should match input signal
    
    ts3 = ts2;
    ts3.data = filter(a,b,ts2.data);
    figure(8)
    spectrogram(ts3.data,hann(256),250,256,ts1.fs,'yaxis');
    title(sprintf('%s - Corrected Channel 2',fname),'interpreter','none')
    set(gca,'clim',[-120 -40])
    colorbar
    colormap jet

    
end

tilefigs(2,3)