clear
clc
close all

%% test cases to try (in this order)
%patt = 'chain27.bin';  idx = 3200:6000;  t0 = 450; N = 4;   % Synthetic echoes
%patt = 'ro_cm_16.bin';  idx = 1400:3100;  t0 = 100; N = 1;%2;  % Strong CM signal
patt = 'ro_cm_9.bin';  idx = 1400:3100;  t0 = 100; N = 1;   % Weak CM signal
%patt = '^ro_cm_.*\.bin$';      % catches all CM data


files = findfiles('.',patt);

for fnum = 1:length(files)
    fname = files{fnum};
    fname = fname(1:end-4);
    
    ts.data = niLoadBin([fname '.bin'], '.', 2);
    ts.params = niLoadPrm([fname '.prm'], '.');
    ts.fs = 5e5;
    ts.time = (1:length(ts.data)).'./ts.fs;
    %ts.data = ts.data./max(ts.data(:));        % normalize data
    ts.data = ts.data - ones(length(ts.data),1)*mean(ts.data);
    
    %% Examine time series and spectrogram of data for a quick-look
%     % plot time series
%     figure(1)
%     subplot(2,1,1)
%     plot(ts.time(idx),ts.data(idx,1))
%     grid on;
%     title(sprintf('%s - Channel 1',fname),'interpreter','none')
%     
%     subplot(2,1,2)
%     plot(ts.time(idx),ts.data(idx,2))
%     grid on;
%     title(sprintf('%s - Channel 2',fname),'interpreter','none')
%     
%     % plot spectrograms
%     figure(2)
%     spectrogram(ts.data(idx,1),hann(256),200,256,ts.fs,'yaxis');
%     title(sprintf('%s - Channel 1',fname),'interpreter','none')
%     set(gca,'clim',[-50 25])
%     colorbar
%     colormap jet
%     
%     figure(3)
%     spectrogram(ts.data(idx,2),hann(256),200,256,ts.fs,'yaxis');
%     title(sprintf('%s - Channel 2',fname),'interpreter','none')
%     set(gca,'clim',[-50 10])
%     colorbar
%     colormap jet


    %% create a rough estimate of first harmonic IF

    % define the fundamental IF used for the HFM pulse
    T = 0.003;
    t = (0:1/ts.fs:T)';
    f0 = 50e3;  f1 = 25e3;
    B = abs(f1-f0);
    a = T*(f0*f1)/B;
    b = T*f1/B;
    IF = a./(t+b);
    
    % extend offset to entire signal
    data = ts.data(idx,2);
    
    % t0 defined above - could use xcorr to determine this
    t1 = t0+length(t);
    IF_est = [f0*ones(t0,1); IF; f1*ones(length(data)-t1,1)];
    
    %% iterate over each component
    % N defined above - could use amplitude threshold to find "trackable" number of components
    IMF = cell(N,1);
    IA = cell(N,1);
    IF = cell(N,1);
    for n=1:N
    
        % filter out component from signal
        res = mca_iffilt(data,n.*IF_est,ts.fs);
        
        fh = figure;
        time = ((1:length(res))-t0)./ts.fs;
        plot(time,n.*IF_est,'--k','linewidth',2)
        hold on;
        grid on;
        
        
        %% Option 1 - perform Hilbert Spectral Analysis (HSA)
        
        % extract IF and IA for each IMF
        [IMF{n},IA{n},IF{n}] = mca_extract(res,ts.fs);
        
        % plot HSA in 3D w/ color
        figure(fh)
        IA{n} = IA{n} ./ max(max(IA{n}));          % normalize IA
        plotcol(time,IF{n}(:,1).',IA{n}(:,1).',IA{n}(:,1).');
        title(sprintf('HSA of Component %d',n))
        view(2)
        
        
        %% Option 2 - reassigned Pseudo-WVD
        
        % Calculate the reassigned PWV
        [PWV,rPWV] = tfrrpwv(res);  clear PWV;
        
        % Find maximum level of rPWV per time column and extract IA/IF
        [IAr,F] = max(rPWV,[],1);
        IFr = F/size(rPWV,1)*ts.fs/2;
        
        % calculate smoothing filter length
        SGfl = ceil(0.1*length(IAr));
        SGfl = SGfl + mod(length(IAr)+1,2);     % force odd filter length (required by SG algorithm)
        
        % filter IF/IA
        IAr = sgolayfilt(IAr,3,SGfl);
        IFr = sgolayfilt(IFr,3,SGfl);
        IAr = IAr ./ max(IAr);          % normalize IAr (post-filtering)
        
        % plot rPWV in 3D w/ color
        figure(fh)
        time = ((1:length(IF_est))-t0)./ts.fs;
        plotcol(time,IFr,IAr,IAr)
        
        
        %% Option 3 - reassigned Spectrogram
        [SPEC,rSPEC] = tfrrsp(res);  clear SPEC;
        
        % Find maximum level of rPWV per time column and extract IA/IF
        [IAr2,F2] = max(rSPEC,[],1);
        IFr2 = F2/size(rSPEC,1)*ts.fs ;
        
        % calculate smoothing filter length
        SGfl = ceil(0.1*length(IAr2));
        SGfl = SGfl + mod(length(IAr2)+1,2);     % force odd filter length (required by SG algorithm)
        
        % filter IF/IA
        IAr2 = sgolayfilt(IAr2,3,SGfl);
        IFr2 = sgolayfilt(IFr2,3,SGfl);
        IAr2 = IAr2 ./ max(IAr2);          % normalize IAr (post-filtering)
        
        % plot rPWV in 3D w/ color
        figure(fh)
        time = ((1:length(IF_est))-t0)./ts.fs;
        plotcol(time,IFr2,IAr2,IAr2)
        
        
    end
    
end

