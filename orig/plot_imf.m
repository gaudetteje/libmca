function plot_imf(varargin)
% PLOT_IMF  displays recorded and post-processed flight data

% plot parameters
PLOT_FLT = 0;
PLOT_WAV = 1;
PLOT_SPEC = 1;
PLOT_PWVD = 0;
PLOT_IMF = 0;
PLOT_HIL = 0;
nfft = 64;
cMin = -80;     % electronic noise floor
cMax = -40;     % loudest signal received on average

CH = 1:23; %6:9;
CALL = 55; %[];

%cwd = '../data/20100128/trial15/processed';
cwd = '../data/20100317/trial06/processed';
%cwd = '../data/20100610/trial22/processed';

% set default figure properties
set(0,'DefaultFigureNumberTitle','off')
set(0,'DefaultFigureColor','white')
set(0,'DefaultFigureMenuBar','none')  % {'figure','none'}
set(0,'DefaultFigureToolBar','figure')  % {'auto','figure','none'}
set(0,'DefaultFigureNextPlot','Replace')

% search for data files
if isempty(CALL)
    pat = '(imf).*(\.mat$)';
else
    pat = sprintf('(imf).*(%2d\\.mat$)',CALL);
end
procdata = findfiles(cwd,pat,Inf,true);
fprintf('Found a total of %d processed calls\n',length(procdata))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% iterate over each processed data file found
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for fNum = 1:length(procdata)
    
    tic
    % load the MAT data file
    load(procdata{fNum})
    fprintf('Loaded data file in %f seconds\n',toc)
    
    % initialize figure windows
    tilefigs(CH,8,3);       % dual displays
    %tilefigs(CH,6,4);       % single display
    %tilefigs(CH,2,2);       % subset of channels
    %tilefigs(CH,3,4);
    
    % length of data samples
    L = size(call,1); %round(diff(flight.audio.absolute.call_intervals(:,cNum)) * fs)
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% display results
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % plot flight path
    if PLOT_FLT
        plotheadaim(flight,'calls',cNum,'figure',100);
    end
    
    % plot original waveform
    if PLOT_WAV
        t = (1:L) .* 1e3 ./ fs;
        amp = max(max(abs(call)));          % find largest signal amplitude for normalization
        for m=CH
            if all(~call(:,m))
                figure(m)
                set(m,'Visible','off')
            else
                figure(m)
                plot(t, call(:,m))
                grid on
                set(gca,'yLim',[-amp amp])
                xlabel('Time (ms)')
                ylabel('Amplitude (normalized)')
                set(m, 'Name', sprintf('BUMP Ch %d - Call #%d',m,cNum));
            end
        end
        
        pause
    end
    
    % plot signal spectrogram
    if PLOT_SPEC
        for m=CH
            if all(~call(:,m))
                figure(m)
                set(m,'Visible','off')
            else
                figure(m)
                [y,f,t,p] = spectrogram(call(:,m),nfft,round(nfft*0.8),[],fs*1e-3,'yaxis');
                surf(t,f,10*log10(abs(p)),'EdgeColor','none');
                axis xy; axis tight; colormap(jet); view(0,90);
                xlabel('Time (ms)')
                ylabel('Frequency (kHz)')
                caxis([cMin cMax])
                colorbar
                set(m, 'Name', sprintf('BUMP Ch %d - Call #%d',m,cNum));
            end
        end
        
        pause
    end
    
    % plot Pseudo Wigner-Ville Distribution
    if PLOT_PWVD
        for m=CH
            if all(~call(:,m))
                figure(m)
                set(m,'Visible','off')
            else
                figure(m)
                
                tfrpwv(hilbert(call(:,m)));
            end
        end
        
        pause
    end

    % plot IMFs of each call
    if PLOT_IMF
        for m=CH
            if all(~call(:,m))
                figure(m)
                set(m,'Visible','off')
            else
                figure(m)
                % plot each IMF on its own subplot
                for n=1 %:size(IMF{m},2)
                    %subplot(size(IMF{m},1),1,n)
                    [y,f,t,p] = spectrogram(IMF{m}(:,n),nfft,round(nfft*0.8),[],fs*1e-3,'yaxis');
                    surf(t,f,10*log10(abs(p)),'EdgeColor','none');
                    axis xy; axis tight; colormap(jet); view(0,90);
                    ylabel('Frequency (kHz)')
                    caxis([cMin cMax])
                    colorbar
                    set(m, 'Name', sprintf('BUMP Ch %d - Call #%d',m,cNum));
                end
                xlabel('Time (ms)')
            end
        end

        pause
    end

    % plot IF/IA
    if PLOT_HIL
        figure

        pause
    end
    
end
