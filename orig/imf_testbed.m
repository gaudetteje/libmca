
clc
clear
close all

% params
MODE = 2;
CH = 16;
CALL = 50;
datapath = '/Volumes/freeflight_backup/freeflight_backups/20100611/MATLAB/flightroom_data/';
flightstruct = 'flt_struc_2-02-10flight_02.mat';


PLOT_PROCESS = 1;
PLOT_SYNTH = 0;

% set default figure properties
set(0,'DefaultAxesFontSize',14)
set(0,'DefaultFigureNumberTitle','on') %{'on','off'}
set(0,'DefaultFigureColor','white')
set(0,'DefaultFigureMenuBar','none')  % {'figure','none'}
set(0,'DefaultFigureToolBar','figure')  % {'auto','figure','none'}
set(0,'DefaultFigureNextPlot','Replace')


% temporarily add 3rd party analysis tools
addpath(genpath('3rd_party'))

switch MODE
    case 1
        if exist(fullfile(datapath,flightstruct),'file')
            load(fullfile(datapath,flightstruct),'flight')
            
            plotheadaim(flight,'figure',99,'calls',CALL,'mics',CH,'viewangle',[0 90]);
            data = flight.audio.channel(CH).call(CALL).wfcut';
            data(data==0)=[];       % prune out noise
            fs = flight.param.audparams.Fs;
            
        else
            error('File not found!!')
        end
    case 2
        if ~exist('flight','var')
            load('battest','flight')
        end
        plotheadaim(flight,'figure',99,'calls',CALL,'mics',CH,'viewangle',[0 90]);
        data = flight.audio.channel(CH).call(CALL).wfcut';
        data(data==0)=[];       % prune out noise
        fs = flight.param.audparams.Fs;

    case 3
        %load('D:\gaudetteje\src\biscat_signals\brown_data\JITSIG250kHz.mat')
        load('/Users/jasongaudette/src/biscat_signals/brown_data/JITSIG250kHz.mat')
        data = ts.data(:);
        fs = ts.fs;
        clear ts;
end

[IF,IA] = process_via_emd(data,fs,PLOT_PROCESS);
if PLOT_PROCESS
    tilefigs([1:12],3,4)
end


% synthesize analytic components and combine
phi0 = -pi/2;
FM1 = gen_ifpulse(fs, IF{1}(:,1), phi0, IA{1}(2:end-1,1));
FM2 = gen_ifpulse(fs, IF{2}(:,1), phi0, IA{2}(2:end-1,1));
Y = FM1+FM2;
t = (0:length(data)-1).'./fs;


if PLOT_SYNTH
    figure
    plot(t,data,t(1:length(Y)),real(Y),'--r')
    title('Time series signals')
    legend('Original signal','Synthesized signal')

    figure;
    plot(IA{1}(:,1)); hold on;
    plot(IA{2}(:,1),'r')
    title('Instantaneous Amplitude')
    legend('FM1','FM2')

    figure;
    plot(IF{1}(:,1)); hold on;
    plot(IF{2}(:,1),'r')
    title('Instantaneous Frequency')
    legend('FM1','FM2')

    figure
    spectrogram(data,64,60,[],fs,'yaxis')
    title('Original signal')

    figure
    spectrogram(real(FM1),64,60,[],fs,'yaxis')
    title('Synthesized FM1')

    figure
    spectrogram(real(FM2),64,60,[],fs,'yaxis')
    title('Synthesized FM2')

    figure
    spectrogram(real(Y),64,60,[],fs,'yaxis')
    title('Combined FM1 and FM2')

    figure
    plot(t(1:length(FM1)), abs(FM1)); hold on;
    plot(t(1:length(IA{1})), IA{1}(:,1),'--r');
    title('Instantaneous Amplitude - FM1')
    legend('Synthesized IA','Estimated IA')

    figure
    idx = 2:size(FM1,1)-1;
    IP  = -FM1(idx+1,:) .* conj(FM1(idx-1,:));
    plot(t(1:length(IP)), fs*(angle(IP)+pi)./(4*pi)); hold on;
    plot(t(1:length(IF{1})), IF{1}(:,1),'--r');
    title('Instantaneous Frequency - FM1')
    legend('Synthesized IF','Estimated IF')
    
    tilefigs(1:24,3,4)
end

% bring flightroom to front
tilefigs(99,2,6)
