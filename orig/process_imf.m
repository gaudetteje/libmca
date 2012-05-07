function process_imf(varargin)
% PROCESS_IMF  decomposes multi-channel acoustic data into intrinsic mode
% functions and generates the Hilbert spectrum (IF/IA) for each mode.
%








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% define local parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tmppath = 'D:\gaudetteje\MATLAB\EMD\3rd_party';
%tmppath = '/Users/jasongaudette/Documents/MATLAB/EMD/src/3rd_party';

% temporarily add 3rd party analysis tools
addpath(genpath(tmppath))

% locate raw and processed data files
%cwd = '../data';
cwd = '../data/20100128/trial15';
%cwd = '../data/20100317/trial06';
%cwd = '../data/20100610/trial22';
procdata = findfiles(cwd, '(flt_struc).*(\.mat$)', Inf, true);
%rawdata = findfiles(,'(audio).*(\.wav$)',Inf,true);        % the flight structure should point to this automatically
fprintf('\nFound a total of %d flight structures\n',length(procdata))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% iterate over each flight structure found
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for fnum = 1:length(procdata)
    fprintf('\nProcessing file %d of %d: "%s"\n',fnum,length(procdata),procdata{fnum})
    
    % load the MAT data file containing a "flight" struct, bail out if necessary
    tic
    load(procdata{fnum},'flight')
    fprintf('Loaded data file in %f seconds\n',toc)
    if ~exist('flight','var')
        warning('IMF:fileNotFound','flight structure not found in file "%s"',procdata{fnum})
        continue
    end
    
    % read flight struct paramaters
    M = flight.param.audparams.nchannels;                       % number of microphones used
    nCalls = flight.audio.absolute.flightstats.ncalls;          % number of calls in flight struct
    fs = flight.param.audparams.Fs;                              % data sampling rate
    pname = fileparts(procdata{fnum});                          % working data path for current file
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% iterate over each call in struct
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    wh = waitbar(0);
    
    tic
    for cNum = 1:nCalls
        waitbar(cNum/nCalls,wh,sprintf('Analyzing call %d of %d calls',cNum,nCalls));
        fprintf('[Call #%d] Processing ',cNum)
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% extract time series waveforms for each call
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % get index from channel map
        cIdx = flight.audio.absolute.call_channel_map(:,cNum);
        
        % length of data samples
        L = 3841; %round(diff(flight.audio.absolute.call_intervals(:,cNum)) * fs)
        
        % preallocate memory
        call = zeros(L,M);
        
        % load raw data directly from flight structure
        if isfield(flight.audio.channel(1).call,'wf')
            for m = 1:M
                if isnan(cIdx(m)), continue, end                                                                    % ignore missing data
                call(:,m) = flight.audio.channel(m).call(cIdx(m)).wf;                                               % read data from struct
            end
        
        % load raw wav data within valid time intervals
        else
            for m = 1:M
                if isnan(cIdx(m)), continue, end                                                                    % ignore missing data
                sIdx = round(flight.audio.channel(m).call(cIdx(m)).call_stats.calltime * fs);                       % get sample index
                wavfile = fullfile(pname, sprintf('Audio-%d-%s.wav', m, flight.param.dataparams.take));
                if ~exist(wavfile,'file'), warning('IMF:fileNotFound','File not found: "%s"',wavfile), continue, end
                call(:,m) = wavread(wavfile, [sIdx sIdx+L-1]);                                                      % read data from file
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% iterate over each microphone signal
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        IMF = cell(1,M);
        HT = IMF;
        IA = IMF;
        IP = IMF;
        IF = IMF;
        for m = 1:M
            if all(~call(:,m)), continue, end
            
            [IF{m}, IA{m}, IMF{m}, COMPS{m}] = Process_via_EMD(call(:,m), fs);
            
%             % isolate IMFs for each call
%             IMF{m} = emd(call(:,m),'maxmodes',3)';
%             
%             % extract IF/IA from the Hilbert spectrum using the analytic approach
%             HT{m} = hilbert(IMF{m});
%             IA{m} = abs(HT{m});
%             idx = 2:size(HT{m},1)-1;
%             IP{m} = -HT{m}(idx+1,:) .* conj(HT{m}(idx-1,:));
%             IF{m} = fs*(angle(IP{m})+pi)./(4*pi);
%             %IP{m} = atan(imag(HT)./real(HT));
%             %IF{m} = fs*diff(IP{m},1)./(2*pi);
            fprintf('.')
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% save data to file
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fname = fullfile(pname, 'processed', sprintf('imf_%s_trial-%s_call-%.2d.mat', flight.param.dataparams.date, flight.param.dataparams.take, cNum));
        if ~exist(fullfile(pname,'processed'),'dir')
            if ~mkdir(fullfile(pname,'processed'))
                warning('IMF:saveData','Could not create directory in folder "%s"',pname)
            end
        end
        save(fname, 'M', 'fs', 'cNum', 'call', 'IF', 'IA', 'IMF', 'COMPS')        % need to resize variables; going from 16-bit fixed to 64-bit float
        fprintf('! (%d channels)\n',sum(~isnan(cIdx)))
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% cleanup environment
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        % cleanup working memory
        clear call I*
    end
    
    fprintf('\nProcessed flight data in %02d:%02d:%02d\n', floor(toc/3600), floor(toc/60), floor(mod(toc,60)))
    
    % close waitbar
    close(wh)
    
    % cleanup working memory
    clear flight
end

