function varargout = mca_process_main(data,fs,PLOT_FLAG)
% MCA_PROCESS_MAIN  This function processes multi-harmonic multi-component
% time series data by separating individual components using the FrFT, then
% looking at the Hilbert estimates of instantaneous amplitude and frequency
% for each component
%
% IF = mca_process_main(data,fs) takes in a column vector, data, and sampling
%   rate, fs, and returns a cell struct of instantaneous frequencies for each
%   component.
% [IF,IA,IMF,COMPS,TIME,YHOLD,SUM_FRAC] = mca_process_main(data,fs) optionally
%   returns the following variables:
%
%A cell structure IF, which holds the instantaneous frequencies of each of
%    the two components - each component is decomposed into 4 modes
%A cell structure IA, which holds the instantaneous amplitudes of each of
%    the two components - each component is decomposed into 4 modes
%A cell structure IMF, which holds the decompositions of each component
%A data-length by 2 matrix COMPS holding the two time series resulting from 
%    the filtering.
%A vector TIME holding the time
%A matrix, YHOLD, holding the fractional FT of the data
%A matrix, SUM_FRAC, holding the summation of all fractional FT used to
%    pick the overlap of the first and second harmonics and the
%    correspondsing frequencies
%
% [...] = mca_process_main(data,fs,true) plots all of the steps in EMD analysis
% mca_process_main(...) simply turns on plotting without return args
% 

% force column vector
data = data(:);

% initialize vector
IMF = cell(2,1);

% turn on plotting if no output args or flag manually set
nout = nargout-2;
if nout<0 || nargin<3
    PLOT_FLAG = true;
end


% plotting params
%cMin = 0;
%cMax = 66;
nfft = 64;

%make signal odd length for algorithms
if ~mod(length(data),2)
    data = data(1:end-1);
end

%calculate the time vector for the sequence
ld = length(data);
time = ((0:ld-1)/fs)';



if PLOT_FLAG
    % plot spectrogram of original data
    fh1 = figure;
    [y,f,t,p] = spectrogram(data,nfft,round(nfft*0.8),[],fs*1e-3,'yaxis');
    surf(t,f,10*log10(abs(p)),'EdgeColor','none');
    axis xy; axis tight; colormap(jet); view(0,90);
    xlabel('Time (ms)')
    ylabel('Frequency (kHz)')
%    caxis([cMin cMax])
    colorbar
    set(fh1, 'Name', sprintf('Original Signal'));
    
end


%% some ideas regarding how to pick starting frequency

ind = 1;
%rotate the Fourier transform through time-frequency plane in fractions
frac_step = 0.5:0.002:1;
lind = length(frac_step);

%preallocate for all the fractional rotations
yhold = zeros(ld,lind);

%if plotting the fractional FT, do it in different colors
% if PLOT_FLAG
%     fh1 = figure;
%     hold on;grid
%     set(gca,'XLim',[0 lind],'YLim',[0 fs/2],'ZLim',[0 2*max(data)],'view',[-60 30])
%     %colors = colormap(jet(lind));
%     set(fh1, 'Name', sprintf('Stepwise Fractional Fourier Transform in incremments of 0.01'),'renderer','OpenGL');
% end

for k = frac_step
    yhold(:,ind) = (abs(frft((data(1:end)),k))).^2;
%     if PLOT_FLAG
%         figure(fh1)
%         %plot((-floor(ld/2):floor(ld/2))/floor(ld/2)*(fs/2),yhold(:,ind),'color',colors(ind,:)),hold on
%         plotcol(ones(floor(ld/2),1)'*ind,(0:floor(ld/2-1))/floor(ld/2)*(fs/2),yhold(ceil(ld/2+1):end,ind)',yhold(ceil(ld/2+1):end,ind)')
%     end
    ind = ind+1;
end

if PLOT_FLAG
    fh2 = figure;
    %surf(1:ind-1,(-floor(ld/2):floor(ld/2))/floor(ld/2)*(fs/2),yhold),shading interp   
    %view(0,-90),set(gca,'YDir','reverse','YLim',[0 fs/2])
    %can do above with imagesc
    imagesc(frac_step,(-floor(ld/2):floor(ld/2))/floor(ld/2)*(fs/2),(yhold),[min(yhold(:)) max(yhold(:))])
    set(gca,'YDir','normal','YLim',[0 fs/2])
    set(fh2, 'Name', sprintf('Fractional Fourier Transform'));
end

%now filter the sum of all frational FT using an LS-optimal (for polynomial fitting)
%filter to preserve high frequency content. The length of the filter is arbitrary 
%but the order fit is not - cubic
golay = sgolayfilt(sum(yhold(ceil(ld/2):end,:),2),3,21);

%get the second half of the summation as the first half is from the
%negative frequencies and should be rejected
%half_gol = golay(ceil(ld/2):end);
half_gol = golay;
%normalize on 0-1
min_gol = min(half_gol);
max_gol = max(half_gol-min_gol);
norm_gol = (half_gol-min_gol)/max_gol;

%set up the frequency vector for the frac FT
xloc = (0:floor(ld/2))/floor(ld/2)*(fs/2);

if PLOT_FLAG
    fh3 = figure;
    plot(xloc,norm_gol)
    set(gca,'XLim',[0 (fs/2)])
    set(fh3, 'Name', sprintf('Summation of Fractional Fourier Transforms'));
end

%find the peaks in the summation.  These will occur at roughly the center
%frequencies of the harmonics.  We are interested in the trough between the
%first and second peaks as this is the crossover between the first and
%second harmonics
[PKS,LOCS] = findpeaks(norm_gol,'minpeakdistance',20);
%crossover  = norm_gol==min(norm_gol(LOCS(1):LOCS(2)));
%fCutoff = xloc(crossover);
fCutoff = 55000;

% remove peaks below ~10kHz
LOCS(LOCS<50) = [];

%% preliminary filtering to get approximate isolation of the first
%% component -  THIS DOES NOT NEED TO BE PRECISE
[b,a] = butter(20,2*fCutoff/fs);
res = filtfilt(b,a,data);

if PLOT_FLAG
    % plot spectrogram of filtered data
    fh4 = figure;
    [y,f,t,p] = spectrogram(res,nfft,round(nfft*0.8),[],fs*1e-3,'yaxis');
    surf(t,f,10*log10(abs(p)),'EdgeColor','none');
    axis xy; axis tight; colormap(jet); view(0,90);
    xlabel('Time (ms)')
    ylabel('Frequency (kHz)')
%    caxis([cMin cMax])
    colorbar
    set(fh4, 'Name', sprintf('Rough Filtered Signal'));
    
end



%% Create a temporary IMF to hold temporary decomposition
[IMF_temp,IA_temp,IF_temp] = mca_extract(res,fs);

% normalize instantaneous amplitude for consistent threshold
IA_temp = IA_temp/max(IA_temp(:,1));

% find the beginning and end of the first harmonic for spline fitting
gamma = 0.10;    % select normalized threshold

% look for monotonicity in IA
[imn,imx] = minmaxloc(IA_temp(:,1));
diff_IA = diff(IA_temp(imx));

%keep only the positive
pos_idx = diff_IA>0;
pos_idx = logical([0;pos_idx]);

%go back and get those peaks in imx that correspond to a positive
%derivative
new_imx = imx(pos_idx);


f = find(IA_temp(new_imx,1)>gamma,1);
l = find(IA_temp(new_imx,1)>2*gamma,1,'last');

first = new_imx(f);
last  = new_imx(l);

%% now get the IF there and fit a quadratic - this should be sufficient
%if the segment to fit the polynomial is less than a third the size of the
%origninal data, just make it a linear fit
if last-first>ld/3
    PP = polyfit(time(first:last),IF_temp(first:last,1),2);
else
   PP = polyfit(time(first:last),IF_temp(first:last,1),1); 
end
%figure,plot(time(first:last),IF_temp(first:last,1)),return
%Can use centering and scaling but it takes much longer and doesn't always
%give better results - it will stop the errors though
%[PP,S,MU] = polyfit(time(first:last),IF_temp(first:last,1),2)
%% Should check goodness of fit and reduce the polynomial to a linear
%% approximation if it is no good.

%% generate a frequency demodulation  - this is similar to a phase function 
PF = polyval(PP,time);
%If you decide to use centering and scaling
%PF = polyval(PP,time,[],MU);
PF = PF(:);
% figure,plot(-diff(PF)),return
if PLOT_FLAG
    figure(fh4),hold on,plot3(time*1000,PF*1e-3,ones(length(time),1)*-50,'k','linewidth',6)
    CLim = get(gca,'CLim');
end
% size(IA_temp)
% figure,plot(time(2:end-1),IF_temp(:,1)),hold on,plot(time,50000*IA_temp,'r')
% figure,plot(time(first:last),IF_temp(first:last))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% filter out first component
first_comp = mca_iffilt(data,PF,fs);

if PLOT_FLAG
    % plot time series waveforms
    fh5 = figure;
    plot(time,real(data))
    hold on;
    plot(time,real(first_comp),'r')
    set(fh5, 'Name', sprintf('Original Signal and Isolated 1^{st} Component Time Series'));
    
    figure,spectrogram(real(first_comp),nfft,round(nfft*0.8),[],fs*1e-3,'yaxis');
    caxis([CLim(1) CLim(2)])
    colorbar
end

%% extract IF and IA for each IMF of 1st component
[IMF{1},IA{1},IF{1}] = mca_extract(first_comp,fs);

minIA = min(IA{1}(:));
maxIA = max(IA{1}(:) - minIA);
IA_norm{1} = IA{1}/maxIA;

if PLOT_FLAG
    fh6 = figure;
    plot(time,IMF{1})
    title('Intrinsic Mode Functions for the 1^{st} Component')
    set(fh6, 'Name', sprintf('IMF-1 Time Series'));
    
    fh7 = figure;
    plot(time*1000,IA_norm{1},'linewidth',2); grid on
    title('Instantaneous Amplitude for 1^{st} Component')
    ylabel('Amplitude (normalized)')
    xlabel('time (ms)')
    axis([0 time(end)*1000 0 1.1])
    set(fh7, 'Name', sprintf('IMF-1 Instantaneous Amplitude'));
    
    fh8 = figure;
    %plot(time(1:length(idx))*1000,IF{1}(:,1:3)*1e-3,'linewidth',2); grid on
    plot(time(1:length(IF{1}))*1e3,IF{1}(:,1:3)*1e-3,'linewidth',2); grid on
    title('Instantaneous Frequency for 1^{st} Component')
    ylabel('Frequency (kHz)')
    xlabel('time (ms)')
    axis([0 time(end)*1000 0 fs*1e-3/2])
    set(fh8, 'Name', sprintf('IMF-1 Instantaneous Frequency'));
    
    % look at spectral content of the first component
    fh9 = figure;
    [y,f,t,p] = spectrogram(first_comp,nfft,round(nfft*0.8),[],fs*1e-3,'yaxis');
    surf(t,f,10*log10(abs(p)),'EdgeColor','none');
    axis xy; axis tight; colormap(jet); view(0,90);
    set(gca,'Ylim',[0 fs*1e-3/2])
    xlabel('Time (ms)')
    ylabel('Frequency (kHz)')
    caxis([CLim(1) CLim(2)])
    colorbar
    set(fh9, 'Name', sprintf('IMF-1 Spectrogram'));
end



%% extract IF and IA for each IMF of 2nd component
second_comp = data-first_comp;
[IMF{2},IA{2},IF{2}] = mca_extract(second_comp,fs);

minIA = min(IA{2}(:));
maxIA = max(IA{2}(:) - minIA);
IA_norm{2} = IA{2}/maxIA;

if PLOT_FLAG
    fh10 = figure;
    [y,f,t,p] = spectrogram(second_comp,nfft,round(nfft*0.8),[],fs*1e-3,'yaxis'); 
    surf(t,f,10*log10(abs(p)),'EdgeColor','none');
    axis xy; axis tight; colormap(jet); view(0,90);
    set(gca,'Ylim',[0 fs*1e-3/2])
    xlabel('Time (ms)')
    ylabel('Frequency (kHz)')
    caxis([CLim(1) CLim(2)])
    colorbar
    set(fh10, 'Name', sprintf('IMF-2 Spectrogram'));
    
    fh11 = figure;
    plot(time,IMF{2})
    title('Intrinsic Mode Functions for the 2^{nd} Component')
    set(fh11, 'Name', sprintf('IMF-2 Time Series'));
    
    fh12 = figure;
    plot(time*1000,IA_norm{2},'linewidth',2); grid on
    title('Instantaneous Amplitudes for the 2^{nd} Component')
    ylabel('Amplitude (normalized)')
    xlabel('time (ms)')
    axis([0 time(end)*1000 0 1.1])
    set(fh12, 'Name', sprintf('IMF-2 Instantaneous Amplitudes'));
    
    fh13 = figure;
    %plot(time(1:length(idx))*1000,IF{2}(:,1:3)*1e-3,'linewidth',2); grid on
    plot(time(1:length(IF{2}))*1000,IF{2}(:,1:3)*1e-3,'linewidth',2); grid on
    title('Instantaneous Frequencies for the 2^{nd} Component')
    ylabel('Frequency (kHz)')
    xlabel('time (ms)')
    axis([0 time(end)*1000 0 fs*1e-3/2])
    set(fh13, 'Name', sprintf('IMF-2 Instantaneous Frequencies'));
    
    %figure(fh4),hold on,plot3(time*1000,IF2*1e-3,ones(length(time),1)*-50,'k','linewidth',6)
end

%%
%Unclear what the below code is doing - fitting spline function to IF/IA?

% %% find the beginning and end of the second harmonic for spline fitting
% %second harmonic is a fair bit noisier so push the IA limit out a bit
% 
% gamma = 0.15;
% first = find(IA_norm{2}(:,1)>gamma,1);
% last  = find(IA_norm{2}(:,1)>3*gamma,1,'last');
% 
% %% now get the IF there and fit a quadratic - second harmonic is less curved
% PP2 = polyfit(time(first:last),IF{2}(first:last,1),2);
% %Can use centering and scaling but it takes much longer and doesn't always
% %give better results - it will stop the errors though
% %[PP2,S2,MU2] = polyfit(time(first:last),IF{2}(first:last,1),2);
% 
% %% generate a frequency demodulation  - this is similar to a phase function 
% PF2 = polyval(PP2,time);
% PF2 = PF2(:);
% %If you decide to use centering and scaling
% %PF = polyval(PP2,time,[],MU2);
% if PLOT_FLAG
% figure(fh4),hold on,plot3(time*1000,PF2*1e-3,ones(length(time),1)*-50,'k','linewidth',6)
% end
% 
% 
% 
% 
% %%
% comps = [real(first_comp) real(second_comp)];
% SUM_FRAC = [xloc' norm_gol];

%% THIS IS FUTURE WORK STUFF
%demod2 = data.*exp(-j*2*pi*cumtrapz(time,PF2-5000));

% %% filter the demodulated signal 
% 
% [b,a] = butter(20,2*35e3/fs);
% resdemod2 = filtfilt(b,a,(demod2));
% 
% %% remodulate the residual
% remod2 = resdemod2.*exp(j*2*pi*PF.*time);



%%assign the optional outputs
switch nargout
    case 1
        varargout(1) = {IF};
    case 2
        varargout(1) = {IF};
        varargout(2) = {IA};
    case 3
        varargout(1) = {IF};
        varargout(2) = {IA};
        varargout(3) = {IMF};
    case 4 
        varargout(1) = {IF};
        varargout(2) = {IA};
        varargout(3) = {IMF};
        varargout(4) = {comps};
    case 5
        varargout(1) = {IF};
        varargout(2) = {IA};
        varargout(3) = {IMF};
        varargout(4) = {comps};
        varargout(5) = time;
    case 6
        varargout(1) = {IF};
        varargout(2) = {IA};
        varargout(3) = {IMF};
        varargout(4) = {comps};
        varargout(5) = {time};
        varargout(6) = {yhold};
    case 7
        varargout(1) = {IF};
        varargout(2) = {IA};
        varargout(3) = {IMF};
        varargout(4) = {comps};
        varargout(5) = {time};
        varargout(6) = {yhold};
        varargout(7) = {SUM_FRAC};
   
end

function [imn,imx]=minmaxloc(y)
%
% Locate indexes of minima and maxima
%
S=diff(y);
S1=S(2:end);
S2=S(1:end-1);
imx=find(S1.*S2 < 0 & S1-S2 < 0)+1;
imn=find(S1.*S2 < 0 & S1-S2 > 0)+1;


