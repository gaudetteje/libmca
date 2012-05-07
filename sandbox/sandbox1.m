fname = 'chain27';

ts.data = niLoadBin([fname '.bin'], '.', 2);
ts.params = niLoadPrm([fname '.prm'], '.');
ts.fs = 5e5;
ts.time = (1:length(ts.data)).'./ts.fs;
%ts.data = ts.data./max(ts.data(:));        % normalize data
ts.data(:,2) = ts.data(:,2) - mean(ts.data(:,2));

idx = (3200:6000);

figure(1)
plot(ts.time(idx),ts.data(idx,2),'b')
hold on;
grid on;

figure
spectrogram(ts.data(idx,2),512,500,512,ts.fs,'yaxis');
title(fname)
set(gca,'clim',[-40 10])
colorbar
colormap jet

%%
figure;
tfrpwv(hilbert(ts.data(idx,2)));

