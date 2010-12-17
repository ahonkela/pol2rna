chr_index=1;

peaks0=read_peak_summit_files('MCF7_No_Treat_PolII_unique_macs_peaks.bed','MCF7_No_Treat_PolII_unique_macs_summits.bed');

peaks5=read_peak_summit_files('MCF7_E2_5min_PolII_unique_macs_peaks.bed','MCF7_E2_5min_PolII_unique_macs_summits.bed');

peaks10=read_peak_summit_files('MCF7_E2_10min_PolII_rep_unique_macs_peaks.bed','MCF7_E2_10min_PolII_rep_unique_macs_summits.bed');

peaks20=read_peak_summit_files('MCF7_E2_20min_PolII_unique_macs_peaks.bed','MCF7_E2_20min_PolII_unique_macs_summits.bed');

peaks40=read_peak_summit_files('MCF7_E2_40min_PolII_unique_macs_peaks.bed','MCF7_E2_40min_PolII_unique_macs_summits.bed');

peaks80=read_peak_summit_files('MCF7_E2_80min_PolII_unique_macs_peaks.bed','MCF7_E2_80min_PolII_unique_macs_summits.bed');

peaks160=read_peak_summit_files('MCF7_E2_160min_PolII_unique_macs_peaks.bed','MCF7_E2_160min_PolII_unique_macs_summits.bed');

peaks320=read_peak_summit_files('MCF7_E2_320min_PolII_unique_macs_peaks.bed','MCF7_E2_320min_PolII_unique_macs_summits.bed');

peaks640=read_peak_summit_files('MCF7_E2_640min_PolII_unique_macs_peaks.bed','MCF7_E2_640min_PolII_unique_macs_summits.bed');

peaks1280=read_peak_summit_files('MCF7_E2_1280min_PolII_unique_macs_peaks.bed','MCF7_E2_1280min_PolII_unique_macs_summits.bed');

figure;
subplot(10,1,1);
plot_pol2data(chr_index,[],peaks0,[],[],[],0,3200,0);
axis tight;

subplot(10,1,2);
plot_pol2data(chr_index,[],peaks5,[],[],[],0,3200,0);
axis tight;

subplot(10,1,3);
plot_pol2data(chr_index,[],peaks10,[],[],[],0,3200,0);
axis tight;

subplot(10,1,4);
plot_pol2data(chr_index,[],peaks20,[],[],[],0,3200,0);
axis tight;

subplot(10,1,5);
plot_pol2data(chr_index,[],peaks40,[],[],[],0,3200,0);
axis tight;

subplot(10,1,6);
plot_pol2data(chr_index,[],peaks80,[],[],[],0,3200,0);
axis tight;

subplot(10,1,7);
plot_pol2data(chr_index,[],peaks160,[],[],[],0,3200,0);
axis tight;

subplot(10,1,8);
plot_pol2data(chr_index,[],peaks320,[],[],[],0,3200,0);
axis tight;

subplot(10,1,9);
plot_pol2data(chr_index,[],peaks640,[],[],[],0,3200,0);
axis tight;

subplot(10,1,10);
plot_pol2data(chr_index,[],peaks1280,[],[],[],0,3200,0);
axis tight;

