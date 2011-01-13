chr_index=1;

addpath ~/synergy_data/PolII/Mapping_results
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


rnaseq0=read_mappingfile('/home/jaakkopeltonen/synergy_data/RNA/Mapping_results/MCF7_t0_unique.bed');
[bins0,binsstart0]=chromosome_mapping_histogram(rnaseq0,1,100,0);
rnaseq5=read_mappingfile('/home/jaakkopeltonen/synergy_data/RNA/Mapping_results/MCF7_t5_unique.bed');
[bins5,binsstart5]=chromosome_mapping_histogram(rnaseq5,1,100,0);
rnaseq10=read_mappingfile('/home/jaakkopeltonen/synergy_data/RNA/Mapping_results/MCF7_t10_unique.bed');
[bins10,binsstart10]=chromosome_mapping_histogram(rnaseq5,1,100,0);
rnaseq20=read_mappingfile('/home/jaakkopeltonen/synergy_data/RNA/Mapping_results/MCF7_t20_unique.bed');
[bins20,binsstart20]=chromosome_mapping_histogram(rnaseq5,1,100,0);
rnaseq40=read_mappingfile('/home/jaakkopeltonen/synergy_data/RNA/Mapping_results/MCF7_t40_unique.bed');
[bins40,binsstart40]=chromosome_mapping_histogram(rnaseq5,1,100,0);
rnaseq80=read_mappingfile('/home/jaakkopeltonen/synergy_data/RNA/Mapping_results/MCF7_t80_unique.bed');
[bins80,binsstart80]=chromosome_mapping_histogram(rnaseq5,1,100,0);
rnaseq160=read_mappingfile('/home/jaakkopeltonen/synergy_data/RNA/Mapping_results/MCF7_t160_unique.bed');
[bins160,binsstart160]=chromosome_mapping_histogram(rnaseq5,1,100,0);


figure;
subplot(17,1,1);
plot_pol2data(chr_index,[],peaks0,[],[],[],0,3200,0);
axis tight;

subplot(17,1,2);
plot_pol2data(chr_index,rnaseq0,[],bins0,binsstart0,100,0,3200,0);
axis tight;

subplot(17,1,3);
plot_pol2data(chr_index,[],peaks5,[],[],[],0,3200,0);
axis tight;

subplot(17,1,4);
plot_pol2data(chr_index,rnaseq5,[],bins5,binsstart5,100,0,3200,0);
axis tight;

subplot(17,1,5);
plot_pol2data(chr_index,[],peaks10,[],[],[],0,3200,0);
axis tight;

subplot(17,1,6);
plot_pol2data(chr_index,rnaseq10,[],bins10,binsstart10,100,0,3200,0);
axis tight;

subplot(17,1,7);
plot_pol2data(chr_index,[],peaks20,[],[],[],0,3200,0);
axis tight;

subplot(17,1,8);
plot_pol2data(chr_index,rnaseq20,[],bins20,binsstart20,100,0,3200,0);
axis tight;

subplot(17,1,9);
plot_pol2data(chr_index,[],peaks40,[],[],[],0,3200,0);
axis tight;

subplot(17,1,10);
plot_pol2data(chr_index,rnaseq40,[],bins40,binsstart40,100,0,3200,0);
axis tight;

subplot(17,1,11);
plot_pol2data(chr_index,[],peaks80,[],[],[],0,3200,0);
axis tight;

subplot(17,1,12);
plot_pol2data(chr_index,rnaseq80,[],bins80,binsstart80,100,0,3200,0);
axis tight;

subplot(17,1,13);
plot_pol2data(chr_index,[],peaks160,[],[],[],0,3200,0);
axis tight;

subplot(17,1,14);
plot_pol2data(chr_index,rnaseq160,[],bins160,binsstart160,100,0,3200,0);
axis tight;

subplot(17,1,15);
plot_pol2data(chr_index,[],peaks320,[],[],[],0,3200,0);
axis tight;

subplot(17,1,16);
plot_pol2data(chr_index,[],peaks640,[],[],[],0,3200,0);
axis tight;

subplot(17,1,17);
plot_pol2data(chr_index,[],peaks1280,[],[],[],0,3200,0);
axis tight;

