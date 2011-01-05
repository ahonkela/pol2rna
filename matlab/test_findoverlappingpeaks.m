peaks0=read_peak_summit_files('PolII_peaks_MCF7_E2_No_treat_PolII_unique.bam_peaks.bed','PolII_peaks_MCF7_E2_No_treat_PolII_unique.bam_summits.bed');
peaks5=read_peak_summit_files('PolII_peaks_MCF7_E2_5min_PolII_unique.bam_peaks.bed','PolII_peaks_MCF7_E2_5min_PolII_unique.bam_summits.bed');
peaks10=read_peak_summit_files('PolII_peaks_MCF7_E2_10min_PolII_unique.bam_peaks.bed','PolII_peaks_MCF7_E2_10min_PolII_unique.bam_summits.bed');
peaks20=read_peak_summit_files('PolII_peaks_MCF7_E2_20min_PolII_unique.bam_peaks.bed','PolII_peaks_MCF7_E2_20min_PolII_unique.bam_summits.bed');
peaks40=read_peak_summit_files('PolII_peaks_MCF7_E2_40min_PolII_unique.bam_peaks.bed','PolII_peaks_MCF7_E2_40min_PolII_unique.bam_summits.bed');
peaks80=read_peak_summit_files('PolII_peaks_MCF7_E2_80min_PolII_unique.bam_peaks.bed','PolII_peaks_MCF7_E2_80min_PolII_unique.bam_summits.bed');
peaks160=read_peak_summit_files('PolII_peaks_MCF7_E2_160min_PolII_unique.bam_peaks.bed','PolII_peaks_MCF7_E2_160min_PolII_unique.bam_summits.bed');
peaks320=read_peak_summit_files('PolII_peaks_MCF7_E2_320min_PolII_unique.bam_peaks.bed','PolII_peaks_MCF7_E2_320min_PolII_unique.bam_summits.bed');
peaks640=read_peak_summit_files('PolII_peaks_MCF7_E2_640min_PolII_unique.bam_peaks.bed','PolII_peaks_MCF7_E2_640min_PolII_unique.bam_summits.bed');
peaks1280=read_peak_summit_files('PolII_peaks_MCF7_E2_1280min_PolII_unique.bam_peaks.bed','PolII_peaks_MCF7_E2_1280min_PolII_unique.bam_summits.bed');

allpeaks={peaks0,peaks5,peaks10,peaks20,peaks40,peaks80,peaks160,peaks320,peaks640,peaks1280};

chr_index=1;
peakextentmultiplier=1.2;
overlaps=find_overlapping_peaksovertime(chr_index,allpeaks,peakextentmultiplier);

% example time series of heights of the overlapping peaks
plot(overlaps{8997}{9})

