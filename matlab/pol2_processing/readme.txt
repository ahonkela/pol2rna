To generate the necessary data files for GP modeling, use the following steps. 
Note: you may need to edit file paths mentioned in the code to correspond to 
locations where you have stored the code and data files.

1. Get the data files transcript_counts.txt and transcript_exonpos.txt

2. Run convert_exonpos.m to generate transcript_exonpos.mat

3. Run check_exonstarts_part1.m and check_exonstarts_part2.m to generate problemgenes.mat

4. Run create_newbininfo.R and concatenate the output files into genes_dec2012_all.txt

5. Run create_newbininfo_part2.m to create bininfo_dec2012.mat

6. Get the data files: 
MCF7_NoTreat_PolII_2012-03_unique.bed.gz
MCF7_5min_PolII_2012-03_unique.bed.gz
MCF7_10min_PolII_2012-03_unique.bed.gz
MCF7_20min_PolII_2012-03_unique.bed.gz
MCF7_40min_PolII_2012-03_unique.bed.gz
MCF7_80min_PolII_2012-03_unique.bed.gz
MCF7_160min_PolII_2012-03_unique.bed.gz
MCF7_320min_PolII_2012-03_unique.bed.gz
MCF7_640min_PolII_2012-03_unique.bed.gz
MCF7_1280min_PolII_2012-03_unique.bed.gz

7. Run compute_pol2_binactivity_corrected_dec2012.m (which requires read_mappingfile.c, compute_pol2activityovergenes_c.c)
and uncomment the first part of the code (comment out the rest), to create
'pol0_2012_03.mat','pol5_2012_03.mat','pol10_2012_03.mat',
'pol20_2012_03.mat','pol40_2012_03.mat','pol80_2012_03.mat',
'pol160_2012_03.mat','pol320_2012_03.mat','pol640_2012_03.mat',
'pol1280_2012_03.mat'

8. Comment out the first part of the code in compute_pol2_binactivity_corrected_dec2012.m and uncomment the rest,
to generate all_gene_pol2bins_2013_01_02.mat and bininfo_dec2012_corrected.mat

9. Get the data file empty_regions.txt

10. Run compute_emptybins_pol2_corrected_sept2012.m (requires read_stringfile.m)
to generate all_empty_pol2bins_2012_09.mat

11. Run compute_pol2_summaryseries_corrected_dec2012.m to generate
pol2_summaryseries_2013_01_02.mat

