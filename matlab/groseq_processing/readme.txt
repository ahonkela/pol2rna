Steps to lift over hg18 wiggle file to hg19 wiggle file:
largely following instructions at https://www.biostars.org/p/81185/


Step 1: Download wiggle files from
http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE27463


Step 2: Download tools from 
http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/
Download the following tools: wigToBigWig , fetchChromSizes.sh , bigWigToBedGraph , liftOver


Step 3: Download liftover data file hg18ToHg19.over.chain from
http://hgdownload.cse.ucsc.edu/goldenPath/hg18/liftOver/


Step 4: Download chromosome sizes for HG18 and HG19 
/bin/sh fetchChromSizes.sh hg18 > hg18.chrom.sizes
/bin/sh fetchChromSizes.sh hg19 > hg19.chrom.sizes


Step 5: Remove unknown chromosome lines from WIG files. The .wig files contain references to 
- CAB chromosome (chromosome aberration)
- LTP6 (a gene???)
- RCBL (a gene???)
- RCA (a gene???)
- RCP1 (a gene???)
- TIM (a gene???)
Manually remove them from the files, name the modified files e.g. E20m_Minus_noCAB.wig


Step 6: convert WIG files to BIGWIG files
./wigToBigWig E20m_Minus_noCAB.wig hg18.chrom.sizes E20m_Minus.bw
./wigToBigWig E20m_Plus_noCAB.wig hg18.chrom.sizes E20m_Plus.bw
./wigToBigWig E210m_Minus_noCAB.wig hg18.chrom.sizes E210m_Minus.bw
./wigToBigWig E210m_Plus_noCAB.wig hg18.chrom.sizes E210m_Plus.bw
./wigToBigWig E2160m_Minus_noCAB.wig hg18.chrom.sizes E2160m_Minus.bw
./wigToBigWig E2160m_Plus_noCAB.wig hg18.chrom.sizes E2160m_Plus.bw
./wigToBigWig E240m_Minus_noCAB.wig hg18.chrom.sizes E240m_Minus.bw
./wigToBigWig E240m_Plus_noCAB.wig hg18.chrom.sizes E240m_Plus.bw


Step 7: convert BIGWIG files to BEDGRAPH files
./bigWigToBedGraph E20m_Minus.bw E20m_Minus.bedGraph
./bigWigToBedGraph E20m_Plus.bw E20m_Plus.bedGraph
./bigWigToBedGraph E210m_Minus.bw E210m_Minus.bedGraph
./bigWigToBedGraph E210m_Plus.bw E210m_Plus.bedGraph
./bigWigToBedGraph E2160m_Minus.bw E2160m_Minus.bedGraph
./bigWigToBedGraph E2160m_Plus.bw E2160m_Plus.bedGraph
./bigWigToBedGraph E240m_Minus.bw E240m_Minus.bedGraph
./bigWigToBedGraph E240m_Plus.bw E240m_Plus.bedGraph


Step 8: liftover BEDGRAPH files from HG18 to HG19
./liftOver E20m_Minus.bedGraph hg18ToHg19.over.chain E20m_Minus_hg19.bedGraph E20m_Minus_hg19_unmapped.bedGraph
./liftOver E20m_Plus.bedGraph hg18ToHg19.over.chain E20m_Plus_hg19.bedGraph E20m_Plus_hg19_unmapped.bedGraph
./liftOver E210m_Minus.bedGraph hg18ToHg19.over.chain E210m_Minus_hg19.bedGraph E210m_Minus_hg19_unmapped.bedGraph
./liftOver E210m_Plus.bedGraph hg18ToHg19.over.chain E210m_Plus_hg19.bedGraph E210m_Plus_hg19_unmapped.bedGraph
./liftOver E2160m_Minus.bedGraph hg18ToHg19.over.chain E2160m_Minus_hg19.bedGraph E2160m_Minus_hg19_unmapped.bedGraph
./liftOver E2160m_Plus.bedGraph hg18ToHg19.over.chain E2160m_Plus_hg19.bedGraph E2160m_Plus_hg19_unmapped.bedGraph
./liftOver E240m_Minus.bedGraph hg18ToHg19.over.chain E240m_Minus_hg19.bedGraph E240m_Minus_hg19_unmapped.bedGraph
./liftOver E240m_Plus.bedGraph hg18ToHg19.over.chain E240m_Plus_hg19.bedGraph E240m_Plus_hg19_unmapped.bedGraph


Step 9: download UTR regions from web based historical ENSEMBL biomart
browser following this explanation: https://www.biostars.org/p/4869/

(go to http://jul2012.archive.ensembl.org/index.html , select biomart,
select ensembl gene id, chromosome name, UTR3 coordinates and gene
coordinates and strand, save as ensemblv68_utrs.txt)


Step 10: remove sequence information from the resulting file, retain
only the header lines.

egrep "^[>]" ensemblv68_utrs_v2.txt > ensemblv68_utrs_headers_v2.txt


Step 11:
in matlab, read gene and UTR3 information using "analyze_utrs.m"


Step 12:
in matlab, read in each BEDGRAPH file and save it in matlab format
tempbed=read_bedgraph('E20m_Minus_hg19.bedGraph');
save E20m_Minus_hg19.mat tempbed -mat
tempbed=read_bedgraph('E20m_Plus_hg19.bedGraph');
save E20m_Plus_hg19.mat tempbed -mat
tempbed=read_bedgraph('E210m_Minus_hg19.bedGraph');
save E210m_Minus_hg19.mat tempbed -mat
tempbed=read_bedgraph('E210m_Plus_hg19.bedGraph');
save E210m_Plus_hg19.mat tempbed -mat
tempbed=read_bedgraph('E240m_Minus_hg19.bedGraph');
save E240m_Minus_hg19.mat tempbed -mat
tempbed=read_bedgraph('E240m_Plus_hg19.bedGraph');
save E240m_Plus_hg19.mat tempbed -mat
tempbed=read_bedgraph('E2160m_Minus_hg19.bedGraph');
save E2160m_Minus_hg19.mat tempbed -mat
tempbed=read_bedgraph('E2160m_Plus_hg19.bedGraph');
save E2160m_Plus_hg19.mat tempbed -mat


Step 13:
in matlab,
load E20m_Minus_hg19.mat
geneprof = make_geneprofiles_v2(finalgenes,tempbed);
save E20m_Minus_hg19_geneprof_v2.mat geneprof -mat
load E20m_Plus_hg19.mat
geneprof = make_geneprofiles_v2(finalgenes,tempbed);
save E20m_Plus_hg19_geneprof_v2.mat geneprof -mat

load E210m_Minus_hg19.mat
geneprof = make_geneprofiles_v2(finalgenes,tempbed);
save E210m_Minus_hg19_geneprof_v2.mat geneprof -mat
load E210m_Plus_hg19.mat
geneprof = make_geneprofiles_v2(finalgenes,tempbed);
save E210m_Plus_hg19_geneprof_v2.mat geneprof -mat

load E240m_Minus_hg19.mat
geneprof = make_geneprofiles_v2(finalgenes,tempbed);
save E240m_Minus_hg19_geneprof_v2.mat geneprof -mat
load E240m_Plus_hg19.mat
geneprof = make_geneprofiles_v2(finalgenes,tempbed);
save E240m_Plus_hg19_geneprof_v2.mat geneprof -mat

load E2160m_Minus_hg19.mat
geneprof = make_geneprofiles_v2(finalgenes,tempbed);
save E2160m_Minus_hg19_geneprof_v2.mat geneprof -mat
load E2160m_Plus_hg19.mat
geneprof = make_geneprofiles_v2(finalgenes,tempbed);
save E2160m_Plus_hg19_geneprof_v2.mat geneprof -mat


Step 14:
If desired, the individual profile files for different time points and
strands can be combined into a single cell array. In matlab,

allprof=cell(size(finalgenes,1),4,2);
load E20m_Minus_hg19_geneprof_v2.mat
t=1;s=1;for k=1:size(finalgenes,1), allprof{k,t,s}=geneprof{k}; end;
load E20m_Plus_hg19_geneprof_v2.mat
t=1;s=2;for k=1:size(finalgenes,1), allprof{k,t,s}=geneprof{k}; end;
load E210m_Minus_hg19_geneprof_v2.mat
t=2;s=1;for k=1:size(finalgenes,1), allprof{k,t,s}=geneprof{k}; end;
load E210m_Plus_hg19_geneprof_v2.mat
t=2;s=2;for k=1:size(finalgenes,1), allprof{k,t,s}=geneprof{k}; end;
load E240m_Minus_hg19_geneprof_v2.mat
t=3;s=1;for k=1:size(finalgenes,1), allprof{k,t,s}=geneprof{k}; end;
load E240m_Plus_hg19_geneprof_v2.mat
t=3;s=2;for k=1:size(finalgenes,1), allprof{k,t,s}=geneprof{k}; end;
load E2160m_Minus_hg19_geneprof_v2.mat
t=4;s=1;for k=1:size(finalgenes,1), allprof{k,t,s}=geneprof{k}; end;
load E2160m_Plus_hg19_geneprof_v2.mat
t=4;s=2;for k=1:size(finalgenes,1), allprof{k,t,s}=geneprof{k}; end;
save allprof.mat allprof finalgenes -mat







