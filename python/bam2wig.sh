#! /bin/zsh

for k in MCF7_t0005min_RNA_2012-03.bam MCF7_t0010min_RNA_2012-03.bam MCF7_t0020min_RNA_2012-03.bam MCF7_t0040min_RNA_2012-03.bam MCF7_t0080min_RNA_2012-03.bam MCF7_t0160min_RNA_2012-03.bam MCF7_t0320min_RNA_2012-03.bam MCF7_t0640min_RNA_2012-03.bam MCF7_t1280min_RNA_2012-03.bam ; do
    $HOME/bin/bamunique $k ${k/_/_unique_}
done

for k in MCF7_unique_t0005min_RNA_2012-03.bam MCF7_unique_t0010min_RNA_2012-03.bam MCF7_unique_t0020min_RNA_2012-03.bam MCF7_unique_t0040min_RNA_2012-03.bam MCF7_unique_t0080min_RNA_2012-03.bam MCF7_unique_t0160min_RNA_2012-03.bam MCF7_unique_t0320min_RNA_2012-03.bam MCF7_unique_t0640min_RNA_2012-03.bam MCF7_unique_t1280min_RNA_2012-03.bam ; do
    $HOME/bin/samtools sort -m 4000000000 $k ${k/.bam/.sorted}
    #echo sort -m 4000000000 $k ${k/.bam/.sorted}
done

for k in MCF7_unique_t0005min_RNA_2012-03.sorted.bam MCF7_unique_t0010min_RNA_2012-03.sorted.bam MCF7_unique_t0020min_RNA_2012-03.sorted.bam MCF7_unique_t0040min_RNA_2012-03.sorted.bam MCF7_unique_t0080min_RNA_2012-03.sorted.bam MCF7_unique_t0160min_RNA_2012-03.sorted.bam MCF7_unique_t0320min_RNA_2012-03.sorted.bam MCF7_unique_t0640min_RNA_2012-03.sorted.bam MCF7_unique_t1280min_RNA_2012-03.sorted.bam ; do
    $HOME/bin/rsem-bam2wig $k ${k/.sorted.bam/.wig} mywiggle --no-fractional-weight
done

