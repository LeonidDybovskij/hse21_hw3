[This README in russian](README.ru.md)

The purpose of this work was comparison of the RNA-seq data between two groups of mouse embryonic fibroblasts using DESeq2.

Three test samples (SRR3414629.fastq, SRR3414630.fastq, SRR3414631.fastq) and three control samples (SRR3414635.fastq, SRR3414636.fastq, SRR3414637.fastq) were taken for analysis.
The quality of the reads was checked using FastQC (the result is in the fastqc folder).

Next, genome mapping was carried out using HISAT2:
hisat2 -p 3 -x mm10/genome -U SRR3414636_1.fastq -S SRR3414636_1.sam  2>  SRR3414636.hisat

After that the number of uniquely mapped reads was calculated:
grep -P '^@|NH:i:1$' SRR3414636_1.sam > SRR3414636.uniq.sam
grep -v '^@' SRR3414636.uniq.sam | wc -l

Using HTSeq, the number of reads that fell on each gene was calculated too:
htseq-count --format=sam --stranded=no SRR3414636.uniq.sam  gencode.vM25.annotation.gtf > SRR3414636.counts (reads folder)

The resulting files show the number of reads corresponding to genome regions where no exons are annotated (__no_feature), as well as the number of reads that may belong to different genes (__ambiguous). Subtracting these reads, the following number of reads was obtained:

| File with reads  | Total reads | Subtracting the above reads |
| ------------- | ------------- | ------------- |
| SRR3414629.uniq.sam  | 18375888  | 16049609 |
| SRR3414630.uniq.sam  | 13186139  | 11465324 |
| SRR3414631.uniq.sam  | 20928945  | 18408851 |
| SRR3414635.uniq.sam  | 18428317  | 16275997 |
| SRR3414636.uniq.sam  | 17825380  | 15757580 |
| SRR3414637.uniq.sam  | 17844858  | 15736978 |

All .counts files were merged into the ALL.counts file (data folder).

Analysis with the DESeq2 package was carried out in python (src folder).
