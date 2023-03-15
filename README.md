# Digenome-detect
Digenome-Seq analysis tool

## Summary
digenome-detect is a tool for detection of Digenome-Seq cleavage site.

## Typical Analysis Workflow

### Preparation

  1. Map FastQ files to Reference Genome
```
usr/local/bwa-0.7.15/bwa mem -t 20 -P -M \
    -R "@RG\tID:1\tSM:Sample_DGS_HAP_NT\tPL:ILLUMINA" \
    Homo_sapiens.GRCh38.dna_sm.primary_assembly_fix_20210506.fa.gz \
    DGS_HAP_NT_FDPL210060334-1a_HVKLYDSXY_L2_1.fq.gz  DGS_HAP_NT_FDPL210060334-1a_HVKLYDSXY_L2_2.fq.gz
```
  2. Sort BAM file 
```
samtools sort -T . -m 12G --threads 4 --output-fmt BAM \
    -o Sample_DGS_HAP_NT.sort.bam \
    DGS_HAP_NT_FDPL210060334-1a_HVKLYDSXY_L2_1.unsort.bam
```
  3. (Option) Run DepthOfCoverage for normalization if you want to normalize read depth for performance evaluation
```
java -X32g -jar GenomeAnalysisTK.jar -T DepthOfCoverage \
   -I Sample_DGS_HAP_NT.sort.bam \
   -L  wgs_calling_regions.hg38.interval_list \
   -R  Homo_sapiens.GRCh38.dna_sm.primary_assembly_fix_20210506.fa \
   -o  DepthOfCoverage.norm \
   -omitIntervals -omitBaseOutput
```
  4. De-duplication
  ```
java -Djava.io.tmpdir=./tmp -Duse_async_io_write_samtools=true -Xmx8g \
      -jar picard-2.26.11.jar MarkDuplicates \
      I=Sample_DGS_HAP_NT.sort.bam \
      O=Sample_DGS_HAP_NT.dedup.bam \
      M=Sample_DGS_HAP_NT.dupe_metrics.txt \
      CREATE_INDEX=true \
      OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
      REMOVE_DUPLICATES=true
```
　5．Realignment
 ```
 java -X32g -jar GenomeAnalysisTK.jar \
   -T RealignerTargetCreator -I \Sample_DGS_HAP_NT.dedup.bam \
   -R Homo_sapiens.GRCh38.dna_sm.primary_assembly_fix_20210506.fa \
   -disable_auto_index_creation_and_locking_when_reading_rods \
   -nt 4 \
   -o Sample_DGS_HAP_NT.interval_list \
   -known Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
   -known dbsnp_146.hg38.vcf.gz
 java -X32g -jar GenomeAnalysisTK.jar \
   -T IndelRealigner -I Sample_DGS_HAP_NT.dedup.bam \
   -R Homo_sapiens.GRCh38.dna_sm.primary_assembly_fix_20210506.fa \
   -disable_auto_index_creation_and_locking_when_reading_rods \
   -targetIntervals Sample_DGS_HAP_NT.interval_list \
   -o Sample_DGS_HAP_NT.realign.bam \
   -filterNoBases
 ```
