# Digenome-Detect
Digenome-seq analysis tool

## Summary
Digenome-Detect is a tool for detection of Digenome-seq cleavage site.

## Install

### Prerequisite
- OS: Linux
  - We tested this system on CentOS7, CentOS8 and Ubuntu 24.04 LTS
- Rust
  - Install Rust programming language from https://www.rust-lang.org/tools/install
- Java17
  - Install Java17 or later from https://jdk.java.net/java-se-ri/17
  - Set JAVA_HOME and PATH environment variables on your environment
### Build and Install
+ Download the source code from GitHub
```
git clone https://github.com/nihsmtgt/Digenome-Detect.git
```
+ Build and install
```
cd Digenome-detect/rust
cargo build
cp target/debug/digenome_seek /usr/local/bin
cd ..
mvn compile
mvn package
mkdir /path/to/install
cp target/digenome_detect-1-jar-with-dependencies.jar /path/to/install
```

## Typical Analysis Workflow

### Preparation
Digenome-Detect can analyze BAM files. Please note that the duplicated reads should be removed from BAM files before this analysis. 

  1. Map FastQ files to Reference Genome
```
/usr/local/bwa-0.7.15/bwa mem -t 20 -P -M \
    -R "@RG\tID:1\tSM:Sample_DGS_HAP_NT\tPL:ILLUMINA" \
    Homo_sapiens.GRCh38.dna_sm.primary_assembly_fix_20210506.fa.gz \
    DGS_HAP_NT_FDPL210060334-1a_HVKLYDSXY_L2_1.fq.gz  DGS_HAP_NT_FDPL210060334-1a_HVKLYDSXY_L2_2.fq.gz  |samtools view -Sb - -o DGS_HAP_NT_FDPL210060334-1a_HVKLYDSXY_L2_1.unsort.bam
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
　5．(Option) Realignment 
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
### analysis
```
java -jar /path/to/install/digenome_detect-1-jar-with-dependencies.jar digenome_detect.Main --bam Sample_DGS_HAP_NT.70x.bam --out Sample_DGS_HAP_NT.70x --threads 24 > log_NT.txt 2>&1
 ```
 This combination of options are used in the published article.

## Result
The result files are generated for each chromosome.
### Example of output
```
chr2    41525655        41525657        CLSCORE=20.56;DP=58.0;CS=0.85;Ratio=1.000;FISHER=48.728;RevHead=0;RevTail=16;FwdHead=6;FwdTail=0;MQ0=0;CLIPS=0;CONTROL=0,0,1,1;CONTROL_CLSCORE=1.20
chr2    41600665        41600667        CLSCORE=102.41;DP=46.0;CS=93.59;Ratio=1.000;FISHER=230.084;RevHead=0;RevTail=35;FwdHead=46;FwdTail=0;MQ0=0;CLIPS=0;CONTROL=0,0,1,0;CONTROL_CLSCORE=1.35
chr2    41629021        41629023        CLSCORE=4.22;DP=59.0;CS=0.00;Ratio=0.750;FISHER=6.854;RevHead=1;RevTail=3;FwdHead=4;FwdTail=1;MQ0=0;CLIPS=0;CONTROL=0,1,0,1;CONTROL_CLSCORE=1.32
chr2    41666711        41666713        CLSCORE=4.17;DP=77.0;CS=0.01;Ratio=0.889;FISHER=13.222;RevHead=0;RevTail=5;FwdHead=3;FwdTail=1;MQ0=0;CLIPS=0;CONTROL=1,1,1,1;CONTROL_CLSCORE=1.04
```
- CLSCORE: Cleavage Likelihood score
- DP: Read depth
- CS: Cleavage score of re-implmented digenome-toolkit method
- Ratio: (reverse tails + forward heads)/(reverse tails + forward heads + reverse heads + forward tails)
  - i.e.  = cleaved ends/all ends 
- FISHER: Strand bias similar to GATK
- RevHead, RevTail, FwdHead, FwdTail: count of these ends
- MQ0: MQ0 read count indicating repeats
- CLIPS: Count of clipped end
- CONTROL: When you use --control [bam file] option, these counts of ends at same position in control's BAM file will be reported.
- CONTROL_CLSCORE: When you use --control [bam file] option, the CLSCORE at same position in control's BAM file will be reported.

----
## Command-Line Options

### Required Options

- `--bam [file]`: Specifies the path to the input BAM file.
- `--out [prefix]`: Specifies the output prefix for generated files.

### Optional Options

- `--control [file]`: Specifies the path to the control BAM file.
- `--threads, -t [number]`: Specifies the number of threads to use. Default is `8`.
- `--regions, --region [comma-separated values]`: Specifies target regions like 'chr19:12345..23456'. Use `GRCm` to specify mouse chromosomes. 
- `--mq [number]`: Specifies the minimum mapping quality filter. Default is `0` (i.e. no filter).
- `--width [number]`: Specifies the detection width, which allows for the inclusion of cleaved sites as well as overhangs. Default is `3`.
- `--debug`: Enables debug mode.
- `--strandbias [true/false]`: Enables or disables Fisher's exact test for evaluating strand bias on cleaved sites.
- `--calc_cleavage_score [true/false]`: Enables or disables cleavage score calculation.
- `--help`: Displays help information.

### Examples

- To run the program with a specific BAM file and output prefix:

  ```bash
  java digenome_detect.Main --bam input.bam --out output_prefix
  ```

- To run the program with 4 threads and specified regions:

  ```bash
  java digenome_detect.Main --bam input.bam --out output_prefix --threads 4 --regions "chr1,chr2,chr3"
  ```

- To enable debug mode:

  ```bash
  java digenome_detect.Main --bam input.bam --out output_prefix --debug
  ```

## Notes

- Make sure to specify both `--bam` and `--out` options; otherwise, the program will terminate with an error.
- If `--regions` is not specified, it defaults to human chromosomes.

