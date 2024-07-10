# ScoreCheck

ScoreCheck is a Java program designed to process and analyze results of Digenome-detect. 
It filters and analyzes data based on various quality metrics and outputs the results.

## Features
- Remove redundant scores that are adjacent.
- Filter false positives based on a number of read with low mapping quality and clipped read.
- Determine the CLSCORE threshold by comparing the data with negative control data.

## Requirements
- Java Development Kit (JDK) 11 or later
- Web Browser for visualization
  
## Build
```
javac ScoreCheck.java
```

## Usage
```
java ScoreCheck [options] <case.bed> <control.bed>
```

## Options
- --out [outputPath]: The directory to output files.
- --threshold [value]: Initial score threshold to filter out low-quality cleavage-likelihood scores (default: 7.0).
- --ratio_threshold [value]: Sample/control ratio are calculated for narrow score ranges and minimum score of the score range with sample/control ratio below this value becomes final CLSCORE threshold (default: 1).
- --debug: Enable debug mode to dump scores to log.

## Input Files
Two bed files for the sample and the negative control resulting from Digenome-detect_Main. 

## Output
The program outputs:

- Filtered data in BED format (case.bed and control.bed).
- Candidate off-target cleavage sites in BED format (result.bed).
- List of sample/control ratio of the score ranges (plotdata.bed).
- HTML visualization of the plot of score distribution (plot.html).
  
