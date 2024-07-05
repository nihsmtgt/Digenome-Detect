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
- --out [outputPath]: Specify the output file path.
- --threshold [value]: Set the score threshold to filter out low-quality cleavage-likelihood scores (default: 7.0).
- --debug: Enable debug mode to dump scores to log.

## Input Files
Two bed files for the sample and the negative control resulting from Digenome-detect_Main. 

## Output
The program outputs:

- Filtered data in BED format.
- HTML visualization of the plot of score distribution.
  
