# ScoreCheck

ScoreCheck is a Java program designed to process and analyze results of Digenome-Detect. 
It filters and analyzes data based on various quality metrics and outputs the results.

## Features
- Filter data based on a score threshold.
- Remove redundant scores that are adjacent.
- Filter out low-quality scores.
- Group scores by their value for both cases and controls.
- Generate graphical visualizations of the score distributions.

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
Two bed files of result of Digenome-Detect for RGEN and CONTROL. 

## Output
The program outputs:

- Filtered data in BED format.
- HTML visualization of the plot of score distribution.
  
