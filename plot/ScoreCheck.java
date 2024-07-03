import java.io.*;
import java.util.*;
import java.util.stream.Collectors;

class ScoreCheck {
    static double scoreThreshold = 7.0;
    String outputPath = null;
    boolean separatedGraph = false;
    public static class Score implements Comparable<Score>{
        String chrom;
        double score = 0.0;
        int start = 0;
        int end = 0;
        int clips = 0;
        int mq0 = 0;
        int fwdHead = 0;
        int revTail = 0;
        boolean isControl = false;
        boolean toRemove = false;
        public Score(String c, int s, int e, double cs, int cl, int mq, int fh, int rt, boolean isControl){
            chrom = c;
            start = s;
            end = e;
            score = cs;
            mq0 = mq;
            clips = cl;
            fwdHead = fh;
            revTail = rt;
            this.isControl = isControl;
        }
        // for sorting by chromosome position order
        public int compareTo(Score s){
            if(this.chrom.equals(s.chrom)){
                return this.start - s.start;
            }else{
                return this.chrom.compareTo(s.chrom);
            }
        }
    }
    ArrayList<Score> cases = new ArrayList<>();
    ArrayList<Score> controls = new ArrayList<>();
    public ScoreCheck(String casePath, String controlPath, String outputPath) throws IOException{
        this.cases = loadBed(casePath, false);
        System.err.println("Cases: " + this.cases.size() + " points");
        this.controls = loadBed(controlPath, true);
        System.err.println("Controls: " + this.controls.size() + " points");
        this.outputPath = outputPath;
    }
    String getFileBasename(String path){
        String[] tokens = path.split("/");
        String filename = tokens[tokens.length-1];
        return filename.replaceAll(".txt", "");
    }
    public static ArrayList<Score> loadBed(String path, boolean isControl) throws IOException{
       ArrayList<Score> scores = new ArrayList<>();
       try ( BufferedReader br = new BufferedReader(new FileReader(path))) {
            String line = "";
            while((line = br.readLine()) != null){
                String[] tokens = line.split("\t");
                String chrom = tokens[0];
                int start = Integer.parseInt(tokens[1]);
                int end = Integer.parseInt(tokens[2]);
                HashMap<String, String> info = new HashMap<>();
                String[] infos = tokens[3].split(";");
                for(int i = 0; i < infos.length; i++){
                    String[] kv = infos[i].split("=");
                    info.put(kv[0], kv[1]);
                }
                try {
                    double clscore = Double.parseDouble(info.get("CLSCORE"));
                    // filter: MQ0+ CLIPS > FwdHead or MQ0 + CLIPS > RevTail or CLSCORE < 6.0
                    int clips = Integer.parseInt(info.get("CLIPS"));
                    int mq0 = Integer.parseInt(info.get("MQ0"));
                    int fwdHead = Integer.parseInt(info.get("FwdHead"));
                    int revTail = Integer.parseInt(info.get("RevTail"));
                    // if(mq0 + clips > fwdHead || mq0 + clips > revTail) continue;

                    Score s = new Score(chrom, start, end, clscore, clips, mq0, fwdHead, revTail, isControl);
                    scores.add(s);
                }catch(Exception e){
                    System.err.println("ERROR with: " +  line);
                    e.printStackTrace();
                    System.exit(1);
                }
            }
        }
        return scores;
    }
    // filter redundant scores between samples
    public static ArrayList<Score> filterRedundantScores2(ArrayList<Score> scores){
        ArrayList<Score> filtered = new ArrayList<>();
        // remove adjuscent scores located within 10 bases
        for(int i = 0; i < scores.size() - 1; i++){
            Score s1 = scores.get(i);
            Score s2 = scores.get(i+1);
            // remove both scores;
            if(s1.chrom.equals(s2.chrom) && Math.abs(s2.start - s1.start) < 10){
                // remove lower score
                System.err.println("removing redundant " + s1.chrom + ":" + s1.start + "-" + s1.end + "\t" + s1.score);
                scores.get(i).toRemove = true;
                scores.get(i+1).toRemove = true;
            }
        }
        for(Score s : scores){
            if(!s.toRemove){
                filtered.add(s);
            }
        }
        return filtered;
    }

    // filter redundant scores within a sample
    public static ArrayList<Score> filterRedundantScores3(ArrayList<Score> scores){
        // sort by chromosome position order
        Collections.sort(scores);
        
        ArrayList<Score> filtered = new ArrayList<>();
        int highestpos = 0;
        boolean redundant = false;
        // remove adjuscent scores located within 10 bases
        for(int i = 0; i < scores.size() - 1; i++){
            Score s1 = scores.get(i);
            Score s2 = scores.get(i+1);
            // remove both scores;
            if(s1.chrom.equals(s2.chrom) && Math.abs(s2.start - s1.start) < 10){
                if (redundant) {
                    if (s2.score > scores.get(highestpos).score) {
                        scores.get(highestpos).toRemove = true;
                        highestpos = i+1;                      
                    }
                    else{
                        scores.get(i+1).toRemove = true;
                    }
                }
                else{
                // remove lower score
                System.err.println("removing redundant " + s1.chrom + ":" + s1.start + "-" + s1.end + "\t" + s1.score);
                if (s1.score > s2.score) {
                    scores.get(i+1).toRemove = true;
                    highestpos =  i;                   
                }
                else{
                    scores.get(i).toRemove = true;
                    highestpos = i+1;
                }
                redundant = true;
                }
            }
            else{
                redundant = false;
            }
        }
        for(Score s : scores){
            if(!s.toRemove){
                filtered.add(s);
            }
        }
        return filtered;
    }

    // filter redundant scores within a sample
    public static ArrayList<Score> filterRedundantScores(ArrayList<Score> scores){
        // sort by chromosome position order
        Collections.sort(scores);

        // remove adjuscent scores located within 10 bases
        boolean redundant = true;
        while(redundant){
            redundant = false;
            System.err.println("checking redundant scores");
            inner_loop:
            for(int i = 0; i < scores.size() - 1; i++){
                Score s1 = scores.get(i);
                Score s2 = scores.get(i+1);
                if(s1.chrom.equals(s2.chrom) && Math.abs(s2.start - s1.start) < 10){
                    // remove lower score
                    System.out.println("removing redundant scores");
                    if (s1.score > s2.score) {
                        System.err.println("removing " + s2.chrom + ":" + s2.start + "-" + s2.end + "\t" + s2.score);                                            
                        scores.remove(i+1);
                    } else {
                        System.err.println("removing " + s1.chrom + ":" + s1.start + "-" + s1.end + "\t" + s1.score);
                        scores.remove(i);
                    }
                    redundant = true;
                    break inner_loop;
                }
            }
        }
        // sort by score desc
        /*
            Collections.sort(scores, new Comparator<Score>(){
            @Override
            public int compare(Score s1, Score s2){
                if(s1.score > s2.score) return -1;
                else if(s1.score < s2.score) return 1;
                else return 0;
            }
        }); */
        return scores;    
    }
    // filter low quality scores
    public static ArrayList<Score> filterLowQualityScores(ArrayList<Score> scores, double threshold){
        ArrayList<Score> filtered = new ArrayList<>();
        for(Score s : scores){
            if(s.score > threshold){
                filtered.add(s);
            }
        }
        return filtered;
    }
    public static Map<Double, List<Score>> groupScoresByValue(ArrayList<Score> scores) {
        Map<Double, List<Score>> groupedScores = new HashMap<>();
        for (Score s : scores) {
            if (!groupedScores.containsKey(s.score)) {
                groupedScores.put(s.score, new ArrayList<>());
            }
            groupedScores.get(s.score).add(s);
        }
        return groupedScores;
    }
    public void run() throws IOException {
        // filter low quality scores
        this.cases = filterLowQualityScores(this.cases, scoreThreshold);
        this.controls = filterLowQualityScores(this.controls, scoreThreshold);
        // filter redundant(adjacent) scores
        this.cases = filterRedundantScores3(this.cases);
        this.controls = filterRedundantScores3(this.controls);
        // merge data
        ArrayList<Score> merged = new ArrayList<>();
        merged.addAll(this.cases);
        merged.addAll(this.controls);
        // sort by position
        Collections.sort(merged);
        // filter redundant(adjacent) scores
        merged = filterRedundantScores2(merged);
        // split and filter by if(mq0 + clips > fwdHead || mq0 + clips > revTail) 
        this.cases = new ArrayList<>();
        this.controls = new ArrayList<>();
        for(Score s : merged){
            if (s.mq0 + s.clips > s.fwdHead || s.mq0 + s.clips > s.revTail) {
                continue;
            }
            if(s.isControl){
                this.controls.add(s);
            }else{
                this.cases.add(s);
            }
        }
        // sort by score desc
        Collections.sort(this.cases, new Comparator<Score>(){
            @Override
            public int compare(Score s1, Score s2){
                if(s1.score > s2.score) return -1;
                else if(s1.score < s2.score) return 1;
                else return 0;
            }
        });
        // output filtered data to files
        // sort by chromosome position order
        Collections.sort(this.cases);
        Collections.sort(this.controls);
        try (
            FileWriter fw = new FileWriter("cases.bed");
            BufferedWriter bw = new BufferedWriter(fw);
            PrintWriter pw = new PrintWriter(bw))
        {
            // print header
            pw.println("#chrom\tstart\tend\tscore\tclips\tmq0\tfwdHead\trevTail");
            for(Score s : this.cases){
                pw.println(s.chrom + "\t" + s.start + "\t" + s.end + "\t" + s.score + "\t" + s.clips + "\t" + s.mq0 + "\t" + s.fwdHead + "\t" + s.revTail);
            }
        }
        try (
            FileWriter fw = new FileWriter("controls.bed");
            BufferedWriter bw = new BufferedWriter(fw);
            PrintWriter pw = new PrintWriter(bw))
        {
            // print header
            pw.println("#chrom\tstart\tend\tscore\tclips\tmq0\tfwdHead\trevTail");
            for(Score s : this.controls){
                pw.println(s.chrom + "\t" + s.start + "\t" + s.end + "\t" + s.score + "\t" + s.clips + "\t" + s.mq0 + "\t" + s.fwdHead + "\t" + s.revTail);
            }
        }

        System.err.println("Remaining Case Number: " + this.cases.size());
        System.err.println("Remaining Control Number: " + this.controls.size());
    /* 
        // sort cases and controls by clscore
        Collections.sort(cases, new Comparator<Score>() {
            @Override
            public int compare(Score s1, Score s2) {
                if (s1.score > s2.score) return -1;
                else if (s1.score < s2.score) return 1;
                else return 0;
            }
        });
        Collections.sort(controls, new Comparator<Score>() {
            @Override
            public int compare(Score s1, Score s2) {
                if (s1.score > s2.score) return -1;
                else if (s1.score < s2.score) return 1;
                else return 0;
            }
        });*/
        BufferedWriter bw = null;
        BufferedWriter bw2 = null;
        StringWriter sw = new StringWriter();
        if (this.outputPath != null) {
            try {
                bw = new BufferedWriter(new FileWriter(this.outputPath));
                bw2 = new BufferedWriter(sw = new StringWriter());
            } catch (IOException e) {
                e.printStackTrace();
                System.exit(1);
            }
        }
        // Group scores by their value for both cases and controls
        Map<Double, List<Score>> groupedCases = groupScoresByValue(this.cases);
        Map<Double, List<Score>> groupedControls = groupScoresByValue(this.controls);

        // Extract unique scores and sort them
        List<Double> uniqueScores = new ArrayList<>(groupedCases.keySet());
        Collections.sort(uniqueScores, Collections.reverseOrder());

        // iterate for each 25 unique scores of cases
        for (int i = 0; i < uniqueScores.size(); i += 1) {
            List<Score> subcases = new ArrayList<>();
            List<Score> subcontrols = new ArrayList<>();

            double max = uniqueScores.get(i);
            double min = (i + 24 < uniqueScores.size()) ? uniqueScores.get(i + 24) : uniqueScores.get(uniqueScores.size() - 1);

            for (int j = i; j < i + 25 && j < uniqueScores.size(); j++) {
                subcases.addAll(groupedCases.get(uniqueScores.get(j)));
            }

            for (Double score : groupedControls.keySet()) {
                if (score >= min && score <= max) {
                    subcontrols.addAll(groupedControls.get(score));
                }
            }
            if(subcontrols.size() == 0) {
                System.err.println("No controls found for " + min + " to " + max);
                // continue;
            }
            // calculate the ratio of subcases and subcontrols
            double ratio = (double) subcontrols.size() / (double) subcases.size();
            double percentage = (double) subcases.size() / (double) (subcases.size() + subcontrols.size());
    
            // output the score intervals of cases and count of cases and controls and its ratio and percentage
            if (bw == null) {
                System.out.println(min + "\t" + max + "\t" + subcases.size() + "\t" + subcontrols.size() + "\t" + ratio + "\t" + percentage);
            } else {
                try {
                    bw.write(min + "\t" + max + "\t" + subcases.size() + "\t" + subcontrols.size() + "\t" + ratio + "\t" + percentage + "\n");
                    bw2.write(min + "\t" + max + "\t" + subcases.size() + "\t" + subcontrols.size() + "\t" + ratio + "\t" + percentage + "\n");
                } catch (IOException e) {
                    e.printStackTrace();
                    System.exit(1);
                }
            }
        }
        if (bw != null) {
            try {
                bw.close();
                bw2.close();
            } catch (IOException e) {
                e.printStackTrace();
                System.exit(1);
            }
        }
        if (this.outputPath != null) {
            String htmlOut = this.outputPath.replaceAll(".txt", ".html");
            try {
                BufferedWriter htmlWriter = new BufferedWriter(new FileWriter(htmlOut));
                if(separatedGraph){
                    htmlWriter.write(createSeparatedHtmlContent(sw.toString(), getFileBasename(this.outputPath)));
                }else {
                    htmlWriter.write(createHtmlContent(sw.toString(), getFileBasename(this.outputPath)));
                }
                htmlWriter.close();
            } catch (IOException e) {
                e.printStackTrace();
                System.exit(1);
            }
        }
    }
    private static String createHtmlContent(String tsvData, String sampleName) {
        String[] lines = tsvData.split("\\n");
        StringBuilder xValues = new StringBuilder();
        StringBuilder ratioYValues = new StringBuilder();
        StringBuilder percentageYValues = new StringBuilder();
    
        for (String line : lines) {
            String[] parts = line.split("\\t");
            double middleScore = (Double.parseDouble(parts[0]) + Double.parseDouble(parts[1])) / 2;
    
            xValues.append("'").append(middleScore).append("',");
            ratioYValues.append(parts[4]).append(",");
            percentageYValues.append(parts[5]).append(",");
        }
    
        return "<html>\n" +
            "<head>\n" +
            "<script src='https://cdn.plot.ly/plotly-latest.min.js'></script>\n" +
            "</head>\n" +
            "<body>\n" +
            "<h2>Graph 1: Control Count / " + sampleName + " Count</h2>\n" +
            "<div id='myDiv1'></div>\n" +
            "<script>\n" +
            "var trace1 = {\n" +
            "  x: [" + xValues + "],\n" +
            "  y: [" + ratioYValues + "],\n" +
            "  mode: 'lines+markers',\n" +
            "  type: 'scatter'\n" +
            "};\n" +
            "var layout1 = { xaxis: { range: [ 0, 50 ] }, yaxis: { range: [0, 3] }, height: 700 };\n" +
            "var data1 = [trace1];\n" +
            "Plotly.newPlot('myDiv1', data1, layout1);\n" +
            "</script>\n" +
//            "<h2>Graph 2: " + sampleName + " Count / (Control Count + " + sampleName + " Count)</h2>\n" +
//            "<div id='myDiv2'></div>\n" +
//            "<script>\n" +
//            "var trace2 = {\n" +
//            "  x: [" + xValues + "],\n" +
//            "  y: [" + percentageYValues + "],\n" +
//            "  mode: 'lines+markers',\n" +
//            "  type: 'scatter'\n" +
//            "};\n" +
//            "var layout2 = { xaxis: { range: [ 0, 50 ] }, yaxis: { range: [0, 1.1] }, height: 700 };\n" +
//            "var data2 = [trace2];\n" +
//            "Plotly.newPlot('myDiv2', data2, layout2);\n" +
//            "</script>\n" +
            "</body>\n" +
            "</html>";
    }
    
    
   private static String createSeparatedHtmlContent(String tsvData, String sampleName) {
        String[] lines = tsvData.split("\\n");
        List<String> lineList = Arrays.asList(lines);
        Collections.reverse(lineList);
        lines = lineList.toArray(new String[0]);

        List<String> labels = new ArrayList<>();
        List<String> ratioData = new ArrayList<>();
        List<String> percentageData = new ArrayList<>();
        List<String> tooltipData = new ArrayList<>();
        List<String> binCenterLabels = new ArrayList<>();
        List<String> dataPointsRatio = new ArrayList<>();
        List<String> dataPointsPercentage = new ArrayList<>();

        for (String line : lines) {
            String[] parts = line.split("\\t");

            // Calculate bin center
            String[] binBounds = { parts[0], parts[1] };
            if (binBounds.length == 1) {
                System.err.println("Error: binBounds length is 1. Check input file format." + parts[0]);
                System.exit(1);
            }
            // Calculate bin center
            double binCenter = (Double.parseDouble(binBounds[0]) + Double.parseDouble(binBounds[1])) / 2;
            // 有効桁数4桁に丸める
            binCenter = Math.round(binCenter * 10000.0) / 10000.0;

            dataPointsRatio.add("{x: " + binCenter + ", y: " + parts[4] + "},");
            dataPointsPercentage.add("{x: " + binCenter + ", y: " + parts[5] + "},");

            labels.add("'" + String.valueOf(binCenter) + "'");
            binCenterLabels.add(String.valueOf(binCenter));
            ratioData.add(parts[4]);
            percentageData.add(parts[5]);
            tooltipData.add("{label: 'Subcases: " + parts[2] + ", Subcontrols: " + parts[3] + "'},");
        }

        String dataPointsRatioStr = String.join(",", dataPointsRatio);
        String dataPointsPercentageStr = String.join(",", dataPointsPercentage);

        Collections.reverse(labels);
        Collections.reverse(binCenterLabels);

        String labelsStr = String.join(",", labels);
        String ratioDataStr = String.join(",", ratioData);
        String percentageDataStr = String.join(",", percentageData);
        String tooltipDataStr = String.join(",", tooltipData);
        String binCenterLabelsStr = String.join(",", binCenterLabels);

        return "<!DOCTYPE html>\n" +
                "<html>\n" +
                "<head>\n" +
                "<script src='https://cdn.jsdelivr.net/npm/chart.js'></script>\n" +
                "</head>\n" +
                "<body>\n" +
                "<canvas id='ratioChart'></canvas>\n" +
                "<canvas id='percentageChart'></canvas>\n" +
                "<canvas id='ratioBinCenterChart'></canvas>\n" +
                "<canvas id='percentageBinCenterChart'></canvas>\n" +
                "<script>\n" +
                "function createChartConfig(title, labels, data) {" +
                "    return {" +
                "        type: 'line'," +
                "        data: {" +
                "            labels: labels," +
                "            datasets: [{" +
                "                label: title," +
                "                data: data," +
                "                borderColor: 'rgba(75, 192, 192)'," +
                "                backgroundColor: 'rgba(75, 192, 192)'," +
                "            }]" +
                "        }," +
                "        options: {" +
                "            responsive: true," +
                "            title: {" +
                "                display: true," +
                "                text: title" +
                "            }," +
                "            scales: {" +
                "                y: {" +
                "                    beginAtZero: true" +
                "                }" +
                "            }" +
                "       }" +
                "    };" +
                "}" +

                "var tooltipData = [" + tooltipDataStr + "];\n" +
                "var labels = [" + labelsStr + "];\n" +
                "var binCenterLabels = [" + binCenterLabelsStr + "];\n" +
                "var ctxRatio = document.getElementById('ratioChart').getContext('2d');\n" +
                "var ctxPercentage = document.getElementById('percentageChart').getContext('2d');\n" +
                "var ctxRatioBinCenter = document.getElementById('ratioBinCenterChart').getContext('2d');\n" +
                "var ctxPercentageBinCenter = document.getElementById('percentageBinCenterChart').getContext('2d');\n" +
                "new Chart(ctxRatio, createChartConfig('Ratio (control / (" + sampleName + ")', labels, [" + dataPointsRatioStr + "]));\n" +
                "new Chart(ctxPercentage, createChartConfig('Percentage (" + sampleName + " / (Case + Control))', labels, [" + dataPointsPercentageStr + "]));\n" +
                "new Chart(ctxRatioBinCenter, createChartConfig('Ratio (control / (" + sampleName + ")', binCenterLabels, [" + ratioDataStr + "]));\n" +
                "new Chart(ctxPercentageBinCenter, createChartConfig('Percentage (" + sampleName + " / (Case + Control))', binCenterLabels, [" + percentageDataStr + "]));\n" +
                "</script>\n" +
                "</body>\n" +
                "</html>";
}


    public static void main(String[] args){
        try {
            String outputPath = null;
            boolean debug = false;
            if(args.length < 2){
                System.err.println("Usage: java ScoreCheck [--mq0 <mq0>] [--clips <clips>] [--debug] <case.bed> <control.bed>");
                System.err.println("     --threshold: filter for low quality scores (default: 7.0)");
                System.exit(1);
            }
            for(int i = 0; i<args.length-1; i++){
                if(args[i].equals("--out")){
                    outputPath = args[i+1];
                    i++;
                }
                if(args[i].equals("--threshold")){
                    scoreThreshold = Double.parseDouble(args[i+1]);
                    i++;
                }
                if(args[i].equals("--debug")){
                    debug = true;
                }
            }
            ScoreCheck sc = new ScoreCheck(args[args.length-2], args[args.length-1], outputPath);
            System.err.println("Effective Case Number: " + sc.cases.size());
            System.err.println("Effective Control Number: " + sc.controls.size());
            if(debug){
                // dump scores to log
                System.err.println("Cases:");
                for(Score s : sc.cases){
                    System.out.println("CASE:" + "\t" + s.chrom + "\t" + s.start + "\t" + s.end + "\t" + s.score);
                }
                System.err.println("Controls:");
                for(Score s : sc.controls){
                    System.out.println("CONTROL:" + "\t" + s.chrom + "\t" + s.start + "\t" + s.end + "\t" + s.score);
                }
            }
            sc.run();

        }catch (IOException e){
            e.printStackTrace();
        }
    }

}
