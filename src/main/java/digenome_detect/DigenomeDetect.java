package digenome_detect;

import java.io.*;
import java.lang.Math;
import java.math.BigInteger;
import java.util.*;
import java.lang.AutoCloseable;
// import java.util.zip.GZIPInputStream;
import org.apache.commons.math3.distribution.ExponentialDistribution;
import org.broadinstitute.hellbender.utils.*;

public class DigenomeDetect implements AutoCloseable{
    public static final double LOG_2 = Math.log(2.0);
    public static final double LOG_10 = Math.log(10.0);
    public boolean is_siteseq = false;
    public static boolean calc_cs = false;
    public static boolean calc_fisher = true;
    public boolean inplace_depth = false;
    public boolean inplace_depth2 = true;

    // numbers greater than 10^MAX_DIGITS_10 or e^MAX_DIGITS_E are considered unsafe ('too big') for floating point operations
    private static final int MAX_DIGITS_10 = 294;
    private static final int MAX_DIGITS_2 = 977; // ~ MAX_DIGITS_10 * LN(10)/LN(2)
    private static final int MAX_DIGITS_E = 677; // ~ MAX_DIGITS_10 * LN(10)

    static int DETECT_WIDTH = 1;
    public final static int READ_LENGTH = 150;
    int[] block_pos = new int[READ_LENGTH*4];
    int[] block_depth = new int[READ_LENGTH*4];
    int[] block_fwd_heads = new int[READ_LENGTH*4];
    int[] block_fwd_tails = new int[READ_LENGTH*4];
    int[] block_rev_heads = new int[READ_LENGTH*4];
    int[] block_rev_tails = new int[READ_LENGTH*4];
    int[] block_mq0 = new int[READ_LENGTH*4];
    int[] block_softclips = new int[READ_LENGTH*4];
    int[] block_body = new int[READ_LENGTH*4];
    int[] block_fwd_depth = new int[READ_LENGTH*4];
    int[] block_rev_depth = new int[READ_LENGTH*4];
    boolean debug = false;
    PrintWriter bed = null;
    ControlChecker checker = null;
    public DigenomeDetect(int detectionWidth, OutputStream ostream) throws IOException{
        DETECT_WIDTH = detectionWidth;
        bed = new PrintWriter(new BufferedWriter(new OutputStreamWriter(ostream)));
    }
    public DigenomeDetect(int detectionWidth) {
        DETECT_WIDTH = detectionWidth;
        bed = new PrintWriter(new BufferedWriter(new OutputStreamWriter(System.out)));
    }
    public void setControlBam(String bamPath){
        if(bamPath != null){
            checker = new ControlChecker(bamPath);
        }
    }
    public void setInplaceDepth(boolean b){
        this.inplace_depth = b;
    }
    public void setInplaceDepth2(boolean b){
        this.inplace_depth2 = b;
    }
    public void setSiteSeq(boolean b){
        this.is_siteseq = b;
    }
    public void setCalcCS(boolean b){
        this.calc_cs = b;
    }
    public void setCalcFisher(boolean b){
        this.calc_fisher = b;
    }
    public void setDebug(boolean b){
        this.debug = b;
    }
    public void close() throws Exception {
        if(cluster.size() > 0){
            handle(cluster);
        }
        bed.close();
    }
    int get_center_genomic(String src){
        // println!("---- len={}, at {}"
        String[] line = src.split(" ");
        if(!line[2].equals("at")){
            System.err.println("bad line:  " +  src);
            System.exit(-1);
        }
        return Integer.parseInt(line[3]);
    }
    int lastPos = 0;
    ArrayList<ArrayList<String>> cluster  = new ArrayList<ArrayList<String>>();
    public void push(ArrayList<String> list){
        // System.out.println("#####" + list.get(0));
        int pos = get_center_genomic(list.get(0));
        if(pos - lastPos < 10 || cluster.size() == 0){
            cluster.add(list);
        }else if(cluster.size() > 0){
            handle(cluster);
            cluster.clear();
            cluster.add(list);
        }
        lastPos = pos;
    }
    private void handle(ArrayList<ArrayList<String>> buffer){
        ArrayList<Result> results = new ArrayList<Result>();
        for(ArrayList<String> block: buffer){
            try{
                ArrayList<Result> b = analyze(block);
                if(b != null){
                    results.addAll(b);
                }
            }catch(ArrayIndexOutOfBoundsException e){
                // near breakpoint
            }
        }
        if(results.size() > 0){
            Result best = chooseBest(results);
            if( best.table[1][1] > 2 && best.table[0][0] > 2){
                double cscore = (best.cscore == null)? 0.0: best.cscore.getScore();
                // bed.printf(best.chr + "\t" + best.start + "\t" + best.end
                //    + "\tCLSCORE="+ format("%.2f", best.phred)+";DP="+ best.median + ";CS=%.2f;"
                //    + "Ratio="+ format("%.3f", best.ratio)+";FISHER="+format("%.3f", best.fisher)+";RevHead="+ best.table[0][1] +";RevTail="+ best.table[1][1]
                //    + ";FwdHead=" + best.table[0][0] + ";FwdTail="+ best.table[1][0]+";MQ0=%d;CLIPS=%d", cscore, best.mq0, best.clips);
                bed.printf(best.chr + "\t" + best.start + "\t" + best.end + "\t");
                bed.printf("CLSCORE=%.2f;", best.phred);
                bed.printf("DP=%.1f;", best.median);
                if(calc_cs){
                    bed.printf("CS=%.2f;", cscore);
                }
                bed.printf("CS=%.2f;", cscore);
                bed.printf("Ratio=%.3f;", best.ratio);
                if(calc_fisher){
                    bed.printf("FISHER=%.4f;", best.fisher);
                }
                bed.printf("RevHead=" + best.table[0][1] + ";");
                bed.printf("RevTail=" + best.table[1][1] + ";");
                bed.printf("FwdHead=" + best.table[0][0] + ";");
                bed.printf("FwdTail=" + best.table[1][0] + ";");
                bed.printf("MQ0=%d;", best.mq0);
                bed.printf("CLIPS=%d;", best.clips);
                
                if(checker != null){
                    int[] control = checker.check(best.chr, best.start, best.end);
                    int cont_depth_start = checker.getDepth(best.chr, best.start-1);
                    int cont_depth_end = checker.getDepth(best.chr, best.end+1);
                    double score = 0.0;
                    if(inplace_depth2){
                        score = calcProb_with_inplaceDepth2(
                            cont_depth_start, cont_depth_end,
                            control[0],
                            control[3],
                            best.end-best.start+1);
                    }else if(inplace_depth){
                        score = calcProb_with_inplaceDepth(
                            cont_depth_start, cont_depth_end,
                            control[0],
                            control[3],
                            best.end-best.start+1);
                    }else {
                        score = calcProb(
                            (cont_depth_start+cont_depth_end)/2,
                            control[0],
                            control[3],
                            best.end-best.start+1);
                    }

                    bed.print(";CONTROL=" + control[0] + "," + control[1] + "," + control[2] + ","+ control[3] +";CONTROL_CLSCORE="+format("%.2f", score));
                }
                bed.print("\n");
                bed.flush();
            }
        }
    }
    public ArrayList<Result> analyze(ArrayList<String> block){
        String chr = null;
        if(block.size() < 100){
            // System.out.println("ok");
            System.err.println("block size was too small: " + block.size() + " at " + get_center_genomic(block.get(0)));
            return null;
        }
        int center_genomic = get_center_genomic(block.get(0));
        int center = -1;
        int block_size = 0;
        for(int i = 1; i<block.size(); i++){
            // System.out.println(block.get(i));
            String[] line = block.get(i).split("\t");
            chr = line[0];
            block_pos[i-1] = Integer.parseInt(line[1]);
            if(block_pos[i-1] == center_genomic){
                center = i-1;
            }
            block_depth[i-1] = Integer.parseInt(line[2]);
            block_fwd_heads[i-1] = Integer.parseInt(line[3]);
            block_fwd_tails[i-1] = Integer.parseInt(line[4]);
            block_rev_heads[i-1] = Integer.parseInt(line[5]);
            block_rev_tails[i-1] = Integer.parseInt(line[6]);
            block_mq0[i-1] = Integer.parseInt(line[7]);
            block_softclips[i-1] = Integer.parseInt(line[8]);
            block_fwd_depth[i-1] = Integer.parseInt(line[9]);
            block_rev_depth[i-1] = Integer.parseInt(line[10]);
            // block_body[i-1] = Integer.parseInt(line[8]);
            block_size++;
        }
        double[] median_mean = calcMedianAndMean(block_depth, block_size, center);
        // check
        boolean skip = false;
        if(is_siteseq){
           skip = false;
           inplace_depth = true;
        }else if(median_mean[0] > 500 || median_mean[1] > 500){ // on repeat
            skip = true;
        }
  //    for(int i = center-3; i<center+3; i++){
  //        if(block_rev_tails[i-1] < block_rev_heads[i-1] || block_rev_tails[i-1] < block_fwd_tails[i]
  //            || block_fwd_heads[i] < block_rev_heads[i-1] || block_fwd_heads[i] < block_fwd_tails[i]){
  //                skip = true;
  //    //  }else if(block_rev_tails[i-1] > 1 && block_fwd_heads[i] > 1 && block_pos[i]-1 == block_pos[i-1]){
  //    //      skip = false;
  //            break;
  //        }
  //    }
        if(debug && !skip){
            System.err.println("-----------------------------------------------------------------------------------");
            System.err.println(block.get(0));
            for(int p = center-5; p<center+5; p++){
                if(p == center){
                    System.err.print("*");
                }
                System.err.println(chr + "\t"+ block_pos[p] + "\tFH:" + block_fwd_heads[p] +"\tRT:" + block_rev_tails[p]
                    + "\tFT:" + block_fwd_tails[p] + "\tRH:" + block_rev_heads[p] + "\tMQ0:" + block_mq0[p] + "\tCLIP:" + block_softclips[p]);
            }
        }
        if(!skip){
            CleavageScore cscore = null;
            if(!is_siteseq && calc_cs){
               cscore =  calcCleavageScore(block_pos, block_depth, block_fwd_heads, block_rev_tails, center);
            }
            ArrayList<Result> results = new ArrayList<Result>();
            for(int width = 0; width < DETECT_WIDTH; width++){
                // calc ratio
                for(int frame = 0; frame < width + 1; frame++){
                    Result result = new Result();
                    result.median  = median_mean[0];
                    result.chr = chr;
                    result.center = center;
                    result.width = width;
                    result.cscore = cscore;
                    result.ratio = calcRatio(center, block_rev_tails, block_fwd_heads, block_fwd_tails, block_rev_heads, width);
                    int rt_start = center-1-width+frame;
                    int rt_end =   center-1-width+frame+width + 1;
                    int fh_start = center - width + frame;
                    int fh_end =   center - width + frame + width + 1;
                    if(inplace_depth){
                        result.inplace = calcProb_with_inplaceDepth(
                            block_depth[center+10], block_depth[center-10],
                            sum(block_fwd_heads, fh_start, fh_end),
                            sum(block_rev_tails, rt_start, rt_end),
                            width);
                        result.phred = result.inplace;
                    }else if(inplace_depth2){
                        result.inplace2 = calcProb_with_inplaceDepth2(
                            block_fwd_depth[center+10], block_rev_depth[center-10],
                            sum(block_fwd_heads, fh_start, fh_end),
                            sum(block_rev_tails, rt_start, rt_end),
                            width);
                        result.phred = result.inplace2;
                    }else {
                        result.phred = calcProb(
                            (int)median_mean[0],
                            sum(block_fwd_heads, fh_start, fh_end),
                            sum(block_rev_tails, rt_start, rt_end),
                            width);
                    }
                    /*
                    if(result.phred > 10000.0 || result.phred < 0){
                        result.phred = 10000.0;
                    }*/
                    result.table[0][0] = sum(block_fwd_heads, fh_start, fh_end);
                    result.table[0][1] = sum(block_rev_heads, fh_start, fh_end);
                    result.table[1][0] = sum(block_fwd_tails, rt_start, rt_end);
                    result.table[1][1] = sum(block_rev_tails, rt_start, rt_end);
                    if(calc_fisher){
                        result.fisher = Math.abs(-10*Math.log10(Fisher.exactTest(result.table)));
                    }
                    result.start = block_pos[fh_start];
                    result.end = block_pos[fh_start+width];
                    result.mq0 =  sum(block_mq0, fh_start-2, fh_end+2);
                    result.clips = sum(block_softclips, fh_start, fh_end);
                    // System.err.println("clips = " + result.clips + " from " + fh_start + " to " + (fh_end));
                    /*
                    if(result.fisher > 1200){
                        result.fisher = 1200;
                    }*/
                    if(debug){
                        System.err.println("fh_start: " + fh_start + ", fh_end: " + fh_end);
                        System.err.println("rt_start: " + rt_start + ", rt_end: " + rt_end);
                        System.err.println("interval: " + block_pos[fh_start] + ".." + block_pos[fh_end]);
                        System.err.println("interval: " + block_pos[rt_start] + ".." + block_pos[rt_end]);
                        System.err.println("phred = " + result.phred
                            + "[" + sum(block_fwd_heads, fh_start, fh_end) + ", "
                            + sum(block_rev_tails, rt_start, rt_end) + "]");
                        System.err.println("fisher = "  + result.fisher);
                        System.err.println("width = " + width + " frame = " + frame + " DETECT_WIDTH: "+ DETECT_WIDTH);
                        System.err.println("for inplace-score:");
                        System.err.println("    " + block_depth[fh_start] + ", " +  block_depth[fh_end] + ", " +
                            sum(block_fwd_heads, fh_start, fh_end) + ", " +
                             sum(block_rev_tails,  rt_start,  rt_end));

                    }
                    results.add(result);
                }
            }
            return results;
        }
        return null;
    }
    private Result chooseBest(ArrayList<Result> candidates){
        // without fisher
        final Comparator<Result> phredComp = new Comparator<Result>(){
            public int compare(Result r1, Result r2){
                if(r1.phred > r2.phred){
                    return -1;
                }else if(r1.phred < r2.phred){
                    return 1;
                }
                return 0;
            }
            public boolean equals(Comparator<Result> r){
                return true;
            }
        };
        // by fisher
//        final Comparator fisherComp = new Comparator<Result>(){
//            public int compare(Result r1, Result r2){
//                if(r1.fisher > r2.fisher){
//                    return -1;
//                }else if(r1.fisher < r2.fisher){
//                    return 1;
//                }else if(r1.phred > r2.phred){
//                    return -1;
//                }else if(r1.phred < r2.phred){
//                    return 1;
//               }
//                return 0;
//            }
//            public boolean equals(Comparator<Result> r){
//                return true;
//            }
//        };

        Result best = candidates.get(0);
        Comparator comp = phredComp;
        for(Result r: candidates){
            if(comp.compare(best, r) > 0){
                best = r;
            }
        }
        return best;
    }
    public class Result {
        String chr;
        int center;
        int width;
        int start;
        int end;
        int mq0;
        int clips;
        CleavageScore cscore;
        int[][] table;
        double phred;
        double background;
        double inplace;
        double inplace2;
        double fisher;
        double ratio;
        double median;
        public Result(){
            table = new int[2][]; // for fisher exact test
            table[0] = new int[2];
            table[1] = new int[2];
        }
    }
    String format(String fmt, double val){
        return String.format(fmt, val);
    }
    String join(double[] buf, String fmt){
        StringBuilder str = new StringBuilder();
        str.append(String.format(fmt, buf[0]));
        for(int i = 1; i<buf.length; i++){
            str.append(",");
            str.append(String.format(fmt, buf[i]));
        }
        return str.toString();
    }
    String join(ArrayList<int[][]> tables, int row, int col){
        StringBuilder buf = new StringBuilder();
        buf.append(tables.get(0)[row][col]);
        for(int i = 1; i<tables.size(); i++){
            buf.append(",");
            buf.append(String.valueOf(tables.get(i)[row][col]));
        }
        return buf.toString();
    }
    public double calcRatio(int idx, int[] rev_tails, int[] fwd_heads, int[] fwd_tails, int[] rev_heads, int width){
        int start = idx - 1 - width;
        int end = idx + width + 1;
        double ratio = (double)(sum(block_rev_tails, start, end) + sum(block_fwd_heads, start, end))
                        /
                        (double)(sum(block_rev_tails, start, end) + sum(block_fwd_heads, start, end) + sum(block_rev_heads, start, end) + sum(block_fwd_tails, start, end));
        return ratio;
    }
    int sum(int[] buf, int start, int end){
        int ret = 0;
        for(int i = start; i<end; i++){
            ret += buf[i];
        }
        return ret;
    }
    public static BigInteger factorialHavingLargeResult(int n) {
        BigInteger result = BigInteger.ONE;
        for (int i = 2; i <= n; i++)
            result = result.multiply(BigInteger.valueOf(i));
        return result;
    }
    public static double logBigInteger(BigInteger val) {
        if (val.signum() < 1)
            return val.signum() < 0 ? Double.NaN : Double.NEGATIVE_INFINITY;
        int blex = val.bitLength() - MAX_DIGITS_2; // any value in 60..1023 works here
        if (blex > 0)
            val = val.shiftRight(blex);
        double res = Math.log(val.doubleValue());
        return blex > 0 ? res + blex * LOG_2 : res;
    }
    public double calcProb_with_inplaceDepth2(int fdepth, int rdepth, int fcount, int rcount, int width){
        double fwd_p = 0;
        double rev_p = 0;
        double flambda = ((double)fdepth)/READ_LENGTH * (width + 1);
        double rlambda = ((double)rdepth)/READ_LENGTH * (width + 1);
        try {
            double denomitor = logBigInteger(factorialHavingLargeResult(fcount));
            double numerator = fcount * Math.log10(flambda) - (flambda * Math.log10(Math.E));
            fwd_p = numerator - denomitor/Math.log(10);

            denomitor = logBigInteger(factorialHavingLargeResult(rcount));
            numerator = rcount * Math.log10(rlambda) - (rlambda * Math.log10(Math.E));
            rev_p = numerator - denomitor/Math.log(10);
        }catch(Exception e){
            System.err.println(e.getMessage());
            e.printStackTrace();
            return 0;
        }
        return -1*(fwd_p + rev_p);
    }

    public double calcProb_with_inplaceDepth(int fdepth, int rdepth, int fcount, int rcount, int width){
        double fwd_p = 0;
        double rev_p = 0;
        double flambda = ((double)fdepth)/READ_LENGTH * (width + 1);
        double rlambda = ((double)rdepth)/READ_LENGTH * (width + 1);
        try {
            double denomitor = logBigInteger(factorialHavingLargeResult(fcount));
            double numerator = fcount * Math.log10(flambda) - (flambda * Math.log10(Math.E));
            // double numerator = Math.log10(Math.pow(flambda, fcount) * Math.exp(-flambda));
            fwd_p = numerator - denomitor/Math.log(10);

            denomitor = logBigInteger(factorialHavingLargeResult(rcount));
            numerator = rcount * Math.log10(rlambda) - (rlambda * Math.log10(Math.E));
            // numerator = Math.log10(Math.pow(rlambda, rcount) * Math.exp(-rlambda));
            rev_p = numerator - denomitor/Math.log(10);
        }catch(Exception e){
            // System.err.println("forward depth=" + fdepth);
            // System.err.println("reverse depth=" + rdepth);
            System.err.println(e.getMessage());
            e.printStackTrace();
            return 0;
        }
        return -1*(fwd_p + rev_p);
    }
    public double calcProb(int depth, int fcount, int rcount, int width){
        double fwd_p = 0;
        double rev_p = 0;
        double lambda = ((double)depth)/READ_LENGTH * (width + 1);
        try {

            double denomitor = logBigInteger(factorialHavingLargeResult(fcount));
            double numerator = Math.log10(Math.pow(lambda, fcount) * Math.exp(-lambda));
            fwd_p = numerator - denomitor/Math.log(10);

            denomitor = logBigInteger(factorialHavingLargeResult(rcount));
            numerator = Math.log10(Math.pow(lambda, rcount) * Math.exp(-lambda));
            rev_p = numerator - denomitor/Math.log(10);
        }catch(Exception e){
            System.err.println("depth=" +depth);
            e.printStackTrace();
            return 0;
        }
        return -1*(fwd_p + rev_p);
    }
    public double calcEfficiency(int f, int r, int b){
        double sum = (f + r)/2.0 + b;
        return (f + r)/2.0/sum;
    }
    public CleavageScore calcCleavageScore(int[] block_gpos, int[] block_depth, int[] block_fwd_heads, int[] blcok_rev_tails, int center){
        int[] D = block_depth;
        int[] F = block_fwd_heads;
        int[] R = block_rev_tails;
        int maxCSpos = 0;
        double max_cleavage_score = 0;
        for(int i = center-10; i < center+10; i++){
            // score at position i
            double cleavage_score = 0.0;
            double score_i = 0.0;
            // forward
            try {
                for(int a = 1; a<=5; a++){
                    score_i += (double)R[i-4+a]/D[i-4+a] * (F[i] + R[i-4+a] - 2);
                }
                cleavage_score += score_i*(F[i]-1)/D[i];
                for(int a = 1; a<=5; a++){
                    score_i += (double)(F[i-3+a]-1)/D[i-3+a] * (R[i-1]+F[i-3+a]-2);
                }
                cleavage_score += score_i * (R[i-1]-1)/D[i-1];
                if(cleavage_score > max_cleavage_score){
                    max_cleavage_score = cleavage_score;
                    maxCSpos = block_gpos[i];
                }
            }catch(ArrayIndexOutOfBoundsException aiobe){
                // caused by neighbouring gap
            }
        }
        return new CleavageScore(maxCSpos, max_cleavage_score);
    }
    public int mean(int[] buf){
        int sum = 0;
        for(int i: buf){
            sum += i;
        }
        return sum/buf.length;
    }
    public double[] calcMedianAndMean(int[] depth, int size, int center){
        double[] result = new double[2];
        int[] _depth = new int[READ_LENGTH*2-20];
        int didx = 0;
        int sum = 0;
        for(int i = 0; didx<READ_LENGTH*2-20 && i<size; i++){
            if(i<center-READ_LENGTH-10 || i >= center + READ_LENGTH + 10){ // avoid ends of stacks of reads
                _depth[didx] = depth[i];
                didx++;
                sum += depth[i];
            }
        }
        int[] pdepth = Arrays.copyOf(_depth, didx);
        Arrays.sort(pdepth);
        double mean = sum/pdepth.length;
        double median = pdepth[pdepth.length/2]; // slightly above of
        result[0] = median;
        result[1] = mean;

        return result;
    }
    public static String join(String[] buf){
        StringBuilder b = new StringBuilder();
        b.append(buf[0]);
        for(int i = 1; i<buf.length; i++){
            b.append("\t");
            b.append(buf[i]);
        }
        return b.toString();
    }
}
