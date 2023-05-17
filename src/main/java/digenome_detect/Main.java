package digenome_detect;

// import org.apache.logging.log4j.core.*;
import org.slf4j.LoggerFactory;
import org.slf4j.Logger;

import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.FileOutputStream;
import java.io.BufferedReader;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

public class Main {
    static String out = null;
    static String bamPath = null;
    static String digenomeDetectPath = "/ldisk1/202205_digenome_mod/digenome-detect_mq10_2depth/rust/target/debug/digenome_seek";
    static int threads = 8;
    static String mqfilter = "0";
    static String[] regions = null;
    static boolean debug = false;
    static int detectWidth = 3;
    static boolean is_siteseq = false;
    static boolean inplaceDepth = false;
    static boolean inplaceDepth2 = false;

    private static final Logger LOGGER = LoggerFactory.getLogger(Main.class);

    public static void main(String[] argv){
        LOGGER.info("Running the main method");
        try {
            for(int i = 0; i<argv.length; i++){
                if(argv[i].equals("--bam")){
                    bamPath = argv[i+1];
                }else if(argv[i].equals("--out")){
                    out = argv[i+1];
                }else if(argv[i].equals("--inplace-depth")){
                    inplaceDepth = true;
                }else if(argv[i].equals("--inplace-depth2")){
                    inplaceDepth2 = true;
                }else if(argv[i].equals("--threads") || argv[i].equals("-t")){
                    threads = Integer.parseInt(argv[i+1]);
                }else if(argv[i].equals("--regions") || argv[i].equals("--region")){
                    if(argv[i+1].startsWith("GRCm")){
                        regions = mouse_chromosomes;
                    }else {
                        regions = argv[i+1].split(",");
                    }
                }else if(argv[i].equals("--debug")){
                    debug = true;
                }else if(argv[i].equals("--siteseq")){
                    is_siteseq = true;
                }else if(argv[i].equals("--mq")){
                    mqfilter = argv[i+1];
                }else if(argv[i].equals("--width")){
                    detectWidth = Integer.parseInt(argv[i+1]);
                }else if(argv[i].equals("--help")){
                    System.out.println("Usage:");
                    System.out.println("     java digenome_detect.Main --out [output prefix] --bam [bam file]");
                    System.out.println("Options:");
                    System.out.println("   -t or --threads: number of threads");
                    System.out.println("   --region       : target regions like 'chr19:12345..23456'");
                    System.exit(0);
                }
            }
            if(regions == null){
                regions = chromosomes;
            }
            // check consistency
            if(bamPath == null){
                System.err.println("Use --bam option for input BAM file");
                System.exit(-1);
            }else {
                checkBamFile(bamPath);
            }
            if(out == null){
                System.err.println("Use --out option for output prefix");
                System.exit(-1);
            }
            run(bamPath, out);
        }catch(Exception e){
            System.err.println(e.getMessage());
            System.exit(-1);
        }
    }
    private static void run(String bamPath, String out) throws IOException, InterruptedException{
        ExecutorService service = Executors.newFixedThreadPool(threads);
        CountDownLatch latch = new CountDownLatch(regions.length);
        List<DetectTask> tasks = new ArrayList<DetectTask>();
        for ( int i = 0; i<regions.length; i++){
            tasks.add(new DetectTask(bamPath, out, regions[i], latch));
        }
        for(DetectTask t: tasks){
            service.submit(t);
        }
        latch.await();
        service.shutdown();
    }
    public static class DetectTask  implements Runnable {
        private String bamPath = null;
        private String out = null;
        private CountDownLatch latch;
        private String chr;
        private DigenomeDetect detect = null;
        public DetectTask(String bamPath_, String out_, String chr_, CountDownLatch latch_) throws IOException{
            bamPath = bamPath_;
            out = out_;
            chr = chr_;
            latch = latch_;
            detect = new DigenomeDetect(detectWidth, new FileOutputStream(out + "." + chr + ".bed"));

            detect.setInplaceDepth(inplaceDepth);
            detect.setInplaceDepth2(inplaceDepth2);
            detect.setSiteSeq(is_siteseq);
            detect.setDebug(debug);
        }

        @Override
        public void run(){
          try {
            // do
            Runtime rt = Runtime.getRuntime();
            // String[] cmd = {digenomeDetectPath, chr, bamPath};
            // Process proc = rt.exec(cmd);
            ProcessBuilder builder = new ProcessBuilder(digenomeDetectPath, chr, bamPath, mqfilter);
            Process proc = builder.start();
            Thread stdThread = new Thread(){
                public void run(){
                  try {
                    BufferedReader br = new BufferedReader(new InputStreamReader(proc.getInputStream()));
                    ArrayList<String> buf = new ArrayList<String>(DigenomeDetect.READ_LENGTH*4+10);
                    String line = null;
                    while(null != (line = br.readLine())){
                        if(line.startsWith("//")){
                            if(buf.size() > DigenomeDetect.READ_LENGTH*3){
                                detect.push(buf);
                            }
                            buf = new ArrayList<String>(DigenomeDetect.READ_LENGTH*4+10);
                        }else {
                            buf.add(line);
                        }
                    }
                    detect.close();
                  }catch(IOException e){
                    e.printStackTrace();
                    throw new RuntimeException(e.getMessage());
                  }
                }
            };
            stdThread.start();
            Thread errThread = new Thread(){
              public void run(){
                try {
                    BufferedReader b = new BufferedReader(new InputStreamReader(proc.getErrorStream()));
                        String eline = null;
                        while(null != (eline = b.readLine())){
                            System.err.println(eline);
                        }
                        b.close();
                    }catch(IOException e){
                        e.printStackTrace();
                    }
                }
            };
            errThread.start();

            int exitCode = proc.waitFor();
            if(exitCode != 0){
                System.err.println("digenome_seek exited with code = " + exitCode + " on " + chr);
                System.exit(-1);
            }
            // finish
            latch.countDown();
            System.err.println("count: " + latch.getCount());
          }catch(IOException ioe){
              throw new RuntimeException(ioe);
          }catch(InterruptedException inte){
              throw new RuntimeException(inte);
          }
        }
    }
    private static void checkBamFile(String path){
        java.io.File file = new java.io.File(path);
        if(!file.exists()){
            throw new RuntimeException("no such file: " + path);
        }
    }
    public static final String[] chromosomes = {
        "chr1",
        "chr2",
        "chr3",
        "chr4",
        "chr5",
        "chr6",
        "chr7",
        "chr8",
        "chr9",
        "chr10",
        "chr11",
        "chr12",
        "chr13",
        "chr14",
        "chr15",
        "chr16",
        "chr17",
        "chr18",
        "chr19",
        "chr20",
        "chr21",
        "chr22",
        "chrX",
        "chrY"
    };
    public static final String[] mouse_chromosomes = {
        "chr1",
        "chr2",
        "chr3",
        "chr4",
        "chr5",
        "chr6",
        "chr7",
        "chr8",
        "chr9",
        "chr10",
        "chr11",
        "chr12",
        "chr13",
        "chr14",
        "chr15",
        "chr16",
        "chr17",
        "chr18",
        "chr19",
        "chrX",
        "chrY"
    };
}
