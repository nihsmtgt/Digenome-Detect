package digenome_detect;

// import org.apache.logging.log4j.core.*;
import org.slf4j.LoggerFactory;
import org.slf4j.Logger;

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
    static String controlPath = null;
    static String digenomeDetectPath = "/usr/local/bin/digenome_seek";
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
                }else if(argv[i].equals("--control")){
                    controlPath = argv[i+1];
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
                }else if(argv[i].equals("--strandbias")){
                    if(argv[i+1].equals("true")){
                        DigenomeDetect.calc_fisher = true;
                    }else {
                        DigenomeDetect.calc_fisher = false;
                    }
                }else if(argv[i].equals("--calc_cleavage_score")){
                    if(argv[i+1].equals("true")){
                        DigenomeDetect.calc_cs = true;
                    }else {
                        DigenomeDetect.calc_cs = false;
                    }
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
                }else if(argv[i].startsWith("--")){
                    System.err.println("Unknown option: " + argv[i]);
                    System.exit(-1);
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
            run(bamPath, controlPath, out);
        }catch(Exception e){
            System.err.println(e.getMessage());
            e.printStackTrace();
            System.exit(-1);
        }
    }
    private static void run(String bamPath, String controlPath, String out) throws IOException, InterruptedException{
        ExecutorService service = Executors.newFixedThreadPool(threads);
        CountDownLatch latch = new CountDownLatch(regions.length);
        List<DetectTask> tasks = new ArrayList<DetectTask>();
        System.err.println("running " + regions.length + " regions");
        for ( int i = 0; i<regions.length; i++){
            System.err.println("bam: " + bamPath + " control: " + controlPath + " out: " + out + " chr: " + regions[i]);
            tasks.add(new DetectTask(bamPath, controlPath, out, regions[i], latch));
        }
        for(DetectTask t: tasks){
            service.submit(t);
        }
        // wait for all tasks to finish
        latch.await();        
        service.shutdown();
    }
    public static class DetectTask  implements Runnable {
        private String bamPath = null;
        private String controlPath = null;
        private String out = null;
        private CountDownLatch latch;
        private String chr;
        private DigenomeDetect detect = null;
        public DetectTask(String bamPath_, String controlPath_, String out_, String chr_, CountDownLatch latch_) throws IOException{
            bamPath = bamPath_;
            controlPath = controlPath_;
            out = out_;
            chr = chr_;
            latch = latch_;
            detect = new DigenomeDetect(detectWidth, new FileOutputStream(out + "." + chr + ".bed"));
            detect.setControlBam(controlPath);

            detect.setInplaceDepth(inplaceDepth);
            detect.setInplaceDepth2(inplaceDepth2);
            detect.setSiteSeq(is_siteseq);
            detect.setDebug(debug);
        }
        public class AutoCloseableThread extends Thread implements AutoCloseable {
            Process proc = null;
            public AutoCloseableThread(Process p){
                proc = p;
            }
            public void run(){
              try (BufferedReader br = new BufferedReader(new InputStreamReader(proc.getInputStream()))){
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
              }catch(IOException e){
                  System.err.println(e.getMessage());
                  e.printStackTrace();
              }
            }
            @Override
            public void close() throws Exception {
                if(detect != null){
                    detect.close();
                }
            }
        }

        @Override
        public void run(){
          try {
            // do
            Runtime rt = Runtime.getRuntime();
            ProcessBuilder builder = new ProcessBuilder(digenomeDetectPath, chr, bamPath, mqfilter);
            Process proc = builder.start();
            AutoCloseableThread stdThread = new AutoCloseableThread(proc);
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
