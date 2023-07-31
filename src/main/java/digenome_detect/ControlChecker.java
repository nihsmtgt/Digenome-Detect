package digenome_detect;

import htsjdk.samtools.*;
import java.io.IOException;
import java.io.File;
import htsjdk.samtools.filter.*;

public class ControlChecker implements AutoCloseable {
    SamReader reader = null;
    public static void main(String[] argv) {
        // 入力ファイルのパス
        String inputFilePath = argv[1];
        ControlChecker checker = new ControlChecker(inputFilePath);

        // 分類する範囲の開始位置と終了位置
        String[] regionParts = argv[0].split(":");
        String contig = regionParts[0];
        String[] positions = regionParts[1].split("-");
        int start = Integer.parseInt(positions[0]) - 1;
        int end = Integer.parseInt(positions[1]) + 1;

        int[] result = checker.check(contig, start, end);

        // カウント結果を表示
        System.out.println("Forward Head Count: " + result[0]);
        System.out.println("Forward Tail Count: " + result[1]);
        System.out.println("Reverse Head Count: " + result[2]);
        System.out.println("Reverse Tail Count: " + result[3]);
    }
    public ControlChecker(String bamPath){
        // BAMファイルを開く
        SamReaderFactory factory =  SamReaderFactory.makeDefault()
                     .enable(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS, SamReaderFactory.Option.VALIDATE_CRC_CHECKSUMS)
                                   .validationStringency(ValidationStringency.SILENT);
        this.reader = factory.open(new File(bamPath));
    }
    public void close() throws IOException{
        reader.close();
    }
    public int[] check(String chr, int start, int end) {
        // カウントを初期化
        int forwardHeadCount = 0;
        int forwardTailCount = 0;
        int reverseHeadCount = 0;
        int reverseTailCount = 0;
        int[] result = new int[4];

        // 読み込むレコードをフィルタリングする条件を設定
        SamRecordFilter filter = new SamRecordFilter() {
            @Override
            public boolean filterOut(SAMRecord record) {
                // 分類する範囲に含まれるレコードのみを通過させる
          //    if (record.getReferenceName().equals(contig) && record.getAlignmentEnd() >= start && record.getAlignmentStart() <= end) {
          //        return false; // レコードを通過させる
          //    }
          return false;
          //return true; // レコードをフィルタリングする
            }

            @Override
            public boolean filterOut(SAMRecord first, SAMRecord second) {
                return false; // ペアリードをフィルタリングしない
            }
        };

        // レコードをイテレートしながら分類とカウントを行う
        SAMRecordIterator iterator = reader.query(chr, start, end, false);
        while (iterator.hasNext()) {
            SAMRecord record = iterator.next();
            // フィルタリングしてカウントする
            if (!filter.filterOut(record)) {
                // System.err.println(record.getAlignmentStart());
                if (record.getReadNegativeStrandFlag()) {
                    if (record.getAlignmentStart() >= start && record.getAlignmentStart() <= end)
                        reverseHeadCount++;
                    else if (record.getAlignmentEnd() >= start && record.getAlignmentEnd() <= end)
                        reverseTailCount++;
                } else {
                    if (record.getAlignmentStart() >= start && record.getAlignmentStart() <= end)
                        forwardHeadCount++;
                    else if (record.getAlignmentEnd() >= start && record.getAlignmentEnd() <= end)
                        forwardTailCount++;
                }
            }
        }
        iterator.close();

        result[0] = forwardHeadCount;
        result[1] = forwardTailCount;
        result[2] = reverseHeadCount;
        result[3] = reverseTailCount;
        return result;
    }

    public int getDepth(String chr, int pos) {
        int depth = 0;

        // レコードをイテレートしながら分類とカウントを行う
        SAMRecordIterator iterator = reader.query(chr, pos, pos + 1, false);
        while (iterator.hasNext()) {
            SAMRecord record = iterator.next();
            // 指定位置がリードに含まれる場合にカウントアップ
            if (record.getAlignmentStart() <= pos && record.getAlignmentEnd() >= pos) {
                depth++;
            }
        }
        iterator.close();
        return depth;
    }
}
