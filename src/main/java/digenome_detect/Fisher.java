package digenome_detect;

import org.broadinstitute.hellbender.utils.*;

public class Fisher {
    private static final double TARGET_TABLE_SIZE = 200.0;
    public static void main(String[] argv){
        int[][] table = new int[2][];
        table[0] = new int[2];
        table[1] = new int[2];
        table[0][0] = Integer.parseInt(argv[0]);
        table[0][1] = Integer.parseInt(argv[1]);
        table[1][0] = Integer.parseInt(argv[2]);
        table[1][1] = Integer.parseInt(argv[3]);
        final int[][] normalizedTable = normalizeContingencyTable(table);
        double p = twoSidedPValue(normalizedTable);
        System.out.println("p = " + p);
    }
    public static double exactTest(int[][] table){
        final int[][] normalizedTable = normalizeContingencyTable(table);
        return twoSidedPValue(normalizedTable);
    }
    public static double twoSidedPValue(final int[][] normalizedTable) {
        double p = FisherExactTest.twoSidedPValue(normalizedTable);
        return p;
    }
    private static int[][] normalizeContingencyTable(final int[][] table) {
        final int sum = addExact(table[0][0], table[0][1], table[1][0], table[1][1]);
        if ( sum <= TARGET_TABLE_SIZE * 2 ) {
            return table;
        }

        final double normFactor = sum / TARGET_TABLE_SIZE;

        return new int[][]{
                {(int) (table[0][0] / normFactor), (int) (table[0][1] / normFactor)},
                {(int) (table[1][0] / normFactor), (int) (table[1][1] / normFactor)},
        };
    }

    //Add a bunch of ints, blows up if there's overflow
    private static int addExact(final int... ints){
        int res = ints[0];
        for (int i = 1; i < ints.length; i++) {
            res = Math.addExact(res, ints[i]);
        }
        return res;
    }
}
