import HUNUS.AlgoHUNUS;

import java.io.IOException;


public class MainTestHUNUS {

    public static void main(String[] args) throws IOException {
//        String dataset = "MSNBC";
//        String dataset = "bible";
        String dataset = "e_shop";
//        String dataset = "e_shop_10k";
//        String dataset = "MicroblogPCU";
//        String dataset = "OnlineRetail_II_all";
//        String dataset = "OnlineRetail_II_best";


        // The algorithm parameters:
        double thresholdratio = 0.055;

        int[][] targetSequence = {
                {356,-1,10,-1,10,-1,10,-1,-2},//bible
                {7,-1,6,-1,-2},//msnbc
                {1,132,192,-1,28,212,-1,-2},//e_shop
                {12766,-1,39356,-1,-2},//MicroblogPCU_two_items
                {39356},//MicroblogPCU_one_item
                {12766,39356,-1,-2},//MicroblogPCU_one_itemset
                {39295,-1,12766,39356,-1,-2},//MicroblogPCU_complex
                {31471,-1,4578,-1,-2},//OnlineRetail_II_all
                {6173,-1,7416,-1,-2},//OnlineRetail_II_best
        };

        String input = "input/" + dataset + ".txt";
        System.out.println(input);

        AlgoHUNUS algo = new AlgoHUNUS();


        int minlen = 1;
        int maxlen = 20;
        algo.setMinLen(minlen);
        algo.setMaxLen(maxlen);


        int mingap = 0;
        int maxgap = 5;
        algo.setMinGap(mingap);
        algo.setMaxGap(maxgap);

        //set target sequence
        algo.setTargetSequence(targetSequence[2]);


        //the path for saving the patterns found
        String  output = "output/HUNUS/" + "HUNUS_" + dataset + '_' + thresholdratio
                + "_len" + minlen + '-' + maxlen
                + "_gap" + mingap + '-' + maxgap + ".txt";
        System.out.println(output);

        algo.runAlgoHUNUS(input, output, thresholdratio, dataset);
    }
}
