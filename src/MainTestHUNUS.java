import HUNUS.AlgoHUNUS;

import java.io.IOException;


public class MainTestHUNUS {

    public static void main(String[] args) throws IOException {
        /**设置数据集和阈值**/
//        String dataset = "test";
//        String dataset = "MSNBC";
        String dataset = "test2";
//        String dataset = "bible";
//        String dataset = "Sign";
//        String dataset = "Leviathan";
//        String dataset = "kosarak";
//        String dataset = "chess_synthetic_3k";
//        String dataset = "e_shop_10";
//        String dataset = "e_shop_1";
//        String dataset = "e_shop";
//        String dataset = "e_shop_10k";
//        String dataset = "MicroblogPCU";
//        String dataset = "OnlineRetail_II_all";
//        String dataset = "OnlineRetail_II_all_300";
//        String dataset = "OnlineRetail_II_best";


        // The algorithm parameters:
        double thresholdratio = 0.055;
//        double thresholdratio = 0.3;//test
//        double thresholdratio = 0.030;//MSNBC
//        double thresholdratio = 0.010;
//        double thresholdratio = 0.014;//BIBLE
//        double thresholdratio = 0.045;
//        double thresholdratio = 0.000230;//相当于bible的0.010
//        double thresholdratio = 130;//e_shop_10
//        double thresholdratio = 0.90;//e_shop_1
//        double thresholdratio = 0.15;//e_shop
//        double thresholdratio = 0.50;//e_shop_10k
//        double thresholdratio = 0.91;//MicroblogPCU
//        double thresholdratio = 0.35;//OnlineRetail_II_all
//        double thresholdratio = 0.4;//OnlineRetail_II_best

        int[][] targetSequence = {
                {356,-1,10,-1,10,-1,10,-1,-2},//bible
                {8,-1,17,-1,8,-1,-2},//leviathan
                {7,-1,6,-1,-2},//msnbc
                {1857,4250,-1,-2},//syn
                {3,-1,5,-1,7,-1,-2},//chess_synthetic_3k [3,5] [3,5,7] [3,5,7,40]
                {8,-1,9,-1,-2},//sign
                {11,-1,218,-1,6,-1,148,-1,-2},//kosarak 10k
                {1,-1,3,5,-1,3,-1,-2},//test
                {1,-1,3,-1,-2},//test2
                {1,-1,-2},//test3
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
        int maxlen = 5;
        algo.setMinLen(minlen);
        algo.setMaxLen(maxlen);


        int mingap = 0;
        int maxgap = 3;
        algo.setMinGap(mingap);
        algo.setMaxGap(maxgap);

        //set target sequence
        algo.setTargetSequence(targetSequence[8]);


        //the path for saving the patterns found
        String  output = "output/HUNUS/" + "HUNUS_" + dataset + '_' + thresholdratio
                + "_len" + minlen + '-' + maxlen
                + "_gap" + mingap + '-' + maxgap + ".txt";
        System.out.println(output);
//        String  output = "output/HUNUS/" + "HUNUS_" + dataset + '_' + thresholdratio
//                + "_len" + minlen + '-' + maxlen
//                + "_gap" + mingap + '-' + maxgap + "complex" + ".txt";
//        System.out.println(output);

        algo.runAlgoHUNUS(input, output, thresholdratio, dataset);
    }
}
