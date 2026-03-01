package HUNUS;

import java.io.*;
import java.util.*;
import java.util.AbstractMap;


public class AlgoHUNUS {

    //the time the algorithm started
    double starTimestamp = 0;
    //the time the algorithm terminated
    double endTimestamp = 0;
    int NumOfCandidate = 0;
    //the number of patterns generated
    int patternCount = 0;

    int minlen = 0;
    int maxlen = 0;

    int mingap = 0;
    int maxgap = 0;

    private String datasetName = "";
    //writer to write the output file
    BufferedWriter writer = null;
    //writer to write the qmatrix debug file
    BufferedWriter qmatrixWriter = null;

    final int BUFFERS_SIZE = 2000;

    //if true, debugging information will be shown in the console
    final int DEBUG = 0;
    int DEBUG_flag = 0;

    //the minUtility threshold
    double minUtility = 0.0;

    //buffer for storing the current pattern that is mined when performing mining
    private int[] patternBuffer = null;

    //the input file path
    String input;

    //Target sequence
    int[] targetSequence;


    private ArrayList<Node>[] reusableNettree = null;
    private int currentNettreeSize = 0;

    ArrayList<ArrayList<Integer>> targetseq;

    //LI_T(LI-Table)
    ArrayList<ArrayList<Integer>> LI_T = new ArrayList<>();

    ArrayList<QMatrix> database = new ArrayList<>();

    private final List<ResultPattern> foundPatterns = new ArrayList<>();

    public AlgoHUNUS() {
    }

    //Inner class representing a node
    static class Node{
        //The corresponding position of a node in sequence
        int name;

        int level;

        //The position of minnum leaf node
        int minLeave;

        //The position of maxnum leaf node
        int maxLeave;

        //The position set of parents
        List<Node> parent = new ArrayList<>();

        //The position set of parents
        List<Node> children = new ArrayList<>();

        //True is having used, false is having not used
        boolean used = false;

        //True is can reach leaves, false is not
        boolean toleave = false;
    }

    static class Instance{
        int sid;
        ArrayList<Integer> path;
        long utility;

        public Instance(int sid, ArrayList<Integer> path, long utility){
            this.sid = sid;
            this.path = path;
            this.utility = utility;
        }

        public int getStartItemsetIdx(){
            return path.getFirst();
        }
        public int getEndItemsetIdx(){
            return path.getLast();
        }
    }


    static class TargetedList {
        public static class UtilityElement {
            public int tid;
            public int acu;
            public int ru;

            public UtilityElement(int tid, int acu, int ru) {
                this.tid=tid;
                this.acu=acu;
                this.ru=ru;
            }
        }


        public int sid;
        public int prel;
        List<UtilityElement> List = new ArrayList<>();
        public int LengthOfUtilityList;
        public int SEU;

        public TargetedList() {
            this.LengthOfUtilityList=0;
            this.SEU = 0;
        }

        public void add(int tid, int acu, int ru) {
            this.List.add(new UtilityElement(tid, acu, ru));
            this.LengthOfUtilityList++;
        }

        public void set_SEU(int seu) {
            this.SEU = seu;
        }

        public void set_sid(int sid) {
            this.sid = sid;
        }

        public int get_sid() {
            return this.sid;
        }

        public void set_prel(int prel) {
            this.prel = prel;
        }

        public int get_prel() {
            return this.prel;
        }
    }


    static class ExtensionKey {
        final int item;
        final String type; // "i" or "s"

        public ExtensionKey(int item, String type) {
            this.item = item;
            this.type = type;
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;
            ExtensionKey that = (ExtensionKey) o;
            return item == that.item && (type != null ? type.equals(that.type) : that.type == null);
        }

        @Override
        public int hashCode() {
            int result = item;
            result = 31 * result + (type != null ? type.hashCode() : 0);
            return result;
        }
    }



    /**
     * @param input the input file path
     * @param output the output file path
     * @param utilityratio the minimum utility threshold ratio
     * @throws IOException exception if error while writing the file
     */
    public void runAlgoHUNUS(String input, String output, double utilityratio, String dataset) throws IOException {
        this.datasetName = dataset;

        //reset maximum
        MemoryLogger.getInstance().reset();

        //input path
        this.input = input;

        //initialize the buffer for staring the current itemset
        patternBuffer = new int[BUFFERS_SIZE];

        //record the start time of the algorithm
        starTimestamp = System.currentTimeMillis();

        //create a writer object to write results to file
        writer = new BufferedWriter(new FileWriter(output));

        //for staring the current sequence number
        int NumberOfSequence = 0;

        //for storing the utility of all sequence
        int totalUtility = 0;

        BufferedReader myInput =null;
        String thisLine;

        //================  First DATABASE SCAN TO CONSTRUCT QMATRIX 	===================

        System.out.println(targetseq);

        try{
            //prepare the object for reading the file
            myInput = new BufferedReader(new InputStreamReader(new FileInputStream(new File(input))));

            //Read each sequence in buffers
            //The first buffer will store the items of a sequence and the -1 between them
            int[] itemBuffer = new int[BUFFERS_SIZE];
            //The following variable will contain the length of the data stored in the previous buffer.
            int itemBufferLength;
            // Finally, we create another buffer for storing the items from a sequence without
            // the -1. This is just used so that we can collect the list of items in that sequence
            // efficiently. We will use this information later to create the number of rows in the
            // QMatrix for that sequence.
            ArrayList<Integer> itemsSequenceBuffer = new ArrayList<>();

            int sequm = 0;
            //for each line (transaction) until the end of file
            while((thisLine = myInput.readLine()) != null){
                sequm++;
                //if the line is a comment, is empty or is a kind of metadata
                if(thisLine.isEmpty() || thisLine.charAt(0) == '#' || thisLine.charAt(0) == '%' || thisLine.charAt(0) == '@'){
                    continue;
                }

                // We reset the two following buffer length to zero because we are reading a new sequence.
                itemBufferLength = 0;
                itemsSequenceBuffer.clear();

                //split the sequence according to the " " separator
                String[] tokens = thisLine.split(" ");

                //get the sequence utility (the last taken on the line)
                String sequenceUtilityString = tokens[tokens.length - 1];
                int positionColons = sequenceUtilityString.indexOf(':');
                int sequenceUtility = Integer.parseInt(sequenceUtilityString.substring(positionColons + 1));


                ArrayList<ArrayList<Integer>> currentSeq = new ArrayList<>();
                ArrayList<ArrayList<Integer>> currentSeqUtils = new ArrayList<>();
                ArrayList<Integer> currentItemset = new ArrayList<>();
                ArrayList<Integer> currentItemsetUtils = new ArrayList<>();
                //This variable will count the number of itemsets
                int nbItemsets = 1;

                for (int i = 0; i < tokens.length - 4; i++) {
                    String currentToken = tokens[i];
                    if (currentToken.isEmpty()) {
                        continue;
                    }

                    if (currentToken.equals("-1")) {
                        currentSeq.add(currentItemset);
                        currentSeqUtils.add(currentItemsetUtils);

                        currentItemset = new ArrayList<>();
                        currentItemsetUtils = new ArrayList<>();
                        nbItemsets++;
                    }else {
                        int positionLeftBracketString = currentToken.indexOf('[');
                        int positionRightBracketString = currentToken.indexOf(']');
                        String itemString = currentToken.substring(0, positionLeftBracketString);
                        Integer item = Integer.parseInt(itemString);
                        String utilityString = currentToken.substring(positionLeftBracketString + 1, positionRightBracketString);
                        Integer itemUtility = Integer.parseInt(utilityString);

                        currentItemset.add(item);
                        currentItemsetUtils.add(itemUtility);

                        itemsSequenceBuffer.add(item);
                    }
                }

                currentSeq.add(currentItemset);
                currentSeqUtils.add(currentItemsetUtils);

                // SPP strategy
                if (sequenceUtility == 0 || !seqContain(currentSeq, targetseq)) {
                    continue;
                }

                ArrayList<Integer> lastInstancePath = findLastInstancePathWithFillLIT(targetseq, currentSeq);

                if(lastInstancePath != null){
                    totalUtility += sequenceUtility;

                    fillLITableForSequence(this.LI_T, NumberOfSequence, lastInstancePath);
                    //Now, we sort the buffer for storing all items from the current sequence in alphabetical order
                    Collections.sort(itemsSequenceBuffer);
                    //but an item may appear multiple times in that buffer so we will loop over the buffer to remove duplicates
                    //This variable remember the last insertion read in that buffer
                    int newItemsPos = 0;
                    //This variable remember the last item read in that buffer
                    int lastItemSeen = -999;
                    //for each position in that buffer
                    for(int i = 0; i < itemsSequenceBuffer.size(); i++){
                        //get the item
                        int item = itemsSequenceBuffer.get(i);
                        //if the item was not seen previously
                        if(item != lastItemSeen){
                            //we copy it at the current insertion position
                            itemsSequenceBuffer.set(newItemsPos++, item);
                            //we remember this item as the last seen item
                            lastItemSeen = item;
                        }
                    }//remove repeating items of itemsSequenceVuffer (length: newItemsPos) and sort it in ascending order
                    while(itemsSequenceBuffer.size() > newItemsPos){
                        itemsSequenceBuffer.remove(itemsSequenceBuffer.size() - 1);
                    }
                    int nbItems = newItemsPos;

                    int[] itemNamesArray = new int[nbItems];
                    for(int i = 0; i < nbItems; i++){
                        itemNamesArray[i] = itemsSequenceBuffer.get(i);
                    }
                    QMatrix matrix = new QMatrix(nbItems, nbItemsets, itemNamesArray, nbItems, sequenceUtility, currentSeq);
                    this.database.add(matrix);

                    int nowsequenceUtility = sequenceUtility;

                    for (int itemsetIdx = 0; itemsetIdx < currentSeq.size(); itemsetIdx++) {
                        ArrayList<Integer> itemset = currentSeq.get(itemsetIdx);
                        ArrayList<Integer> itemsetUtils = currentSeqUtils.get(itemsetIdx);

                        Map<Integer, Integer> utilityMapForCurrentItemset = new HashMap<>();
                        for (int i = 0; i < itemset.size(); i++) {
                            utilityMapForCurrentItemset.put(itemset.get(i), itemsetUtils.get(i));
                        }

                        for (int rowIdx = 0; rowIdx < matrix.itemNames.length; rowIdx++) {
                            int item = matrix.itemNames[rowIdx];

                            if (utilityMapForCurrentItemset.containsKey(item)) {
                                int utility = utilityMapForCurrentItemset.get(item);
                                nowsequenceUtility -= utility;
                                matrix.registerItem(rowIdx, itemsetIdx, utility, nowsequenceUtility);
                            } else {
                                matrix.registerItem(rowIdx, itemsetIdx, 0, nowsequenceUtility);
                            }
                        }
                    }
                }
                // we update the number of transactions；
                NumberOfSequence++;
            }
            System.out.println("totalUtility:" + totalUtility);
            System.out.println("num of q-seq that contain target sequence:" + NumberOfSequence);


        } catch (Exception e) {
            //catches exception if error while reading the input file
            e.printStackTrace();
        }finally{
            if(myInput != null){
                //close the input file
                myInput.close();
            }
        }
        System.out.println("time cost of loading data:" + (System.currentTimeMillis() - starTimestamp)/1000 + "s");
        //check the memory usage
        MemoryLogger.getInstance().checkMemory();


        this.minUtility = (totalUtility * utilityratio);
        System.out.println("数据库总效用(Database Total Utility):" + totalUtility);
        System.out.println("效用阈值(Minimum Utility):" + minUtility);


        Map<Integer, Long> mapItemToUtility = new HashMap<>();
        Map<Integer, Long> mapItemToSEU = new HashMap<>();

        for(int sid = 0; sid < this.database.size(); sid++){
            QMatrix qm =  this.database.get(sid);

            Map<Integer, Long> mapItemToMaxPotentionalInSeq = new HashMap<>();

            for(int itemsetIdx = 0; itemsetIdx < qm.getNbItemsets(); itemsetIdx++){
                for(Integer item : qm.getPresentItemsInItemset(itemsetIdx)){

                    long utility = qm.getUtility(itemsetIdx, item);
                    mapItemToUtility.put(item, mapItemToUtility.getOrDefault(item, 0L) + utility);

                    int firstSetPos = 0;
                    int itemFirstPos = itemsetIdx;

                    if(sid < this.LI_T.size() && this.LI_T.get(sid) != null &&
                            !this.LI_T.get(sid).isEmpty()){
                        if(this.LI_T.get(sid).get(firstSetPos) > itemFirstPos ||
                                (this.LI_T.get(sid).get(firstSetPos) == itemFirstPos &&
                                        isSextFeasibleForNextTargetSet_S(item, firstSetPos)))
                        {
                            long restUtility = qm.getRemainingUtility(itemsetIdx, item);
                            long currentPotentional = utility + restUtility;

                            long maxPotentional = mapItemToMaxPotentionalInSeq.getOrDefault(item, 0L);
                            if(currentPotentional > maxPotentional){
                                mapItemToMaxPotentionalInSeq.put(item, currentPotentional);
                            }
                        }
                    }
                }
            }

            for(Map.Entry<Integer, Long> entry : mapItemToMaxPotentionalInSeq.entrySet()){
                mapItemToSEU.put(entry.getKey(),
                        mapItemToSEU.getOrDefault(entry.getKey(), 0L) + entry.getValue());
            }
        }

        Set<Integer> promisingOneSequences = new HashSet<>();

        for (Map.Entry<Integer, Long> entry : mapItemToSEU.entrySet()) {
            if (entry.getValue() >= this.minUtility) {
                promisingOneSequences.add(entry.getKey());
                NumOfCandidate++;
            }
        }


        Map<Integer, ArrayList<TargetedList>> mapItemUC = new HashMap<>();
        for (int sid = 0; sid < this.database.size(); sid++) {
            QMatrix qm = this.database.get(sid);

            for (int item : qm.getItemNames()) {
                if (!promisingOneSequences.contains(item)) {
                    continue;
                }

                TargetedList ul = new TargetedList();
                ul.set_sid(sid);
                long localSEU = 0;

                int row = qm.getItemIndex(item);
                if (row == -1) continue;

                for (int k = 0; k < qm.getNbItemsets(); k++) {
                    int itemUtility = qm.getUtility(k, item);
                    if (itemUtility > 0) {
                        int prel = 0;
                        int newprel = updateprel(item, prel, 1);

                        int firstSetPos = getFirstSetPos(newprel, this.targetSequence);
                        if (firstSetPos == 99999 || this.LI_T.get(sid).get(firstSetPos) > k ||
                                (this.LI_T.get(sid).get(firstSetPos) == k && isSextFeasibleForNextTargetSet_S(item, firstSetPos))) {
                            ul.set_prel(newprel);
                            int acu = itemUtility;
                            int ru = qm.getRemainingUtility(k, item);
                            ul.add(k, acu, ru);

                            if (acu + ru > localSEU) localSEU = acu + ru;
                        }
                    }
                }

                if (ul.LengthOfUtilityList > 0) {
                    ul.set_SEU((int) localSEU);

                    ArrayList<TargetedList> chain = mapItemUC.getOrDefault(item, new ArrayList<>());
                    chain.add(ul);
                    mapItemUC.put(item, chain);
                }
            }
        }

        for(Integer item : promisingOneSequences){
            if (!mapItemUC.containsKey(item)) {
                continue;
            }

            ArrayList<ArrayList<Integer>> pattern = new ArrayList<>();
            ArrayList<Integer> itemset = new ArrayList<>();
            itemset.add(item);
            pattern.add(itemset);

            ArrayList<TargetedList> initialChain = mapItemUC.get(item);

            long u_no_1_seq = mapItemToUtility.getOrDefault(item, 0L);

            if (u_no_1_seq >= this.minUtility && patternContainsTarget(pattern, this.targetseq)) {
                foundPatterns.add(new ResultPattern(pattern, u_no_1_seq));
                patternCount++;
            }

            patternGrowth(pattern, initialChain);
        }

        endTimestamp = System.currentTimeMillis();
        double runtime = endTimestamp - starTimestamp;
        StringBuilder builder = new StringBuilder();
        builder.append(this.datasetName).append("\n");
        builder.append("totalUtility: ").append(totalUtility).append("\n");
        builder.append("utilityratio: ").append(utilityratio).append("\n");
        builder.append("minUtility: ").append(this.minUtility).append("\n");
        builder.append("target sequence: ").append(targetseq).append("\n");
        builder.append("[minlen,maxlen]: ").append(this.minlen).append("~").append(this.maxlen).append("\n");
        builder.append("[mingap,maxgap]: ").append(this.mingap).append("~").append(this.maxgap).append("\n");
        builder.append("pattern count: ").append(patternCount).append("\n\n");
        builder.append("runtime: ").append(runtime / 1000.0).append(" s\n");
        builder.append("MemoryLogger: ").append(MemoryLogger.getInstance().getMaxMemory()).append(" MB\n");
        builder.append("NumOfCandidate: ").append(NumOfCandidate).append("\n");


        builder.append("====");

        writer.write(builder.toString());

        writer.close();

        if(qmatrixWriter != null){
            qmatrixWriter.close();
        }

        printStatistics(runtime);
    }

    /**
     *
     * @param pattern
     * @param utilitychain
     */
    private void patternGrowth(ArrayList<ArrayList<Integer>> pattern,
                               ArrayList<TargetedList> utilitychain)
    {
        ArrayList<Integer> lastItemsetOfPattern = pattern.getLast();
        int lastItemOfPattern = lastItemsetOfPattern.getLast();


        Map<ExtensionKey, Long> extensionsMap = new HashMap<>();


        for (TargetedList utilitylist : utilitychain) {
            if (utilitylist.List.isEmpty()) continue;

            int sid = utilitylist.get_sid();
            QMatrix qm = this.database.get(sid);
            long localSEU = utilitylist.SEU;
            int prel = utilitylist.get_prel();


            for (int j = 0; j < utilitylist.LengthOfUtilityList; j++) {
                int column = utilitylist.List.get(j).tid;

                for (int item : qm.getPresentItemsInItemset(column)) {
                    if (item > lastItemOfPattern) {
                        int newprel = updateprel(item, prel, 0); // 0 for i-ext
                        int firstSetPos = getFirstSetPos(newprel, this.targetSequence);
                        if (firstSetPos == 99999 || this.LI_T.get(sid).get(firstSetPos) > column ||
                                (this.LI_T.get(sid).get(firstSetPos) == column &&
                                        isIextFeasibleForNextTargetSet_I(item, prel))) {
                            ExtensionKey key = new ExtensionKey(item, "i");
                            extensionsMap.put(key, extensionsMap.getOrDefault(key, 0L) + localSEU);
                        }
                    }
                }
            }

            if (utilitylist.List.isEmpty()) continue;
            int startColumn = utilitylist.List.getFirst().tid + 1;

            for (int item : qm.getItemsAfter(startColumn)) {
                int newprel = updateprel(item, prel, 1); // 1 for s-ext
                int itemFirstPos = qm.getFirstOccurrence(item, startColumn);

                if (itemFirstPos != -1) {
                    int firstSetPos = getFirstSetPos(newprel, this.targetSequence);
                    if(firstSetPos == 99999 || this.LI_T.get(sid).get(firstSetPos) > itemFirstPos ||
                            (this.LI_T.get(sid).get(firstSetPos) == itemFirstPos &&
                                    isSextFeasibleForNextTargetSet_S(item, firstSetPos)))
                    {
                        ExtensionKey key = new ExtensionKey(item, "s");
                        extensionsMap.put(key, extensionsMap.getOrDefault(key, 0L) + localSEU);
                    }
                }
            }
        }


        for (Map.Entry<ExtensionKey, Long> entry : extensionsMap.entrySet()) {
            int item = entry.getKey().item;
            String type = entry.getKey().type;
            long seu = entry.getValue();

            if (seu < this.minUtility) {
                continue;
            }

            NumOfCandidate++;
            ArrayList<ArrayList<Integer>> childPattern = extendPattern(pattern, item, type);

            if (childPattern.size() > this.maxlen) {
                continue;
            }

            Map.Entry<ArrayList<TargetedList>, Long> projection = constructProjection(
                    item, type, utilitychain);

            ArrayList<TargetedList> childChain = projection.getKey();
            long childSEU = projection.getValue();

            if (childChain.isEmpty()) {
                continue;
            }

            boolean containsTarget = patternContainsTarget(childPattern, this.targetseq);

            if (containsTarget) {
                long u_no_valid = calculateNonOverlappingUtility_Projected(childPattern, childChain);

                if (u_no_valid >= this.minUtility) {
                    foundPatterns.add(new ResultPattern(childPattern, u_no_valid));
                    patternCount++;
                }
            }

            if (childSEU >= this.minUtility) {
                patternGrowth(childPattern, childChain);
            }

            childChain.clear();
        }
        extensionsMap.clear();
    }

    public void setTargetSequence(int[] targetSequence){
        this.targetSequence = targetSequence;
        this.targetseq = convertT(targetSequence);
    }

    public void setMinLen(int minLen){
        this.minlen = minLen;
    }

    public void setMaxLen(int maxLen){
        this.maxlen = maxLen;
    }

    public void setMinGap(int mingap){
        this.mingap = mingap;
    }

    public void setMaxGap(int maxGap){
        this.maxgap = maxGap;
    }


    public int updateprel(int exitem, int prel, int kind){
        if(prel >= this.targetSequence.length - 2) return prel;

        if(kind==0){ // i-extension
            if(this.targetSequence[prel]==-1) return prel;
            else {
                if (this.targetSequence[prel] < exitem) {
                    while (prel > 0 && this.targetSequence[prel - 1] != -1) {
                        prel--;
                    }
                    return prel;
                }
                else if (this.targetSequence[prel] > exitem) return prel;
                else return prel + 1;
            }
        }
        else{ // s-extension
            if(this.targetSequence[prel]==-1){
                prel++;
                return updateprel(exitem, prel, 0);
            }
            else{
                while(prel>0 && this.targetSequence[prel-1]!=-1) { prel--; }
                return updateprel(exitem, prel, 0);
            }
        }
    }


    public int getFirstSetPos(int prel, int []T){
        if(prel >= T.length - 2) return 99999;
        else{
            int curitemset = 0;
            for(int i=0; i<=prel; i++){
                if(T[i]==-1) curitemset++;
            }
            return curitemset;
        }
    }

    private boolean isSextFeasibleForNextTargetSet_S(int item, int firstSetPos){
        if(firstSetPos == 99999){
            return true;
        }
        ArrayList<Integer> targetSet = this.targetseq.get(firstSetPos);

        for(Integer y : targetSet){
            if(y < item){
                return false;
            }
        }
        return true;
    }


    private boolean isIextFeasibleForNextTargetSet_I(int itemToExtend, int prel) {
        if (prel >= this.targetSequence.length - 2) {
            return true;
        }

        int nextTargetItem = this.targetSequence[prel];

        if (nextTargetItem == -1) {
            return true;
        }
        return itemToExtend <= nextTargetItem;
    }


    private AbstractMap.SimpleEntry<ArrayList<TargetedList>, Long> constructProjection(
            int item, String type,
            ArrayList<TargetedList> utilitychain)
    {
        ArrayList<TargetedList> newUChain = new ArrayList<>();
        long globalSEU = 0;

        for (TargetedList utilitylist : utilitychain) {
            int sid = utilitylist.get_sid();
            QMatrix qm = this.database.get(sid);

            int row = qm.getItemIndex(item);
            if (row == -1) continue;

            TargetedList newUL = new TargetedList();
            newUL.set_sid(sid);

            long localSEU = 0;

            if ("i".equals(type)) {
                for (int j = 0; j < utilitylist.LengthOfUtilityList; j++) {
                    TargetedList.UtilityElement ue = utilitylist.List.get(j);
                    int column = ue.tid;

                    int itemUtility = qm.getUtility(column, item);

                    if (itemUtility != 0) {
                        int prel = utilitylist.get_prel();
                        int newprel = updateprel(item, prel, 0); // 0 for i-ext

                        int firstSetPos = getFirstSetPos(newprel, this.targetSequence);
                        if (firstSetPos == 99999 || this.LI_T.get(sid).get(firstSetPos) > column ||
                                (this.LI_T.get(sid).get(firstSetPos) == column && isIextFeasibleForNextTargetSet_I(item, prel))) {
                            newUL.set_prel(newprel);
                            int acu = ue.acu + itemUtility;
                            int ru = qm.getRemainingUtility(column, item);
                            newUL.add(column, acu, ru);

                            if (acu + ru > localSEU)
                                localSEU = acu + ru;
                        }
                    }
                }
            } else {
                if(!utilitylist.List.isEmpty()){
                    TargetedList.UtilityElement firstUe = utilitylist.List.getFirst();
                    int startColumn = firstUe.tid + 1;

                    for (int k = startColumn; k < qm.getNbItemsets(); k++) {
                        int itemUtility = qm.getUtility(k, item);

                        if (itemUtility != 0) {
                            int prel = utilitylist.get_prel();
                            int newprel = updateprel(item, prel, 1); // 1 for s-ext

                            int firstSetPos = getFirstSetPos(newprel, this.targetSequence);
                            if (firstSetPos == 99999 || this.LI_T.get(sid).get(firstSetPos) > k ||
                                    (this.LI_T.get(sid).get(firstSetPos) == k && isSextFeasibleForNextTargetSet_S(item, firstSetPos))) {
                                newUL.set_prel(newprel);
                                int acu = firstUe.acu + itemUtility;
                                int ru = qm.getRemainingUtility(k, item);
                                newUL.add(k, acu, ru);

                                if (acu + ru > localSEU) localSEU = acu + ru;
                            }
                        }
                    }
                }


            }

            if (newUL.LengthOfUtilityList > 0) {
                newUL.set_SEU((int) localSEU);
                newUChain.add(newUL);
                globalSEU += localSEU;
            }
        }

        return new AbstractMap.SimpleEntry<>(newUChain, globalSEU);
    }


    private long calculateNonOverlappingUtility_Projected(
            ArrayList<ArrayList<Integer>> pattern,
            ArrayList<TargetedList> utilitychain)
    {
        long u_no_valid = 0;
        int patternSize = pattern.size();

        checkAndResizeNettree(patternSize);

        for (TargetedList ul : utilitychain) {
            int sid = ul.get_sid();
            QMatrix qm = this.database.get(sid);
            ArrayList<ArrayList<Integer>> sequence = qm.getSequenceAsItemsets();

            cleanNettree(patternSize);

            createNetTree(this.reusableNettree, sequence, pattern);

            updateNetTreeWithLeaves(this.reusableNettree, patternSize);

            ArrayList<Instance> nonOverlappingInSeq = findNonOverlappingInSequence(this.reusableNettree, pattern, sid, qm);

            for (Instance inst : nonOverlappingInSeq) {
                u_no_valid += inst.utility;
            }
        }
        return u_no_valid;
    }



    public boolean seqContain(ArrayList<ArrayList<Integer>> sequenceA, ArrayList<ArrayList<Integer>> subSequence) {
        if (subSequence.isEmpty()) {
            return true;
        }
        if (sequenceA.isEmpty()) {
            return false;
        }

        int subIdx = 0;
        int mainIdx = 0;

        while (subIdx < subSequence.size() && mainIdx < sequenceA.size()) {
            if (sequenceA.get(mainIdx).containsAll(subSequence.get(subIdx))) {
                subIdx++;
            }
            mainIdx++;
        }

        return subIdx == subSequence.size();
    }

    /**
     * transfer array A to a list
     * @param T target sequence
     * @return list T.
     */
    public ArrayList<ArrayList<Integer>> convertT(int[] T) {
        if (T == null) {
            throw new IllegalArgumentException("T can not null");
        }
        ArrayList<ArrayList<Integer>> B = new ArrayList<>();
        int i = 0;
        ArrayList<Integer> templist = new ArrayList<>();
        while(i < T.length && T[i] != -2){
            if(T[i] == -1){
                B.add(templist);
                templist = new ArrayList<>();
                i++;
            }else{
                templist.add(T[i]);
                i++;
            }
        }
        return B;
    }


    void createNetTree(ArrayList<Node>[] nettree, ArrayList<ArrayList<Integer>> sequence,
                       ArrayList<ArrayList<Integer>> pattern){
        int patternLayerCount = pattern.size();

        int[] start = new int[patternLayerCount];

        for(int itemsetIdx = 0; itemsetIdx < sequence.size(); itemsetIdx++){
            ArrayList<Integer> currentSeqItemset = sequence.get(itemsetIdx);

            List<Integer> patternItem0 = pattern.getFirst();
            if(currentSeqItemset.containsAll(patternItem0)){
                Node newNode = new Node();
                newNode.name = itemsetIdx;
                newNode.level = 0;
                newNode.minLeave = itemsetIdx;
                newNode.maxLeave = itemsetIdx;
                nettree[0].add(newNode);
            }

            for(int j = 1; j < patternLayerCount; j++){
                ArrayList<Integer> patterItem_j = pattern.get(j);

                if(currentSeqItemset.containsAll(patterItem_j)){
                    ArrayList<Node> parentLayer = nettree[j - 1];
                    if(parentLayer.isEmpty()){
                        continue;
                    }

                    while(start[j-1] < parentLayer.size() &&
                            (itemsetIdx - parentLayer.get(start[j-1]).name -1) > this.maxgap){
                        start[j-1]++;
                    }

                    if(start[j-1] >= parentLayer.size() ||
                            (itemsetIdx - parentLayer.getFirst().name - 1  < this.mingap)){
                        continue;
                    }

                    Node childNode = new Node();
                    childNode.name = itemsetIdx;
                    childNode.level = j;
                    childNode.minLeave = itemsetIdx;
                    childNode.maxLeave = itemsetIdx;
                    nettree[j].add(childNode);

                    for(int k = start[j-1]; k < parentLayer.size(); k++){
                        Node parentNode = parentLayer.get(k);
                        int currentGap = itemsetIdx - parentNode.name - 1;

                        if(currentGap < this.mingap){
                            break;
                        }

                        parentNode.children.add(childNode);
                        childNode.parent.add(parentNode);
                    }
                }
            }
        }
    }


    private ArrayList<Integer> findLastInstancePathWithFillLIT(
            ArrayList<ArrayList<Integer>> target,
            ArrayList<ArrayList<Integer>> sequence) {
        int startScanIdx = sequence.size() - 1;

        while (startScanIdx >= 0) {
            int[] lastpos = new int[target.size()];
            Arrays.fill(lastpos, -1);

            int targetIdx = target.size() - 1;

            for (int seqIdx = startScanIdx; seqIdx >= 0; seqIdx--) {
                if (targetIdx < 0) {
                    break;
                }

                if (sequence.get(seqIdx).containsAll(target.get(targetIdx))) {
                    lastpos[targetIdx] = seqIdx;
                    targetIdx--;
                }
            }


            if (targetIdx >= 0) {
                return null;
            }

            boolean isValid = true;
            int instanceLength = lastpos[lastpos.length - 1] - lastpos[0] + 1;
            if (instanceLength > this.maxlen) {
                isValid = false;
            }

            if (isValid) {
                for (int i = 0; i < lastpos.length - 1; i++) {
                    int gap = lastpos[i + 1] - lastpos[i] - 1;
                    if (gap < this.mingap) {
                        isValid = false;
                        break;
                    }
                }
            }

            if (isValid) {
                ArrayList<Integer> path = new ArrayList<>();
                for (int pos : lastpos) {
                    path.add(pos);
                }
                return path;
            } else {
                startScanIdx = lastpos[lastpos.length - 1] - 1;
            }
        }

        return null;
    }

    private ArrayList<Node> findFirstValidInstance(ArrayList<Node>[] nettree,
                                                   ArrayList<ArrayList<Integer>> pattern){
        if(nettree[0].isEmpty()){
            return null;
        }

        for(Node rootNode : nettree[0]){
            //CHECK: Skip if root node is already used
            if(rootNode.used){
                continue;
            }

            ArrayList<Node> path = new ArrayList<>();

            ArrayList<Node> foundPath = findFullPathRecursiveFromFirst(rootNode, rootNode,
                    1, pattern.size(), path);

            if(foundPath != null){
                return foundPath;
            }
        }
        return null;
    }


    private ArrayList<Node> findFullPathRecursiveFromFirst(Node rootNode, Node currentNode,  int currentLevel,
                                                           int targetLevel, ArrayList<Node> path){
        //CHECK: Skip if current node is already used
        if(currentNode.used){
            return null;
        }

        path.add(currentNode);

        if(currentLevel == targetLevel){
            int span = currentNode.name - rootNode.name + 1;
            if(span >= minlen && span <= maxlen){
                ArrayList<Node> result = new ArrayList<>(path);
                path.removeLast();
                return result;
            }else {
                path.removeLast();
                return null;
            }
        }

        for(Node childNode : currentNode.children){
            int currentSpan = childNode.name - rootNode.name + 1;
            if (currentSpan > this.maxlen) {
                break;
            }

            ArrayList<Node> foundPath = findFullPathRecursiveFromFirst(rootNode, childNode,
                    currentLevel + 1, targetLevel, path);
            if(foundPath != null){
                path.removeLast();
                return foundPath;
            }
        }

        path.removeLast();

        return null;
    }


    void updateNetTreeToRemovePath(ArrayList<Node> pathNodes){
        if(pathNodes == null || pathNodes.isEmpty()){
            return;
        }

        for(Node node : pathNodes){
            node.used = true;
        }


        for(int i = pathNodes.size() - 2; i >= 0; i--){
            Node parentNode = pathNodes.get(i);

            if(parentNode.used){
                continue;
            }

            boolean allChildernUsed = true;
            if(parentNode.children.isEmpty()){
                allChildernUsed = true;
            }else{
                for(Node child : parentNode.children){
                    if(!child.used){
                        allChildernUsed = false;
                        break;
                    }
                }
            }

            if(allChildernUsed){
                parentNode.used = true;
            }
        }


    }


    private ArrayList<Instance> findNonOverlappingInSequence(
            ArrayList<Node>[] nettree,
            ArrayList<ArrayList<Integer>> pattern,
            int sid,
            QMatrix qm) {
        ArrayList<Instance> instancesInSeq = new ArrayList<>();


        if (nettree.length > 0 && nettree[0] != null) {
            for (Node rootNode : nettree[0]) {
                if (!rootNode.toleave) {
                    rootNode.used = true;
                    continue;
                }
                int shortestPossibleSpan = rootNode.minLeave - rootNode.name + 1;
                int longestPossibleSpan = rootNode.maxLeave - rootNode.name + 1;

                if (longestPossibleSpan < this.minlen || shortestPossibleSpan > this.maxlen) {
                    rootNode.used = true;
                }
            }
        }

        while (true) {
            ArrayList<Node> nodepath = findFirstValidInstance(nettree, pattern);

            if (nodepath == null) {
                break;
            }

            long instanceUtility = calculateUtilityOfInstance(nodepath, qm, pattern);

            ArrayList<Integer> indexPath = new ArrayList<>();
            for (Node node : nodepath) {
                indexPath.add(node.name);
            }
            instancesInSeq.add(new Instance(sid, indexPath, instanceUtility));

            updateNetTreeToRemovePath(nodepath);
        }
        return instancesInSeq;
    }

    public void fillLITableForSequence(ArrayList<ArrayList<Integer>> LI_T,int NumberOfSequence,
                                       ArrayList<Integer> lastInstancePath){
        while(LI_T.size() <= NumberOfSequence){
            LI_T.add(new ArrayList<>());
        }

        LI_T.get(NumberOfSequence).clear();

        LI_T.get(NumberOfSequence).addAll(lastInstancePath);
    }


    static String arrTostr(int[] seq){
        StringBuilder stringBuilder = new StringBuilder();
        for (int j : seq) {
            stringBuilder.append(j).append(" ");
        }
        return stringBuilder.toString();
    }

    /**
     * Print statistics about the latest execution to System.out.
     */
    public void printStatistics(double runtime){
        System.out.println("============= ALGOHUNUS  v1.0 - STATS =============\n");
        System.out.println("Target sequence: " + arrTostr(targetSequence));
        System.out.println(" Threshold :"+this.minUtility);
        System.out.println(" Total time ~ " + runtime/1000 + " s");
        System.out.println(" Max Memory ~ " + MemoryLogger.getInstance().getMaxMemory() + " MB");
        System.out.println(" High-utility sequential pattern count : " + patternCount);
        System.out.println(" Number of candidates : " + NumOfCandidate + " \n");
        System.out.println("========================================================");

    }


    private long calculateUtilityOfInstance(ArrayList<Node> nodepath,
                                            QMatrix qm,
                                            ArrayList<ArrayList<Integer>> pattern){
        long totalUtility = 0;

        if(nodepath.size() != pattern.size()){
            throw new IllegalArgumentException("path and pattern size do not match");
        }

        for(int i = 0; i < nodepath.size(); i++){
            int itemIdx = nodepath.get(i).name;
            ArrayList<Integer> patternItemset = pattern.get(i);

            for(int patternItem : patternItemset){
                totalUtility += qm.getUtility(itemIdx, patternItem);
            }
        }
        return totalUtility;
    }


    private ArrayList<ArrayList<Integer>> extendPattern(ArrayList<ArrayList<Integer>> parentPattern,
                                                        int item, String type){
        ArrayList<ArrayList<Integer>> newPattern = new ArrayList<>();
        for(ArrayList<Integer> itemset : parentPattern){
            newPattern.add(new ArrayList<>(itemset));
        }

        if("i".equals(type)){
            ArrayList<Integer> lastItemset = newPattern.getLast();
            lastItemset.add(item);
            Collections.sort(lastItemset);
        }else{
            ArrayList<Integer> newItemset = new ArrayList<>();
            newItemset.add(item);
            newPattern.add(newItemset);
        }
        return newPattern;
    }


    private boolean patternContainsTarget(ArrayList<ArrayList<Integer>> pattern,
                                          ArrayList<ArrayList<Integer>> target){
        int targetIdx = 0;
        int patternIdx = 0;
        while(targetIdx < target.size() && patternIdx < pattern.size()){
            if(pattern.get(patternIdx).containsAll(target.get(targetIdx))){
                targetIdx++;
            }
            patternIdx++;
        }
        return targetIdx == target.size();
    }


    static class ResultPattern {
        ArrayList<ArrayList<Integer>> pattern;
        long utility;

        public ResultPattern(ArrayList<ArrayList<Integer>> pattern, long utility) {
            this.pattern = new ArrayList<>();
            for (ArrayList<Integer> itemset : pattern) {
                this.pattern.add(new ArrayList<>(itemset));
            }
            this.utility = utility;
        }
    }


    private void updateNetTreeWithLeaves(ArrayList<Node>[] nettree, int patternSize) {
        if (patternSize == 0) return;


        for (int i = 0; i < patternSize - 1; i++) {
            if (nettree.length > i && nettree[i] != null) {
                for (Node node : nettree[i]) {
                    node.toleave = false;
                }
            }
        }
        if (nettree.length > patternSize - 1 && nettree[patternSize - 1] != null) {
            for (Node node : nettree[patternSize - 1]) {
                node.toleave = true;
            }
        }


        for (int i = patternSize - 2; i >= 0; i--) {
            for (Node node : nettree[i]) {
                if (node.children.isEmpty()) {
                    continue;
                }

                int minL = Integer.MAX_VALUE;
                int maxL = Integer.MIN_VALUE;
                boolean canReachLeaf = false;


                for (Node child : node.children) {
                    if (child.toleave) {
                        canReachLeaf = true;
                        if (child.minLeave < minL) {
                            minL = child.minLeave;
                        }
                        if (child.maxLeave > maxL) {
                            maxL = child.maxLeave;
                        }
                    }
                }

                if (canReachLeaf) {
                    node.toleave = true;
                    node.minLeave = minL;
                    node.maxLeave = maxL;
                }
            }
        }
    }


    private void checkAndResizeNettree(int requiredSize) {
        if (reusableNettree == null || currentNettreeSize < requiredSize) {

            reusableNettree = null;

            reusableNettree = (ArrayList<Node>[]) new ArrayList[requiredSize];

            for (int i = 0; i < requiredSize; i++) {
                reusableNettree[i] = new ArrayList<>();
            }
            currentNettreeSize = requiredSize;
        }
    }


    private void cleanNettree(int layersToClean) {
        int layers = Math.min(layersToClean, this.currentNettreeSize);
        for (int i = 0; i < layers; i++) {
            reusableNettree[i].clear();
        }
    }

}
