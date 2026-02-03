package HUNUS;

import java.io.*;
import java.util.*;
import java.util.AbstractMap;


public class AlgoHUNUS {

    //the time the algorithm started
    double starTimestamp = 0;
    //the time the algorithm terminated
    double endTimestamp = 0;
    //候选者数量
    int NumOfCandidate = 0;
    //the number of patterns generated
    int patternCount = 0;

    int minlen = 0;
    int maxlen = 0;

    int mingap = 0;
    int maxgap = 0;

    private String datasetName = ""; // 用于存储数据集的名称
    //writer to write the output file
    BufferedWriter writer = null;
    //writer to write the qmatrix debug file
    BufferedWriter qmatrixWriter = null;

    final int BUFFERS_SIZE = 2000;

    //if true, debugging information will be shown in the console
    //1.q-seq不包含target seq
    //2.未找到一个实例
    //3.QMatrix
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

    // 用于重用的 Nettree 结构
    private ArrayList<Node>[] reusableNettree = null;
    // 记录当前 reusableNettree 分配的大小
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

        //节点层级信息
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

    /**
     * 内部类，用于封装一个模式实例的信息
     */
    static class Instance{
        int sid;//所在序列的ID
        ArrayList<Integer> path;//实例的项集索引路径
        long utility;

        public Instance(int sid, ArrayList<Integer> path, long utility){
            this.sid = sid;
            this.path = path;
            this.utility = utility;
        }

        //方便获取起止位置
        public int getStartItemsetIdx(){
            return path.getFirst();
        }
        public int getEndItemsetIdx(){
            return path.getLast();
        }
    }

    /**
     * TUSQ 算法中的 TargetedList 结构。
     * 用于存储一个模式在特定序列中的投影信息（一个"效用链"）。
     *
     */
    static class TargetedList {
        /**
         * 内部类，代表效用链表中的一个元素（一个实例的结束点）。
         *
         */
        public static class UtilityElement {
            /** itemset  **/
            public int tid; // 实例结束项集的索引
            /** 效用 **/
            public int acu; // 实例的累积效用 (actual utility)
            /** 剩余效用 **/
            public int ru;  // 该实例结束后的剩余效用

            public UtilityElement(int tid, int acu, int ru) {
                this.tid=tid;
                this.acu=acu;
                this.ru=ru;
            }
        }

        // 序列ID
        public int sid;
        // prel：当前模式包含目标序列的前缀的长度
        public int prel;
        // 效用元素列表
        List<UtilityElement> List = new ArrayList<>();
        // 列表长度
        public int LengthOfUtilityList;
        // SEU
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

    /**
     * [新增] 内部类，用作 extensionsMap 的键。
     * 它可以唯一标识一个扩展（通过 项 + 类型）。
     * 必须正确实现 equals 和 hashCode。
     */
    static class ExtensionKey {
        final int item;
        final String type; // "i" 或 "s"

        public ExtensionKey(int item, String type) {
            this.item = item;
            this.type = type;
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;
            ExtensionKey that = (ExtensionKey) o;
            // Java 1.7+ 可以使用 Objects.equals(type, that.type)
            return item == that.item && (type != null ? type.equals(that.type) : that.type == null);
        }

        @Override
        public int hashCode() {
            // Java 1.7+ 可以使用 Objects.hash(item, type)
            int result = item;
            result = 31 * result + (type != null ? type.hashCode() : 0);
            return result;
        }
    }

    /**
     * 内部类，用于在 TUSQ 逻辑中收集扩展项及其 TDU。
     *
     */
    static class ExtensionTDU {
        long TDU; // Terminated Descendants Utility
        String type; // "i" or "s"

        // i-extension 优先
        public void update(long newTDU, String newType) {
            this.TDU += newTDU;
            if (this.type.equals("s") && newType.equals("i")) {
                this.type = "i";
            }
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

        //create a writer object to write qmatrix debug information to file
        if(DEBUG == 3){
            String qmatrixOutputPath = output.replace(".txt", "_qmatrix_debug.txt");
            qmatrixWriter = new BufferedWriter(new FileWriter(qmatrixOutputPath));
            qmatrixWriter.write("=== QMatrix调试信息 ===\n");
            qmatrixWriter.write("数据集: " + dataset + "\n");
            qmatrixWriter.write("目标序列: " + targetseq + "\n");
            qmatrixWriter.write("阈值比例: " + utilityratio + "\n");
            qmatrixWriter.write("长度约束: [" + minlen + "," + maxlen + "]\n");
            qmatrixWriter.write("间隔约束: [" + mingap + "," + maxgap + "]\n");
            qmatrixWriter.write("========================================\n\n");
        }

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
            // [修复] 使用 ArrayList 代替固定大小的数组，避免数组越界
            ArrayList<Integer> itemsSequenceBuffer = new ArrayList<>();
//            int[] utilitiesSequenceBuffer = new int[BUFFERS_SIZE];

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
                // [修复] 重置 itemsSequenceBuffer 以供新序列使用
                itemsSequenceBuffer.clear();

                //split the sequence according to the " " separator
                String[] tokens = thisLine.split(" ");

                //get the sequence utility (the last taken on the line)
                String sequenceUtilityString = tokens[tokens.length - 1];
                int positionColons = sequenceUtilityString.indexOf(':');
                int sequenceUtility = Integer.parseInt(sequenceUtilityString.substring(positionColons + 1));

                // 累加所有序列的效用（包括不包含目标序列的序列）
//                totalUtility += sequenceUtility;

                // 1. 直接初始化最终目标的数据结构
                ArrayList<ArrayList<Integer>> currentSeq = new ArrayList<>();
                ArrayList<ArrayList<Integer>> currentSeqUtils = new ArrayList<>(); // 同样为效用值创建一个对应结构
                ArrayList<Integer> currentItemset = new ArrayList<>();
                ArrayList<Integer> currentItemsetUtils = new ArrayList<>();
                //This variable will count the number of itemsets
                int nbItemsets = 1;

                // 2. 遍历tokens，直接填充我们的目标结构
                for (int i = 0; i < tokens.length - 4; i++) {
                    String currentToken = tokens[i];
                    if (currentToken.isEmpty()) {
                        continue;
                    }

                    if (currentToken.equals("-1")) {
                        // 当遇到-1时，说明一个项集结束
                        // 将构建好的当前项集添加到主列表中
                        currentSeq.add(currentItemset);
                        currentSeqUtils.add(currentItemsetUtils);

                        // 创建新的空列表以准备下一个项集
                        currentItemset = new ArrayList<>();
                        currentItemsetUtils = new ArrayList<>();
                        nbItemsets++;
                    }else {
                        // 解析item和其效用
                        int positionLeftBracketString = currentToken.indexOf('[');
                        int positionRightBracketString = currentToken.indexOf(']');
                        String itemString = currentToken.substring(0, positionLeftBracketString);
                        Integer item = Integer.parseInt(itemString);
                        String utilityString = currentToken.substring(positionLeftBracketString + 1, positionRightBracketString);
                        Integer itemUtility = Integer.parseInt(utilityString);

                        // 3. 直接将解析出的item和utility加入当前项集列表
                        currentItemset.add(item);
                        currentItemsetUtils.add(itemUtility);

                        // itemsSequenceBuffer 仍然需要，用于后续为QMatrix确定唯一的itemNames
                        // [修复] 使用 ArrayList 的 add 方法，避免数组越界
                        itemsSequenceBuffer.add(item);
                    }
                }

                // 4. 循环结束后，不要忘记添加最后一个正在构建的项集
                currentSeq.add(currentItemset);
                currentSeqUtils.add(currentItemsetUtils);

                // 策略一：根据T限制来过滤
                if (sequenceUtility == 0 || !seqContain(currentSeq, targetseq)) {
                    if (DEBUG == 1) {
                        // ... (debug 代码) ...
                    }
                    continue;
                }

                //策略二：直接自右向左遍历序列，定位目标序列各项集的最后出现位置来填充 LI_T
                ArrayList<Integer> lastInstancePath = findLastInstancePathWithFillLIT(targetseq, currentSeq);

                if(lastInstancePath != null){
                    totalUtility += sequenceUtility;

                    fillLITableForSequence(this.LI_T, NumberOfSequence, lastInstancePath);
                    // *** 为有效序列构建并存储Q-Matrix ***
                    //Now, we sort the buffer for storing all items from the current sequence in alphabetical order
                    // [修复] 使用 Collections.sort 对 ArrayList 进行排序
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
                    // [修复] 移除重复项后，截取有效部分
                    while(itemsSequenceBuffer.size() > newItemsPos){
                        itemsSequenceBuffer.remove(itemsSequenceBuffer.size() - 1);
                    }
                    int nbItems = newItemsPos;

                    // [修复] 将 ArrayList 转换为数组以传递给 QMatrix 构造函数
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

                        // 为当前项集创建一个临时的 item->utility 映射，方便查找
                        Map<Integer, Integer> utilityMapForCurrentItemset = new HashMap<>();
                        for (int i = 0; i < itemset.size(); i++) {
                            utilityMapForCurrentItemset.put(itemset.get(i), itemsetUtils.get(i));
                        }

                        // 遍历QMatrix的所有行（即序列中的所有不重复项）
                        for (int rowIdx = 0; rowIdx < matrix.itemNames.length; rowIdx++) {
                            int item = matrix.itemNames[rowIdx];

                            // 检查该item是否存在于当前项集中
                            if (utilityMapForCurrentItemset.containsKey(item)) {
                                int utility = utilityMapForCurrentItemset.get(item);
                                nowsequenceUtility -= utility;
                                matrix.registerItem(rowIdx, itemsetIdx, utility, nowsequenceUtility);
                            } else {
                                // 如果不存在，效用为0，剩余效用为进入此项集之前的值
                                matrix.registerItem(rowIdx, itemsetIdx, 0, nowsequenceUtility);
                            }
                        }
                    }

                    //if in debug mode, we print the QMatrix that we just built
                    if(DEBUG == 3){
                        System.out.println(matrix.toString());
                        System.out.println();

                        // 同时将QMatrix信息写入调试文件
                        if(qmatrixWriter != null){
                            qmatrixWriter.write("序列 " + NumberOfSequence + ":\n");
                            qmatrixWriter.write(matrix.toString());
                            qmatrixWriter.write("\n");
                            qmatrixWriter.flush(); // 确保数据立即写入文件
                        }
                    }
                }else{
                    if(DEBUG == 2){
                        System.out.print("序列" + sequm + ":没有有效的实例:");
                        for(int i = 0; i < itemBufferLength; i++){
                            System.out.print(itemBuffer[i] + " ");
                        }
                        System.out.println();
                    }
                    //如果返回null，说明不存在任何有效实例，直接过滤掉该序列
                    continue;
                }
                // we update the number of transactions；序列数++，准备处理下一个序列
                NumberOfSequence++;
            }
            System.out.println("totalUtility:" + totalUtility);
            System.out.println("num of q-seq that contain target sequence:" + NumberOfSequence);


            System.out.println("LI_T表：");
            for(int i = 0; i < LI_T.size(); i++){
                System.out.println("序列" + i + LI_T.get(i));
            }

        } catch (Exception e) {
            //catches exception if error while reading the input file
            e.printStackTrace();
        }finally{
            if(myInput != null){
                //close the input file
                myInput.close();
            }
        }//完成读取数据
        System.out.println("time cost of loading data:" + (System.currentTimeMillis() - starTimestamp)/1000 + "s");
        //check the memory usage
        MemoryLogger.getInstance().checkMemory();

        // ================ 第二阶段：计算1-序列全局指标并筛选 =================

        //1.根据阈值比例，计算出具体的最小效用阈值
        this.minUtility = (totalUtility * utilityratio);
//        this.minUtility = 528867;
        System.out.println("数据库总效用(Database Total Utility):" + totalUtility);
        System.out.println("效用阈值(Minimum Utility):" + minUtility);

        //2.准备用于存储全局信息的Map
        //Map<项目，全局真实效用>
        Map<Integer, Long> mapItemToUtility = new HashMap<>();
        //Map<项目，全局SEU>
        Map<Integer, Long> mapItemToSEU = new HashMap<>();

        //3.遍历QMatrix数据库，计算1-序列的全局u_no和SEU
        for(int sid = 0; sid < this.database.size(); sid++){
            QMatrix qm =  this.database.get(sid);

            //用于存储当前序列中，每一个1-序列所能贡献的最大潜力值
            Map<Integer, Long> mapItemToMaxPotentionalInSeq = new HashMap<>();

            //遍历当前序列(qm)中的每一个项集
            for(int itemsetIdx = 0; itemsetIdx < qm.getNbItemsets(); itemsetIdx++){
                //这个循环现在只会遍历那些效用不为0的、真正在当前项集中的item。
                for(Integer item : qm.getPresentItemsInItemset(itemsetIdx)){

                    //--- a. 累加真实效用 u_no ---
                    long utility = qm.getUtility(itemsetIdx, item);
                    mapItemToUtility.put(item, mapItemToUtility.getOrDefault(item, 0L) + utility);

                    //--- b. 计算并更新 SEU ---(通过LI_T)
                    //直接使用更精细的S-extension剪枝逻辑
                    int firstSetPos = 0; //对应目标序列的第一个项集 T_1
                    int itemFirstPos = itemsetIdx;//当前项的出现位置

                    //检查LI_T表是否存在该序列的信息
                    if(sid < this.LI_T.size() && this.LI_T.get(sid) != null &&
                            !this.LI_T.get(sid).isEmpty()){
                        //应用S-extension的剪枝条件
                        if(this.LI_T.get(sid).get(firstSetPos) > itemFirstPos ||
                                (this.LI_T.get(sid).get(firstSetPos) == itemFirstPos &&
                                        isSextFeasibleForNextTargetSet_S(item, firstSetPos)))
                        {
                            long restUtility = qm.getRemainingUtility(itemsetIdx, item);
                            long currentPotentional = utility + restUtility;

                            //更新当前序列中该item的最大潜力值
                            long maxPotentional = mapItemToMaxPotentionalInSeq.getOrDefault(item, 0L);
                            if(currentPotentional > maxPotentional){
                                mapItemToMaxPotentionalInSeq.put(item, currentPotentional);
                            }
                        }
                    }

//                    boolean isPromising = checkIsPromising(item, itemsetIdx, sid, this.targetseq, this.LI_T);
//
//                    if(isPromising){
//                        long restUtility = qm.getRemainingUtility(itemsetIdx, item);
//                        long currentPotentional = utility + restUtility;
//
//                        //更新当前序列中该item的最大潜力值
//                        long maxPotentional = mapItemToMaxPotentionalInSeq.getOrDefault(item, 0L);
//                        if(currentPotentional > maxPotentional){
//                            mapItemToMaxPotentionalInSeq.put(item, currentPotentional);
//                        }
//                    }
                }
            }

            // --- c. 将当前序列中算出的最大潜力值累加到全局SEU中
            for(Map.Entry<Integer, Long> entry : mapItemToMaxPotentionalInSeq.entrySet()){
                mapItemToSEU.put(entry.getKey(),
                        mapItemToSEU.getOrDefault(entry.getKey(), 0L) + entry.getValue());
            }
        }

        //4. 筛选有希望的1-序列，并为启动递归做准备
        Set<Integer> promisingOneSequences = new HashSet<>();

        for (Map.Entry<Integer, Long> entry : mapItemToSEU.entrySet()) {
            //如果SEU不小于阈值，则为有希望的1-序列
            if (entry.getValue() >= this.minUtility) {
                promisingOneSequences.add(entry.getKey());
                NumOfCandidate++;
            }
        }


        System.out.println("有希望的1-序列 (种子): " + promisingOneSequences);

        System.out.println("=== 阶段三：启用递归挖掘 ===");
        // ======================= 第三阶段：启动递归挖掘 =========================
        Map<Integer, ArrayList<TargetedList>> mapItemUC = new HashMap<>();
        for (int sid = 0; sid < this.database.size(); sid++) {
            QMatrix qm = this.database.get(sid);

            for (int item : qm.getItemNames()) {
                // 只为 "有希望的" 1-序列构建初始投影链
                if (!promisingOneSequences.contains(item)) {
                    continue;
                }

                // 获取该 item 在此序列的投影
                TargetedList ul = new TargetedList();
                ul.set_sid(sid);
                long localSEU = 0;

                int row = qm.getItemIndex(item);
                if (row == -1) continue;

                // 查找该 item 的所有出现
                for (int k = 0; k < qm.getNbItemsets(); k++) {
                    int itemUtility = qm.getUtility(k, item);
                    if (itemUtility > 0) {
                        // 检查后缀约束
                        //prel：当前模式包含目标序列的前缀的长度
                        int prel = 0; // 1-sequence 的 prel
                        int newprel = updateprel(item, prel, 1); // 假设1-seq是s-ext

                        int firstSetPos = getFirstSetPos(newprel, this.targetSequence);
                        if (firstSetPos == 99999 || this.LI_T.get(sid).get(firstSetPos) > k ||
                                (this.LI_T.get(sid).get(firstSetPos) == k && isSextFeasibleForNextTargetSet_S(item, firstSetPos))) {
                            ul.set_prel(newprel);
                            int acu = itemUtility;
                            int ru = qm.getRemainingUtility(k, item); // 使用原始RU
                            ul.add(k, acu, ru);

                            if (acu + ru > localSEU) localSEU = acu + ru;
                        }
                    }
                }

                if (ul.LengthOfUtilityList > 0) {
                    ul.set_SEU((int) localSEU);

                    // 添加到全局 mapItemUC
                    ArrayList<TargetedList> chain = mapItemUC.getOrDefault(item, new ArrayList<>());
                    chain.add(ul);
                    mapItemUC.put(item, chain);
                }
            }
        }

        // 2. 对每一个有希望的1-序列，启动递归挖掘
        for(Integer item : promisingOneSequences){
            // 检查该 1-sequence 的投影是否存在
            if (!mapItemUC.containsKey(item)) {
                continue;
            }

            // 创建 1-sequence 模式
            ArrayList<ArrayList<Integer>> pattern = new ArrayList<>();
            ArrayList<Integer> itemset = new ArrayList<>();
            itemset.add(item);
            pattern.add(itemset);

            // 获取该 1-sequence 的投影链
            ArrayList<TargetedList> initialChain = mapItemUC.get(item);

            // 获取该 1-sequence 的真实非重叠效用 (u_no)
            // (我们已经在"阶段二"计算了 1-sequence 的真实效用 mapItemToUtility)
            long u_no_1_seq = mapItemToUtility.getOrDefault(item, 0L);

            // 检查 1-sequence 本身是否为解
            if (u_no_1_seq >= this.minUtility && patternContainsTarget(pattern, this.targetseq)) {
                foundPatterns.add(new ResultPattern(pattern, u_no_1_seq));
                patternCount++;
            }

            // 启用递归，传入模式和它的投影链
            patternGrowth(pattern, initialChain);
        }

        // 记录结束时间并计算总运行时间
        endTimestamp = System.currentTimeMillis();
        double runtime = endTimestamp - starTimestamp;
        // ======================= 第四阶段：格式化输出并写入文件 =========================
        StringBuilder builder = new StringBuilder();
        builder.append(this.datasetName).append("\n");
        builder.append("数据库总效用:").append(totalUtility).append("\n");
        builder.append("阈值:").append(utilityratio).append("\n");
        builder.append("效用阈值:").append(this.minUtility).append("\n");
        builder.append("目标序列序列:").append(targetseq).append("\n");
        builder.append("长度约束[minlen,maxlen]:").append(this.minlen).append("~").append(this.maxlen).append("\n");
        builder.append("间隔约束[mingap,maxgap]:").append(this.mingap).append("~").append(this.maxgap).append("\n");
        builder.append("高效用序列有: ").append(patternCount).append("个\n");
        builder.append("输出高效用序列的数量: ").append(patternCount).append("\n\n");
        builder.append("运行时间为: ").append(runtime / 1000.0).append(" s\n");
        builder.append("最大内存使用量: ").append(MemoryLogger.getInstance().getMaxMemory()).append(" MB\n");
        builder.append("候选集数量为: ").append(NumOfCandidate).append("\n");

        builder.append("输出高效用序列:\n");

        for (ResultPattern res : foundPatterns) {
            builder.append("序列: ").append(res.pattern.toString())
                    .append(" 效用: ").append(res.utility).append("\n");
        }

        builder.append("==打印完毕==");

// 将构建好的字符串一次性写入文件
        writer.write(builder.toString());

// 关闭writer
        writer.close();

        // 关闭qmatrix调试文件writer
        if(qmatrixWriter != null){
            qmatrixWriter.close();
        }

// 调用原有的打印统计函数，使其仍在控制台输出一份信息
        printStatistics(runtime);
    }

    /**
     * [核心递归函数 - TUSQ 投影 + HUNUS Nettree 融合版]
     *
     * @param pattern       当前要处理的候选模式
     * @param utilitychain  父模式的投影数据库 (Targeted Chain)
     */
    private void patternGrowth(ArrayList<ArrayList<Integer>> pattern,
                               ArrayList<TargetedList> utilitychain)
    {
        //获取父模式的最后一个项，用于I-Extension的有序性检查
        ArrayList<Integer> lastItemsetOfPattern = pattern.getLast();
        int lastItemOfPattern = lastItemsetOfPattern.getLast();
        // === 1. 宽度剪枝：查找扩展并计算 SEU ===

        // [修改] Key 不再是 Integer，而是 ExtensionKey (item + type)
        // [修改] Value 不再是 ExtensionTDU，而是 Long (SEU)
        Map<ExtensionKey, Long> extensionsMap = new HashMap<>();

        // 遍历父模式投影链中的 *每个序列*
        for (TargetedList utilitylist : utilitychain) {
            if (utilitylist.List.isEmpty()) continue;

            int sid = utilitylist.get_sid();
            QMatrix qm = this.database.get(sid);
            long localSEU = utilitylist.SEU; // SEU 是基于父模式的 SEU 计算的
            int prel = utilitylist.get_prel();

            // [新修改] 使用 Set 跟踪当前序列中已处理过的扩展项，避免重复累加 SEU
//            Set<Integer> processedExtensionsInSeq = new HashSet<>();

            // --- 1a. 查找 I-Extensions ---
            // 遍历父模式在该序列的 *每个实例*
            for (int j = 0; j < utilitylist.LengthOfUtilityList; j++) {
                int column = utilitylist.List.get(j).tid; // 父实例结束项集

                // 查找 QMatrix 中可用于 i-extension 的项
                for (int item : qm.getPresentItemsInItemset(column)) {
                    //确保只扩展比父模式最后一个项更大的项，以保证有序性
                    if (item > lastItemOfPattern) {
                        // 检查后缀约束
                        int newprel = updateprel(item, prel, 0); // 0 for i-ext
                        int firstSetPos = getFirstSetPos(newprel, this.targetSequence);//获取后缀的第一个项集在 T 中的位置（索引）
                        if (firstSetPos == 99999 || this.LI_T.get(sid).get(firstSetPos) > column ||
                                (this.LI_T.get(sid).get(firstSetPos) == column &&
                                        isIextFeasibleForNextTargetSet_I(item, prel))) {
//                            ExtensionTDU ext = extensionsMap.getOrDefault(item, new ExtensionTDU(0, "s"));
//                            ext.update(localSEU, "i"); // SEU 累加父模式的 SEU
//                            extensionsMap.put(item, ext);
//                            // 标记该 item 在此序列中已处理
//                            processedExtensionsInSeq.add(item);
                            // [修改] 使用 ExtensionKey
                            ExtensionKey key = new ExtensionKey(item, "i");
                            // [修改] 累加 SEU
                            extensionsMap.put(key, extensionsMap.getOrDefault(key, 0L) + localSEU);
                        }
                    }
                }
            }

            // --- 1b. 查找 S-Extensions ---
            // S-Extension 只需考虑父模式的 *第一个* 实例（TUSQ 优化）
            // [安全保护] 检查 List 是否为空，虽然在循环开始时已检查
            if (utilitylist.List.isEmpty()) continue;
            int startColumn = utilitylist.List.getFirst().tid + 1;

            // 查找 QMatrix 中可用于 s-extension 的项
            for (int item : qm.getItemsAfter(startColumn)) {
                // 检查后缀约束
                int newprel = updateprel(item, prel, 1); // 1 for s-ext
                // 找到该 item 在此序列中 *首次* 出现的位置 >= startColumn
                int itemFirstPos = qm.getFirstOccurrence(item, startColumn);

                if (itemFirstPos != -1) {
                    //获取后缀的第一个项集在 T 中的位置（索引）。
                    int firstSetPos = getFirstSetPos(newprel, this.targetSequence);
//                    if ((firstSetPos == 99999 || this.LI_T.get(sid).get(firstSetPos) >= itemFirstPos)
//                            && isSextFeasibleForNextTargetSet(item, firstSetPos))
                    if(firstSetPos == 99999 || this.LI_T.get(sid).get(firstSetPos) > itemFirstPos ||
                            (this.LI_T.get(sid).get(firstSetPos) == itemFirstPos &&
                                    isSextFeasibleForNextTargetSet_S(item, firstSetPos)))
                    {
//                        ExtensionTDU ext = extensionsMap.getOrDefault(item, new ExtensionTDU(0, "s"));
//                        ext.update(localSEU, "s"); // SEU 累加父模式的 SEU
//                        extensionsMap.put(item, ext);
                        // [修改] 使用 ExtensionKey
                        ExtensionKey key = new ExtensionKey(item, "s");
                        // [修改] 累加 SEU
                        extensionsMap.put(key, extensionsMap.getOrDefault(key, 0L) + localSEU);
                    }
                }
            }
        }

        // === 2. 递归：处理有希望的扩展 ===

        // 遍历所有找到的扩展
        for (Map.Entry<ExtensionKey, Long> entry : extensionsMap.entrySet()) {
            int item = entry.getKey().item;
            String type = entry.getKey().type;
            long seu = entry.getValue();

            // --- 宽度剪枝 (Width Pruning) ---
            if (seu < this.minUtility) {
                continue;
            }

            NumOfCandidate++;
            // a. 生成子模式
            ArrayList<ArrayList<Integer>> childPattern = extendPattern(pattern, item, type);

            // [新优化] 长度剪枝：如果子模式的项集数量已经超过maxlen，则直接跳过
            // 这是因为任何基于此模式的进一步扩展都会更长，同样不满足约束。
            if (childPattern.size() > this.maxlen) {
                continue;
            }

            // b. [核心] 构建子模式的投影 (Targeted Chain) 和 SEU
            Map.Entry<ArrayList<TargetedList>, Long> projection = constructProjection(
                    item, type, utilitychain);

            ArrayList<TargetedList> childChain = projection.getKey();
            long childSEU = projection.getValue();

            // 如果投影为空，则无需继续
            if (childChain.isEmpty()) {
                continue;
            }

            // c. [新优化] 先进行廉价的检查
            boolean containsTarget = patternContainsTarget(childPattern, this.targetseq);

            // 只有当模式可能成为最终结果时，才计算其真实效用
            if (containsTarget) {
                long u_no_valid = calculateNonOverlappingUtility_Projected(childPattern, childChain);

                // 检查是否输出
                if (u_no_valid >= this.minUtility) {
                    foundPatterns.add(new ResultPattern(childPattern, u_no_valid));
                    patternCount++;
                }
            }

            // e. [深度剪枝 (Depth Pruning)]
            // 只要子模式的效用上界有希望，就继续递归
            if (childSEU >= this.minUtility) {
                // 递归调用
                System.out.println("递归模式: " + childPattern + " SEU: " + childSEU);
                // [新修改] 不再传递效用参数
                patternGrowth(childPattern, childChain);
            }

            // [内存优化]
            // 在对一个子模式的递归探索完成后，我们手动清空其投影链并解除引用。
            // 这是为了明确地告诉垃圾回收器(GC)，这块巨大的内存可以被回收了。
            // 如果不这样做，这部分内存会一直被占用，直到当前`patternGrowth`函数的所有循环结束，
            // 在处理分支众多的搜索树时，这很容易导致内存溢出。
            childChain.clear();
        }
        // [内存优化]
        // 同样，在函数退出前，清空本次递归层级使用的扩展项map
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

    /**
     * [TUSQ 辅助方法] 更新前缀长度 (prel)。
     * 此版本使用 AlgoHUNUS 中的 int[] targetSequence 变量。
     *
     * @param exitem 扩展项
     * @param prel 当前 prel
     * @param kind 0:i-Concatenation, 1:s-Concatenation
     * @return 新的 prel
     */
    public int updateprel(int exitem, int prel, int kind){
        // TUSQ 的 T 数组以 -1, -2 结尾
        if(prel >= this.targetSequence.length - 2) return prel; // 已经匹配完成

        if(kind==0){ // i-extension
            if(this.targetSequence[prel]==-1) return prel; // 目标需要 s-extension，但当前是 i-extension
            else {
                // TUSQ 的逻辑是基于字典序比较
                if (this.targetSequence[prel] < exitem) { // 扩展项更大，此路不通，回溯 prel
                    while (prel > 0 && this.targetSequence[prel - 1] != -1) {
                        prel--;
                    }
                    return prel;
                }
                else if (this.targetSequence[prel] > exitem) return prel; // 扩展项更小，prel 不变
                else return prel + 1; // 匹配成功，prel+1
            }
        }
        else{ // s-extension
            if(this.targetSequence[prel]==-1){ // 目标需要 s-extension，匹配成功
                prel++;
                // 递归检查 s-extension 后的新项集是否匹配 (TUSQ 逻辑)
                return updateprel(exitem, prel, 0);
            }
            else{ // 目标需要 i-extension，但当前是 s-extension，此路不通，回溯 prel
                while(prel>0 && this.targetSequence[prel-1]!=-1) { prel--; }
                return updateprel(exitem, prel, 0);
            }
        }
    }

    /**
     * 获取后缀的第一个项集在 T 中的位置（索引）。
     *
     * @param prel 当前 prel
     * @param T 目标序列 (TUSQ 的 int[] 格式)
     * @return 项集索引。如果后缀为空，返回 99999。
     */
    public int getFirstSetPos(int prel, int []T){
        // T 数组以 -1, -2 结尾
        if(prel >= T.length - 2) return 99999;
        else{
            int curitemset = 0;
            for(int i=0; i<=prel; i++){
                if(T[i]==-1) curitemset++;
            }
            return curitemset;
        }
    }

    /**
     * 更严格的 s-拓展可行性判断：
     * 检查通过 s-拓展添加 item 后，新项集是否有可能通过后续的 i-拓展完整匹配目标项集。
     *
     * 核心约束：由于项集内的字典序（仅允许递增的 i-拓展），如果目标项集中存在比 item 更小的项，
     * 那么这些项将无法在之后通过 i-拓展补上（被顺序约束阻止），因而当前 s-拓展没有希望导向完整的目标后缀。
     *
     * 直观解释：
     * - 如果 item 在目标项集中：必须保证目标项集中不存在比 item 更小的项
     * - 如果 item 不在目标项集中：必须保证目标项集中不存在比 item 更小的项（否则无法通过 i-拓展添加这些更小的项）
     *
     * @param item         当前尝试的 s-拓展项
     * @param firstSetPos  新后缀的第一个目标项集在 targetseq 中的索引（若为 99999 表示后缀为空）
     * @return true 表示可行；false 表示应剪枝
     */
    private boolean isSextFeasibleForNextTargetSet_S(int item, int firstSetPos){
        if(firstSetPos == 99999){
            return true; // 后缀为空，无约束
        }
        // 与 getFirstSetPos 对齐：firstSetPos 对应 targetseq 的项集索引
        ArrayList<Integer> targetSet = this.targetseq.get(firstSetPos);

        // 关键修复：无论 item 是否在目标项集中，都要检查目标项集中是否存在比 item 更小的项
        // 如果存在，则无法通过 i-拓展添加这些更小的项（因为字典序约束），因此应该剪枝
        for(Integer y : targetSet){
            if(y < item){
                // 目标项集中存在比 item 更小的项，无法通过 i-拓展补齐，应剪枝
                return false;
            }
        }
        // 目标项集中所有项都 >= item，可以通过 i-拓展逐步添加，可行
        return true;
    }

    /**
     * 检查一个 i-extension 在必须于当前项集完成下一个目标匹配时是否可行。
     * 这个检查是针对一个非常特殊且关键的剪枝场景：
     * 当我们知道下一个目标项集必须在当前数据库项集（column）内完成匹配时，
     * 我们进行的任何 I-Extension 都不能“越过”下一个需要匹配的目标项。
     *
     * @param itemToExtend 要用于扩展的项。
     * @param prel         当前模式已经匹配到的目标序列 T 的前缀长度。因此，下一个要匹配的目标项就是 T[prel]。
     * @return 如果扩展项小于下一个目标项，则返回 true，表示可行；否则返回 false。
     */
    private boolean isIextFeasibleForNextTargetSet_I(int itemToExtend, int prel) {
        // 如果 prel 已经到达 T 的末尾（由-2标记），说明目标已完全匹配，任何扩展都是可行的。
        if (prel >= this.targetSequence.length - 2) {
            return true;
        }

        // 获取下一个需要匹配的目标项。
        int nextTargetItem = this.targetSequence[prel];

        // 如果下一个目标“项”实际上是一个项集分隔符（-1），
        // 这意味着目标序列期望一个 S-Extension，而我们正在进行 I-Extension。
        // 在这种情况下，字典序检查不适用，我们认为它是可行的（由其他逻辑处理）。
        if (nextTargetItem == -1) {
            return true;
        }

        // 核心逻辑：扩展项必须小于下一个目标项。
        // 因为在同一个项集内，项是按字典序处理的，如果扩展了一个更大的项，就无法再回头匹配那个更小的目标项了。
        return itemToExtend <= nextTargetItem;
    }

    /**
     * 构建投影数据库（Targeted Chain）。
     * 此方法根据父模式的效用链，构建子模式的效用链，并计算子模式的全局 SEU。
     *
     * @param item 扩展项
     * @param type 扩展类型 ("i" 或 "s")
     * @param utilitychain 父模式的效用链
     * @return 子模式的投影信息 (新的效用链 和 全局SEU)
     */
    private AbstractMap.SimpleEntry<ArrayList<TargetedList>, Long> constructProjection(
            int item, String type,
            ArrayList<TargetedList> utilitychain)
    {
        ArrayList<TargetedList> newUChain = new ArrayList<>(); // 子模式的新效用链
        long globalSEU = 0; // 子模式的全局 SEU

        // 遍历父模式效用链中的每个序列
        for (TargetedList utilitylist : utilitychain) {
            int sid = utilitylist.get_sid();
            QMatrix qm = this.database.get(sid);

            // 查找扩展项 item 在 QMatrix 中的行索引
            int row = qm.getItemIndex(item);
            if (row == -1) continue; // 该序列不包含扩展项，跳过

            TargetedList newUL = new TargetedList(); // 子模式在该序列的新效用链
            newUL.set_sid(sid);

            long localSEU = 0; // 子模式在该序列的 SEU

            if ("i".equals(type)) {
                // --- I-Extension ---
                // 遍历父模式在该序列的每个实例
                for (int j = 0; j < utilitylist.LengthOfUtilityList; j++) {
                    TargetedList.UtilityElement ue = utilitylist.List.get(j);
                    int column = ue.tid; // 父实例的结束项集索引

                    // 检查 QMatrix，看扩展项 item 是否在 *同一* 项集
                    int itemUtility = qm.getUtility(column, item);

                    if (itemUtility != 0) { // 找到了 i-extension
                        int prel = utilitylist.get_prel();
                        int newprel = updateprel(item, prel, 0); // 0 for i-ext

                        // 检查后缀约束
                        int firstSetPos = getFirstSetPos(newprel, this.targetSequence);
                        if (firstSetPos == 99999 || this.LI_T.get(sid).get(firstSetPos) > column ||
                                (this.LI_T.get(sid).get(firstSetPos) == column && isIextFeasibleForNextTargetSet_I(item, prel))) {
                            newUL.set_prel(newprel);
                            int acu = ue.acu + itemUtility;
                            int ru = qm.getRemainingUtility(column, item); // 使用原始RU
                            newUL.add(column, acu, ru);

                            // 更新 localSEU
                            if (acu + ru > localSEU)
                                localSEU = acu + ru;
                        }else {
//                            System.out.println("剪枝 i-ext: item=" + item + " sid=" + sid + " column=" + column +
//                                    " firstSetPos=" + firstSetPos + " LI_T=" + this.LI_T.get(sid).get(firstSetPos));
                        }
                    }
                }
            } else {
                // --- S-Extension ---
                // 遍历父模式在该序列的每个实例
//                for (TargetedList.UtilityElement ue : utilitylist.List) {
//                    int startColumn = ue.tid + 1; // S-extension 必须在父实例之后
//
//                    // 遍历 QMatrix 中该实例之后的所有项集
//                    for (int k = startColumn; k < qm.getNbItemsets(); k++) {
//                        int itemUtility = qm.getUtility(k, item);
//
//                        if (itemUtility != 0) { // 找到了 s-extension
//                            int prel = utilitylist.get_prel();
//                            int newprel = updateprel(item, prel, 1); // 1 for s-ext
//
//                            int firstSetPos = getFirstSetPos(newprel, this.targetSequence);
//                            // 优化：当 LI_T > k 时，说明有足够空间，不需要可行性检查
//                            // 只有当 LI_T == k 时（边界情况），才需要进行可行性检查
//                            if (firstSetPos == 99999 || this.LI_T.get(sid).get(firstSetPos) > k ||
//                                    (this.LI_T.get(sid).get(firstSetPos) == k && isSextFeasibleForNextTargetSet_S(item, firstSetPos))) {
//                                newUL.set_prel(newprel);
//                                int acu = ue.acu + itemUtility;
//                                int ru = qm.getRemainingUtility(k, item); // 使用原始RU
//                                newUL.add(k, acu, ru);
//
//                                if (acu + ru > localSEU) localSEU = acu + ru;
//                            }else {
//                                System.out.println("剪枝 s-ext: item=" + item + " sid=" + sid + " k=" + k +
//                                        " firstSetPos=" + firstSetPos + " LI_T=" + this.LI_T.get(sid).get(firstSetPos));
//                            }
//                        }
//                    }
//                }
                if(!utilitylist.List.isEmpty()){
                    TargetedList.UtilityElement firstUe = utilitylist.List.getFirst();
                    int startColumn = firstUe.tid + 1; // S-extension 必须在父实例之后

                    // 遍历 QMatrix 中该实例之后的所有项集
                    for (int k = startColumn; k < qm.getNbItemsets(); k++) {
                        int itemUtility = qm.getUtility(k, item);

                        if (itemUtility != 0) { // 找到了 s-extension
                            int prel = utilitylist.get_prel();
                            int newprel = updateprel(item, prel, 1); // 1 for s-ext

                            int firstSetPos = getFirstSetPos(newprel, this.targetSequence);
                            // 优化：当 LI_T > k 时，说明有足够空间，不需要可行性检查
                            // 只有当 LI_T == k 时（边界情况），才需要进行可行性检查
                            if (firstSetPos == 99999 || this.LI_T.get(sid).get(firstSetPos) > k ||
                                    (this.LI_T.get(sid).get(firstSetPos) == k && isSextFeasibleForNextTargetSet_S(item, firstSetPos))) {
                                newUL.set_prel(newprel);
                                int acu = firstUe.acu + itemUtility;
                                int ru = qm.getRemainingUtility(k, item); // 使用原始RU
                                newUL.add(k, acu, ru);

                                if (acu + ru > localSEU) localSEU = acu + ru;
                            }
                        }
                    }
                }


            }

            // 如果子模式在该序列中找到了实例
            if (newUL.LengthOfUtilityList > 0) {
                newUL.set_SEU((int) localSEU);
                newUChain.add(newUL);
                globalSEU += localSEU; // 累加全局 SRU
            }
        }

        // 使用 Java 8 Pair 来返回两个值。如果环境不支持，可以创建一个简单的 Pair 类。
//         return new Pair<>(newUChain, globalSEU);
        // 为确保兼容性，我们使用 Map.Entry
        return new AbstractMap.SimpleEntry<>(newUChain, globalSEU);
    }



//    /**
//     * [HUNUS 融合方法] 计算投影链的真实非重叠效用。
//     * 此方法遍历投影链，仅对包含该模式的序列运行 Nettree 算法。
//     *
//     * @param pattern      要计算效用的模式
//     * @param utilitychain 该模式的投影效用链
//     * @return 真实的、满足约束的非重叠效用 (u_no_valid)
//     */
//    private long calculateNonOverlappingUtility_Projected(
//            ArrayList<ArrayList<Integer>> pattern,
//            ArrayList<TargetedList> utilitychain)
//    {
//        long u_no_valid = 0;
//
//        // 仅遍历投影链中的序列，而不是整个数据库
//        for (TargetedList ul : utilitychain) {
//            int sid = ul.get_sid();
//            QMatrix qm = this.database.get(sid);
//            ArrayList<ArrayList<Integer>> sequence = qm.getSequenceAsItemsets();
//
//            // 1. 构建 Nettree
//            // 注意：Nettree 结构应为实例变量或重用，以避免重复创建
//            // 为简化，这里暂时重新创建
//            ArrayList<Node>[] nettree = (ArrayList<Node>[]) new ArrayList[pattern.size()];
//            createNetTree(nettree, sequence, pattern);
//
//            // [新增步骤] 更新 Nettree, 计算 min/max leave, 为剪枝做准备
//            updateNetTreeWithLeaves(nettree, pattern.size());
//
//            // 2. 查找非重叠实例
//            ArrayList<Instance> nonOverlappingInSeq = findNonOverlappingInSequence(nettree, pattern, sid, qm);
//
//            // 3. 累加满足约束的效用
//            for (Instance inst : nonOverlappingInSeq) {
//                // 实例是在findNonOverlappingInSequence中找到的，并且已经通过了minlen/maxlen检查
//                u_no_valid += inst.utility;
//            }
//        }
//        return u_no_valid;
//    }

    /**
     * [HUNUS 融合方法] 计算投影链的真实非重叠效用。
     * 此方法遍历投影链，仅对包含该模式的序列运行 Nettree 算法。
     * (已修改为重用 Nettree)
     *
     * @param pattern      要计算效用的模式
     * @param utilitychain 该模式的投影效用链
     * @return 真实的、满足约束的非重叠效用 (u_no_valid)
     */
    private long calculateNonOverlappingUtility_Projected(
            ArrayList<ArrayList<Integer>> pattern,
            ArrayList<TargetedList> utilitychain)
    {
        long u_no_valid = 0;
        int patternSize = pattern.size(); // 获取模式大小

        // *** 1. (新增) 检查并确保可重用 Nettree 足够大 ***
        checkAndResizeNettree(patternSize);

        // 仅遍历投影链中的序列，而不是整个数据库
        for (TargetedList ul : utilitychain) {
            int sid = ul.get_sid();
            QMatrix qm = this.database.get(sid);
            ArrayList<ArrayList<Integer>> sequence = qm.getSequenceAsItemsets();

            // *** 2. (新增) 在处理每个序列之前，清理 Nettree ***
            cleanNettree(patternSize);

            // *** 3. (修改) 使用 this.reusableNettree 代替局部变量 ***

            // 1. 构建 Nettree (传入 reusableNettree)
            // [原代码] ArrayList<Node>[] nettree = (ArrayList<Node>[]) new ArrayList[pattern.size()];
            createNetTree(this.reusableNettree, sequence, pattern); // [修改]

            // [新增步骤] 更新 Nettree, 计算 min/max leave, 为剪枝做准备
            updateNetTreeWithLeaves(this.reusableNettree, patternSize); // [修改]

            // 2. 查找非重叠实例
            ArrayList<Instance> nonOverlappingInSeq = findNonOverlappingInSequence(this.reusableNettree, pattern, sid, qm); // [修改]

            // 3. 累加满足约束的效用
            for (Instance inst : nonOverlappingInSeq) {
                // 实例是在findNonOverlappingInSequence中找到的，并且已经通过了minlen/maxlen检查
                u_no_valid += inst.utility;
            }
        }
        return u_no_valid;
    }


    /**
     * (重载版本) 判断序列A是否包含子序列B
     * @param sequenceA 序列 A (标准格式)
     * @param subSequence 子序列 B (标准格式)
     * @return true 如果包含，否则 false
     */
    public boolean seqContain(ArrayList<ArrayList<Integer>> sequenceA, ArrayList<ArrayList<Integer>> subSequence) {
        if (subSequence.isEmpty()) {
            return true;
        }
        if (sequenceA.isEmpty()) {
            return false;
        }

        int subIdx = 0; // 指向子序列的指针
        int mainIdx = 0; // 指向主序列的指针

        while (subIdx < subSequence.size() && mainIdx < sequenceA.size()) {
            // 如果主序列的当前项集包含了子序列的当前项集
            if (sequenceA.get(mainIdx).containsAll(subSequence.get(subIdx))) {
                subIdx++; // 匹配成功，子序列指针前进
            }
            mainIdx++; // 无论是否匹配，主序列指针总是前进
        }

        // 如果子序列的指针走到了最后，说明所有部分都找到了匹配
        return subIdx == subSequence.size();
    }

    /**
     * transfer array A to a list
     * @param T target sequence
     * @return list T.
     */
    public ArrayList<ArrayList<Integer>> convertT(int[] T) {
        if (T == null) {
            throw new IllegalArgumentException("输入数组T不能为null");
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

    /**
     * 为一个给定的模式P，在序列S中构建Nettree
     * @param nettree 用于存储结果的Nettree结构
     * @param sequence S, 表示为一个项集的列表, e.g., List<List<Integer>>
     * @param pattern P, 同样是项集的列表, e.g., List<List<Integer>>
     */

    void createNetTree(ArrayList<Node>[] nettree, ArrayList<ArrayList<Integer>> sequence,
                       ArrayList<ArrayList<Integer>> pattern){
        int patternLayerCount = pattern.size();//targetseq的项集数，即Nettree的层数

        //1.初始化Nettree和start优化数组
//        for(int i = 0; i < patternLayerCount; i++){
//            if(nettree[i] == null){
//                nettree[i] = new ArrayList<>();
//            }else{
//                nettree[i].clear();
//            }
//        }
        int[] start = new int[patternLayerCount];//全部默认为0

        //2.主循环，遍历itemsets
        for(int itemsetIdx = 0; itemsetIdx < sequence.size(); itemsetIdx++){
            ArrayList<Integer> currentSeqItemset = sequence.get(itemsetIdx);

            //处理第0层（根节点）
            List<Integer> patternItem0 = pattern.getFirst();
            //3.匹配逻辑：包含
            if(currentSeqItemset.containsAll(patternItem0)){
                Node newNode = new Node();
                newNode.name = itemsetIdx;//4.Node.name存储的是"项集索引"
                newNode.level = 0;
                newNode.minLeave = itemsetIdx;
                newNode.maxLeave = itemsetIdx;
                nettree[0].add(newNode);
            }

            //处理第1层到最后一层
            //遍历模式的每个部分，看当前itemset能否作为该部分的结点
            for(int j = 1; j < patternLayerCount; j++){
                ArrayList<Integer> patterItem_j = pattern.get(j);

                if(currentSeqItemset.containsAll(patterItem_j)){
                    ArrayList<Node> parentLayer = nettree[j - 1];
                    if(parentLayer.isEmpty()){
                        continue;//如果父层是空的，不可能有父节点，跳过
                    }

                    //5.start数组
                    //更新start指针，跳过那些因为maxgap约束而不可能成为父节点的旧节点
                    while(start[j-1] < parentLayer.size() &&
                            (itemsetIdx - parentLayer.get(start[j-1]).name -1) > this.maxgap){
                        start[j-1]++;
                    }

                    //6.检查是否存在可能的父节点
                    //如果start指针已经越界，或者当前项集离第一个父节点都太近（不满足mingap），则不可能有连接
                    if(start[j-1] >= parentLayer.size() ||
                            (itemsetIdx - parentLayer.getFirst().name - 1  < this.mingap)){
                        continue;
                    }

                    //7.创建子节点并建立连接
                    //如果走到这里，说明至少有一个潜在的父节点，我们创建当前子节点
                    Node childNode = new Node();
                    childNode.name = itemsetIdx;
                    childNode.level = j;
                    childNode.minLeave = itemsetIdx;
                    childNode.maxLeave = itemsetIdx;
                    nettree[j].add(childNode);

                    //遍历所有可能的父节点（从start[j]开始），建立连接
                    for(int k = start[j-1]; k < parentLayer.size(); k++){
                        Node parentNode = parentLayer.get(k);
                        int currentGap = itemsetIdx - parentNode.name - 1;

                        //因为上面的检查，这里的gap肯定 <= maxgap
                        //如果gap太小，后续的父节点只会更小，直接跳出
                        if(currentGap < this.mingap){
                            break;
                        }

                        //满足所有约束，建立父子双向连接
                        parentNode.children.add(childNode);
                        childNode.parent.add(parentNode);
                    }
                }
            }
        }
    }

    /**
     * [新增方法]
     * 使用反向贪心搜索策略，在单个序列(sequence)中寻找目标(target)的最后一个实例路径。
     * 该方法逻辑源自 TUSQ 算法的 fillLOT，并进行了适配。
     * @param target 目标序列 T 的标准格式 (ArrayList<ArrayList<Integer>>)
     * @param sequence 当前数据库序列 S 的标准格式 (ArrayList<ArrayList<Integer>>)
     * @return 如果找到完整的最后实例，则返回一个包含项集索引的路径列表；否则返回 null。
     */
    private ArrayList<Integer> findLastInstancePathWithFillLIT(
            ArrayList<ArrayList<Integer>> target,
            ArrayList<ArrayList<Integer>> sequence) {
        // 该方法现在通过迭代寻找最后一个有效实例，而不仅仅是字面上的最后一个。
        // 从序列的末尾开始扫描
        int startScanIdx = sequence.size() - 1;

        while (startScanIdx >= 0) {
            // 1. 初始化一个数组，用于存储T的每个项集在S中最后出现的位置
            int[] lastpos = new int[target.size()];
            Arrays.fill(lastpos, -1); // 初始化为-1，表示尚未找到

            // 2. 核心逻辑：反向贪心搜索
            // 从目标序列T的最后一个项集开始匹配
            int targetIdx = target.size() - 1;

            // 从当前的 startScanIdx 开始在序列S中反向扫描
            for (int seqIdx = startScanIdx; seqIdx >= 0; seqIdx--) {
                if (targetIdx < 0) {
                    // T的所有项集都已匹配完
                    break;
                }

                // 如果S的当前项集 包含了 T的当前项集
                if (sequence.get(seqIdx).containsAll(target.get(targetIdx))) {
                    // 记录下这个位置
                    lastpos[targetIdx] = seqIdx;
                    // 继续向前匹配T的前一个项集
                    targetIdx--;
                }
            }

            // 3. 检查是否找到了完整的实例
            if (targetIdx >= 0) {
                // 如果targetIdx不为-1，说明T的前缀部分没有找到匹配。
                // 因为我们是从右向左扫描的，所以不可能再找到其他实例了。
                return null;
            }

            // 4. 找到了一个完整实例。现在，检查它是否满足mingap和maxlen约束。
            boolean isValid = true;
            // 4.1 检查maxlen约束
            int instanceLength = lastpos[lastpos.length - 1] - lastpos[0] + 1;
            if (instanceLength > this.maxlen) {
                isValid = false;
            }

            // 4.2 检查mingap约束
            if (isValid) {
                for (int i = 0; i < lastpos.length - 1; i++) {
                    int gap = lastpos[i + 1] - lastpos[i] - 1;
                    if (gap < this.mingap) {
                        isValid = false;
                        break;
                    }
                }
            }

            // 5. 如果实例有效，任务完成。返回路径。
            if (isValid) {
                ArrayList<Integer> path = new ArrayList<>();
                for (int pos : lastpos) {
                    path.add(pos);
                }
                return path;
            } else {
                // 如果实例无效，我们需要继续搜索一个更早的实例。
                // 下一次搜索应该从我们找到的这个无效实例的最后一个部分的前一个位置开始。
                startScanIdx = lastpos[lastpos.length - 1] - 1;
            }
        }

        // 如果while循环结束，说明在整个序列中没有找到有效的实例。
        return null;
    }

    /**
     * 专门用于寻找第一个有效实例的函数
     * @return 返回一个代表实例路径的整数列表(项集索引)，如果找不到则返回null
     */
    private ArrayList<Node> findFirstValidInstance(ArrayList<Node>[] nettree,
                                                   ArrayList<ArrayList<Integer>> pattern){
        if(nettree[0].isEmpty()){
            return null;
        }

        //从最左边的根节点开始，遍历
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

    /**
     * 递归辅助函数，寻找从currentNode出发并满足约束的、最靠左的一条完整路径
     * @param rootNode 起始根节点，用于计算最终的span
     * @param currentNode 当前正在访问的节点
     * @param currentLevel 当前所在的层级
     * @param targetLevel 目标层级（即模式的长度）
     * @param path 用于在递归中构建路径的列表
     * @return 如果找到有效路径，则返回该路径的Node列表；否则返回null
     */
    private ArrayList<Node> findFullPathRecursiveFromFirst(Node rootNode, Node currentNode,  int currentLevel,
                                                           int targetLevel, ArrayList<Node> path){
        //CHECK: Skip if current node is already used
        if(currentNode.used){
            return null;
        }

        //1.将当前节点加入路径（进入当前递归层）
        path.add(currentNode);

        //2.如果已经达到目标层数说明找到了一条完整路径
        if(currentLevel == targetLevel){
            int span = currentNode.name - rootNode.name + 1;
            if(span >= minlen && span <= maxlen){
                //找到一个有效的完整路径，返回它的一个副本
                ArrayList<Node> result = new ArrayList<>(path);
                //在返回前，执行回溯，将当前节点从路径中移除
                path.removeLast();
                return result;
            }else {
                //路径完整，但不满足len约束，同样需要回溯
                path.removeLast();
                return null;
            }
        }

        //3.递归步骤：从左向右遍历所有节点
        for(Node childNode : currentNode.children){
            // 在深入递归之前，立刻检查从根到这个子节点的span
            int currentSpan = childNode.name - rootNode.name + 1;
            if (currentSpan > this.maxlen) {
                // 如果路径到这里已经超长，那么这条分支以及它右边的所有兄弟分支
                // (因为它们的name更大) 都不可能满足maxlen约束了。
                // 因此我们可以直接终止对当前父节点的所有子节点的搜索。
                break; // 剪枝！
            }
            // *** 剪枝逻辑结束 ***

            //向下一层递归
            ArrayList<Node> foundPath = findFullPathRecursiveFromFirst(rootNode, childNode,
                    currentLevel + 1, targetLevel, path);
            if(foundPath != null){
                //如果成功了，说明我们已经找到了最右边的那条可行路径，无需再试其他兄弟节点
                //将结果向上传递之前，同样需要回溯
                path.removeLast();
                return foundPath;
            }
            //如果返回null,说明从这个子节点出发是死胡同，循环将继续，尝试左边的下一个兄弟节点
        }

        //4.回溯：如果所有子节点都无法形成完整路径，则将当前节点移出路径
        path.removeLast();

        return null;//从当前节点出发的所有分支都走不通
    }

    /**
     * 移除一个已经找到的实例路径，并级联清理因此产生的孤立节点。
     * @param pathNodes     要移除的、由Node对象组成的实例路径
     */
    void updateNetTreeToRemovePath(ArrayList<Node> pathNodes){
        if(pathNodes == null || pathNodes.isEmpty()){
            return;
        }

        //步骤一：直接将路径上的所有节点标记为"已使用"
        for(Node node : pathNodes){
            node.used = true;
        }

        //步骤二：反向传播"已使用"状态，清理因此产生的"死胡同"父节点
        //从路径的倒数第二个节点（即路径上的第一个父节点）开始，向前回溯到根
        for(int i = pathNodes.size() - 2; i >= 0; i--){
            Node parentNode = pathNodes.get(i);

            //如果父节点已经被（之前的传播）标记为used，则无需重复处理
            if(parentNode.used){
                continue;
            }

            //检查该父节点的所有子节点是否都已被使用
            boolean allChildernUsed = true;
            //如果父节点没有子节点了（不太可能但做保护），也视为"已用完"
            if(parentNode.children.isEmpty()){
                allChildernUsed = true;
            }else{
                for(Node child : parentNode.children){
                    if(!child.used){
                        //只要有一个子节点还没被使用，就说明父节点还有其他路径可走
                        allChildernUsed = false;
                        break;
                    }
                }
            }

            //如果所有子节点确实都用完了，那么这个父节点本身也变成了死胡同
            if(allChildernUsed){
                parentNode.used = true;
            }
        }


    }

    /**
     * [新辅助函数] 对于单个序列，基于已构建的Nettree，找出所有非重叠实例。
     * @param nettree       为当前序列和模式构建好的Nettree
     * @param pattern       当前模式
     * @param sid           当前序列的ID
     * @param qm            当前序列的QMatrix
     * @return              该序列中所有非重叠实例的列表
     */
    private ArrayList<Instance> findNonOverlappingInSequence(
            ArrayList<Node>[] nettree,
            ArrayList<ArrayList<Integer>> pattern,
            int sid,
            QMatrix qm) {
        ArrayList<Instance> instancesInSeq = new ArrayList<>();

        // [新增步骤] 预处理和剪枝根节点
        // 在进入循环之前，我们可以先过滤掉那些不可能形成有效实例的根节点
        if (nettree.length > 0 && nettree[0] != null) {
            for (Node rootNode : nettree[0]) {
                if (!rootNode.toleave) { // 如果从该根节点出发，无法到达最终层，则标记为已使用
                    rootNode.used = true;
                    continue;
                }
                // 检查长度约束
                int shortestPossibleSpan = rootNode.minLeave - rootNode.name + 1;
                int longestPossibleSpan = rootNode.maxLeave - rootNode.name + 1;

                if (longestPossibleSpan < this.minlen || shortestPossibleSpan > this.maxlen) {
                    // 如果从该根节点出发的最长路径都达不到minlen，
                    // 或者最短路径都已经超过maxlen，那么它是不可能满足约束的
                    rootNode.used = true;
                }
            }
        }

        // 循环：不断从nettree中寻找第一个有效实例，然后"销毁"它，再找下一个
        while (true) {
            // 1. 寻找当前Nettree状态下的第一个（最左）有效实例路径
            ArrayList<Node> nodepath = findFirstValidInstance(nettree, pattern);

            // 2. 如果再也找不到实例了，就跳出循环
            if (nodepath == null) {
                break;
            }

            // 3. 计算这个实例的效用
            long instanceUtility = calculateUtilityOfInstance(nodepath, qm, pattern);

            // 4. 将Node路径转换为项集索引路径，并创建Instance对象
            ArrayList<Integer> indexPath = new ArrayList<>();
            for (Node node : nodepath) {
                indexPath.add(node.name);
            }
            instancesInSeq.add(new Instance(sid, indexPath, instanceUtility));

            // 5. "销毁"这个实例路径上的节点及其关联，为寻找下一个非重叠实例做准备
            updateNetTreeToRemovePath(nodepath);
        }
        return instancesInSeq;
    }

    //填充LI_T表
    public void fillLITableForSequence(ArrayList<ArrayList<Integer>> LI_T,int NumberOfSequence,
                                       ArrayList<Integer> lastInstancePath){
        //确保LI_T有足够行
        while(LI_T.size() <= NumberOfSequence){
            LI_T.add(new ArrayList<>());
        }

        // 清空当前行（可选，防止重复添加）
        LI_T.get(NumberOfSequence).clear();

        LI_T.get(NumberOfSequence).addAll(lastInstancePath);
    }

//    /**
//     * 检查一个1-序列的出现是否"有希望"
//     * (此版本专门用于处理1-序列的SEU计算)
//     * @param item 当前1-序列的项
//     * @param itemsetIdx 当前项所在的项集索引
//     * @param sid 当前序列的ID
//     * @param target 目标序列T
//     * @param liTable LI-T
//     * @return 是否有希望
//     */
//    //true表示该item有潜力拓展
//    //检查一个1-序列的出现是否"有希望"
//    private boolean checkIsPromising(int item, int itemsetIdx, int sid,
//                                     ArrayList<ArrayList<Integer>> target,
//                                     ArrayList<ArrayList<Integer>> liTable){
//        // 检查 liTable 是否已初始化或 sid 是否合法
//        if (sid >= liTable.size() || liTable.get(sid) == null || liTable.get(sid).isEmpty()) {
//            System.out.println("sid:" + sid + "不合法");
//
//            return false; // 或无数据时保守返回 true
//        }
//
//        //确定后缀开始的位置（在target中的索引）
//        int suffixStartIndexInTarget;
//
//        if(!target.getFirst().contains(item)){
//            suffixStartIndexInTarget = 0;
//        }else{
//            if(target.size() > 1){
//                suffixStartIndexInTarget = 1;
//            }else{
//                //若目标序列只有一个项集
//                return true;
//            }
//        }
//
//        //查询LI_T
//        //LI_T的行索引sid，列索引是后缀开始位置在T中的索引-1（因为LI_T的列也是从0开始）
//        //无论如何都和第一个目标序列的第一个项集比较位置
//        int lastPosOfSuffixStart = liTable.get(sid).get(suffixStartIndexInTarget);
//
//        //比较位置
//        return lastPosOfSuffixStart > itemsetIdx;
//    }

//    /**
//     * 通用版本：检查一个候选模式 P 的实例是否 "有希望" (Promising)。
//     * @param instance P 的一个具体实例
//     * @param sid 当前序列的ID
//     * @return 是否有希望
//     */
//    private boolean checkIsPromising(Instance instance, int sid){
//        if(sid >= this.LI_T.size() || this.LI_T.get(sid).isEmpty()){
//            return false;
//        }
//
//        //1.找到 P 包含了 T 的多长的前缀
//        int prefixLen = findLongestPrefixContained(this.targetseq, this.targetseq);
//
//        //2.如果 P 已经包含了整个 T ，那么任何拓展都是有希望的
//        if(prefixLen == this.targetseq.size()){
//            return true;
//        }
//
//        //3. T 的后缀(Suf(T, P)从 targetseq 的第 prefixLen 个项集开始
//        //4. 从LI_T 中查询后缀头部的最后一个出现位置
//        int lastPosOfSuffixStart = this.LI_T.get(sid).get(prefixLen);
//
//        //5.比较位置：后缀的最后一次出现必须在当前实例结束之后，或者在同一个项集开始。
//        return lastPosOfSuffixStart >= instance.getEndItemsetIdx();
//    }



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

    /**
     * 计算一个模式实例的总效用值
     * @param nodepath 实例路径，包含每个项集的索引位置
    //     * @param sequence 序列数据，以项集列表形式表示
     * @param qm 对应的QMatrix对象
     * @return 该实例的总效用值
     */
    private long calculateUtilityOfInstance(ArrayList<Node> nodepath,
                                            QMatrix qm,
                                            ArrayList<ArrayList<Integer>> pattern){
        long totalUtility = 0;

        //检查路径长度是否与模式长度一致
        if(nodepath.size() != pattern.size()){
            throw new IllegalArgumentException("path and pattern size do not match");
        }

        //遍历路径中的每个项集位置和模式中的对应项集
        for(int i = 0; i < nodepath.size(); i++){
            int itemIdx = nodepath.get(i).name;//实例中当前项集的索引
            ArrayList<Integer> patternItemset = pattern.get(i);//模式中对应的项集

            //计算当前项集中模式对应项的效用
            for(int patternItem : patternItemset){
                totalUtility += qm.getUtility(itemIdx, patternItem);
            }
        }
        return totalUtility;
    }

    /**
     * [新增方法] 穷尽地、递归地寻找一个模式在Nettree中的所有合法路径。
     * @param rootNode      当前路径的起始（根）节点
     * @param currentNode   当前递归正在访问的节点
     * @param currentLevel  当前节点的层级
     * @param targetLevel   模式的总层数（递归终点）
     * @param allPaths      用于收集所有找到的合法路径的列表
     * @param currentPath   当前正在构建的路径
     */
    private void findAllPathsRecursive(Node rootNode, Node currentNode, int currentLevel, int targetLevel,
                                       List<ArrayList<Node>> allPaths, ArrayList<Node> currentPath) {
        // 将当前节点加入正在构建的路径
        currentPath.add(currentNode);

        // 如果到达了目标层级，说明一条完整的路径找到了
        if (currentLevel == targetLevel) {
            // 检查路径是否满足maxlen约束
            if (currentNode.name - rootNode.name + 1 <= this.maxlen) {
                // 将这条合法路径的副本加入结果集
                allPaths.add(new ArrayList<>(currentPath));
            }
        } else {
            // 如果还没到最后一层，则继续向子节点递归
            for (Node childNode : currentNode.children) {
                // 剪枝：如果路径到下一个节点已经超长，则跳过该分支
                if (childNode.name - rootNode.name + 1 <= this.maxlen) {
                    findAllPathsRecursive(rootNode, childNode, currentLevel + 1, targetLevel, allPaths, currentPath);
                }
            }
        }

        // 回溯：在返回上一层之前，将当前节点从路径中移除
        currentPath.remove(currentPath.size() - 1);
    }

    /**
     * 找到目标序列 target 在 候选模式 pattern 中匹配上的最长前缀的长度。
     * @param target 目标序列 T
     * @param pattern 候选模式 P
     * @return 匹配上的 T 的前缀的项集数量 (长度)
     */
    private int findLongestPrefixContained(ArrayList<ArrayList<Integer>> target,
                                           ArrayList<ArrayList<Integer>> pattern){
        int targetIdx = 0;
        int patternIdx = 0;
        while(targetIdx < target.size() && patternIdx < pattern.size()){
            //判断模式项集是否包含目标项集
            if(pattern.get(patternIdx).containsAll(target.get(targetIdx))){
                targetIdx++;
            }
            patternIdx++;
        }
        return targetIdx;//返回匹配上的目标前缀长度
    }

    /**
     * 根据给定的模式、拓展项和拓展类型，生成新的子模式。
     * @param parentPattern 父模式
     * @param item 拓展项
     * @param type 拓展类型 ("i" 或 "s")
     * @return 新生成的子模式
     */
    private ArrayList<ArrayList<Integer>> extendPattern(ArrayList<ArrayList<Integer>> parentPattern,
                                                        int item, String type){
        //深拷贝父模式，避免修改原始对象
        ArrayList<ArrayList<Integer>> newPattern = new ArrayList<>();
        for(ArrayList<Integer> itemset : parentPattern){
            newPattern.add(new ArrayList<>(itemset));
        }

        if("i".equals(type)){
            //I-拓展：将item加入最后一个项集并排序
            ArrayList<Integer> lastItemset = newPattern.getLast();
            lastItemset.add(item);
            Collections.sort(lastItemset);
        }else{
            //S-拓展：创建一个包含item的新项集并追加到模式结尾
            ArrayList<Integer> newItemset = new ArrayList<>();
            newItemset.add(item);
            newPattern.add(newItemset);
        }
        // 移除了重复的 NumOfCandidate++，因为在 patternGrowth 中已经增加过了
        return newPattern;
    }

    /**
     * 通用版本：判断模式pattern是否包含目标target (作为子序列)
     * @param pattern 候选模式
     * @param target 目标序列
     * @return true如果包含，否则false
     */
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

    // 用于存储最终解的内部类
    static class ResultPattern {
        ArrayList<ArrayList<Integer>> pattern;
        long utility;

        public ResultPattern(ArrayList<ArrayList<Integer>> pattern, long utility) {
            // 深拷贝，防止后续的递归修改pattern内容
            this.pattern = new ArrayList<>();
            for (ArrayList<Integer> itemset : pattern) {
                this.pattern.add(new ArrayList<>(itemset));
            }
            this.utility = utility;
        }
    }

    /**
     * [新增方法, 借鉴自TALENT]
     * 更新Nettree，从叶子节点反向传播minLeave和maxLeave信息，并标记无法到达最终层级的节点。
     * @param nettree 要更新的Nettree
     * @param patternSize 模式的长度 (Nettree的总层数)
     */
    private void updateNetTreeWithLeaves(ArrayList<Node>[] nettree, int patternSize) {
        if (patternSize == 0) return;

        // 1. 初始化所有节点的toleave状态
        // 只有最后一层的节点能直接到达"终点"
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

        // 2. 从倒数第二层开始，向上遍历每一层
        for (int i = patternSize - 2; i >= 0; i--) {
            // 遍历当前层的每一个节点
            for (Node node : nettree[i]) {
                // 如果一个节点没有子节点，它肯定无法到达终点
                if (node.children.isEmpty()) {
                    continue;
                }

                int minL = Integer.MAX_VALUE;
                int maxL = Integer.MIN_VALUE;
                boolean canReachLeaf = false;

                // 遍历所有子节点，更新min/max leave信息
                for (Node child : node.children) {
                    // 只考虑那些自身能到达终点的子节点
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

                // 如果至少有一个子节点能到达终点，则当前节点也能
                if (canReachLeaf) {
                    node.toleave = true;
                    node.minLeave = minL;
                    node.maxLeave = maxL;
                }
            }
        }
    }

    /**
     * 检查并根据需要调整 reusableNettree 的大小。
     * @param requiredSize 当前模式所需的 Nettree 层数 (pattern.size())
     */
    private void checkAndResizeNettree(int requiredSize) {
        // 如果 Nettree 未初始化，或者当前分配的大小不足以容纳所需层数
        if (reusableNettree == null || currentNettreeSize < requiredSize) {

            // 释放旧的（如果存在）
            reusableNettree = null; // 帮助 GC

            // 分配新的数组
            reusableNettree = (ArrayList<Node>[]) new ArrayList[requiredSize];

            // 关键：必须初始化内部的 ArrayList
            for (int i = 0; i < requiredSize; i++) {
                reusableNettree[i] = new ArrayList<>();
            }
            currentNettreeSize = requiredSize;
        }
        // 如果 currentNettreeSize >= requiredSize，我们什么也不做，
        // 直接重用现有的、足够大的结构。
    }

    /**
     * 清理 Nettree 结构以供重用。
     * 只清理到指定的层数。
     * @param layersToClean 要清理的层数 (通常是 pattern.size())
     */
    private void cleanNettree(int layersToClean) {
        // 只清理我们需要用到的层数
        int layers = Math.min(layersToClean, this.currentNettreeSize);
        for (int i = 0; i < layers; i++) {
            // 清空列表，移除所有 Node 对象引用，但保留 ArrayList 对象本身
            reusableNettree[i].clear();
        }
    }

}
