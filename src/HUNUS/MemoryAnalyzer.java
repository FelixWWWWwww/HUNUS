package HUNUS;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

/**
 * 内存分析工具类，用于分析算法中各个数据结构的内存占用
 */
public class MemoryAnalyzer {
    
    /**
     * 分析 QMatrix 数据库的内存占用
     */
    public static MemoryInfo analyzeDatabase(List<QMatrix> database) {
        long totalBytes = 0;
        int sequenceCount = database.size();
        long totalMatrixBytes = 0;
        long totalSequenceDataBytes = 0;
        long totalMapBytes = 0;
        long totalItemNamesBytes = 0;
        
        for (QMatrix qm : database) {
            // 矩阵内存：matrixItemUtility + matrixItemRemainingUtility
            int nbItem = qm.getItemNames().length;
            int nbItemset = qm.getNbItemsets();
            // 二维数组开销：数组对象 + 行数组 + 数据
            long matrixOverhead = 24 + nbItem * 8; // 外层数组对象 + 行引用数组
            long matrixData = (long) nbItem * nbItemset * 4; // 实际数据（int 4字节）
            long matrixBytes = (matrixOverhead + matrixData) * 2; // 两个矩阵
            totalMatrixBytes += matrixBytes;
            
            // itemNames 数组
            long itemNamesBytes = 24 + (long) nbItem * 4; // 数组对象 + 数据
            totalItemNamesBytes += itemNamesBytes;
            
            // sequenceData (估算)
            long sequenceDataBytes = estimateSequenceDataSize(qm.getSequenceAsItemsets());
            totalSequenceDataBytes += sequenceDataBytes;
            
            // itemToIndexMap (估算：HashMap开销 + 每个entry)
            long mapOverhead = 48; // HashMap基础开销
            long mapEntries = (long) nbItem * 40; // 每个entry约40字节（key + value + 链表节点）
            long mapBytes = mapOverhead + mapEntries;
            totalMapBytes += mapBytes;
            
            // QMatrix对象本身的开销（约32字节）
            totalBytes += 32;
        }
        
        totalBytes = totalMatrixBytes + totalSequenceDataBytes + totalMapBytes + totalItemNamesBytes;
        
        return new MemoryInfo("QMatrix数据库", totalBytes, sequenceCount, 
            String.format("矩阵: %.2f MB, 序列数据: %.2f MB, 映射表: %.2f MB, 项名: %.2f MB [总字节: %d]",
                totalMatrixBytes / 1024.0 / 1024.0,
                totalSequenceDataBytes / 1024.0 / 1024.0,
                totalMapBytes / 1024.0 / 1024.0,
                totalItemNamesBytes / 1024.0 / 1024.0,
                totalBytes));
    }
    
    /**
     * 估算序列数据的内存占用
     */
    private static long estimateSequenceDataSize(ArrayList<ArrayList<Integer>> sequence) {
        long bytes = 0;
        // ArrayList 对象开销：约24字节 + 引用数组
        bytes += 24;
        bytes += sequence.size() * 8; // 引用数组（每个引用8字节）
        // 每个项集的 ArrayList
        for (ArrayList<Integer> itemset : sequence) {
            bytes += 24; // ArrayList对象开销
            bytes += itemset.size() * 8; // 引用数组（每个引用8字节）
            // Integer对象：在Java中，小整数可能被缓存，但这里保守估算
            bytes += itemset.size() * 16; // Integer对象（每个约16字节）
        }
        return bytes;
    }
    
    /**
     * 分析 LI_T 表的内存占用
     */
    public static MemoryInfo analyzeLITable(ArrayList<ArrayList<Integer>> LI_T) {
        long totalBytes = 0;
        int sequenceCount = LI_T.size();
        
        // 外层 ArrayList
        totalBytes += 24; // ArrayList对象
        totalBytes += sequenceCount * 8; // 引用数组
        
        // 每个内层 ArrayList
        for (ArrayList<Integer> row : LI_T) {
            if (row != null && !row.isEmpty()) {
                totalBytes += 24; // ArrayList对象
                totalBytes += row.size() * 8; // 引用数组
                // 注意：如果存储的是Integer对象（装箱），每个约16字节
                // 如果存储的是int原始值，每个4字节
                // 这里假设是Integer对象（更保守的估算）
                totalBytes += row.size() * 16; // Integer对象
            }
        }
        
        return new MemoryInfo("LI_T表", totalBytes, sequenceCount, 
            String.format("平均每序列: %.2f KB [总字节: %d]", (totalBytes / 1024.0) / Math.max(1, sequenceCount), totalBytes));
    }
    
    /**
     * 分析投影链的内存占用
     */
    public static MemoryInfo analyzeProjectionChains(Map<Integer, ArrayList<AlgoHUNUS.TargetedList>> mapItemUC) {
        long totalBytes = 0;
        int itemCount = mapItemUC.size();
        int totalTargetedLists = 0;
        int totalUtilityElements = 0;
        
        // HashMap 本身的开销
        totalBytes += 48; // HashMap对象开销
        totalBytes += mapItemUC.size() * 32; // 每个entry的开销
        
        for (Map.Entry<Integer, ArrayList<AlgoHUNUS.TargetedList>> entry : mapItemUC.entrySet()) {
            // Integer key
            totalBytes += 16;
            
            // ArrayList<TargetedList>
            ArrayList<AlgoHUNUS.TargetedList> chain = entry.getValue();
            totalBytes += 24; // ArrayList对象
            totalBytes += chain.size() * 8; // 引用数组
            
            totalTargetedLists += chain.size();
            
            // 每个 TargetedList
            for (AlgoHUNUS.TargetedList tl : chain) {
                totalBytes += 32; // TargetedList对象开销
                totalBytes += 24; // List<UtilityElement> 对象
                totalBytes += tl.LengthOfUtilityList * 8; // 引用数组
                
                totalUtilityElements += tl.LengthOfUtilityList;
                
                // 每个 UtilityElement
                totalBytes += tl.LengthOfUtilityList * 24; // UtilityElement对象（3个int字段）
            }
        }
        
        return new MemoryInfo("投影链(mapItemUC)", totalBytes, itemCount,
            String.format("TargetedList数: %d, UtilityElement数: %d, 平均每项: %.2f KB [总字节: %d]",
                totalTargetedLists, totalUtilityElements, 
                (totalBytes / 1024.0) / Math.max(1, itemCount),
                totalBytes));
    }
    
    /**
     * 分析结果模式的内存占用
     */
    public static MemoryInfo analyzeFoundPatterns(List<AlgoHUNUS.ResultPattern> foundPatterns) {
        long totalBytes = 0;
        int patternCount = foundPatterns.size();
        
        // ArrayList 开销
        totalBytes += 24;
        totalBytes += foundPatterns.size() * 8; // 引用数组
        
        // 每个 ResultPattern
        for (AlgoHUNUS.ResultPattern rp : foundPatterns) {
            totalBytes += 32; // ResultPattern对象
            totalBytes += 24; // pattern ArrayList
            totalBytes += rp.pattern.size() * 8; // 项集引用数组
            
            // 每个项集
            for (ArrayList<Integer> itemset : rp.pattern) {
                totalBytes += 24; // ArrayList对象
                totalBytes += itemset.size() * 4; // int值
                totalBytes += itemset.size() * 16; // Integer对象（估算）
            }
        }
        
        return new MemoryInfo("结果模式(foundPatterns)", totalBytes, patternCount,
            String.format("平均每个模式: %.2f KB [总字节: %d]", (totalBytes / 1024.0) / Math.max(1, patternCount), totalBytes));
    }
    
    /**
     * 内存信息类
     */
    public static class MemoryInfo {
        public final String name;
        public final long bytes;
        public final int count;
        public final String details;
        
        public MemoryInfo(String name, long bytes, int count, String details) {
            this.name = name;
            this.bytes = bytes;
            this.count = count;
            this.details = details;
        }
        
        public double getMB() {
            return bytes / 1024.0 / 1024.0;
        }
        
        @Override
        public String toString() {
            double mb = getMB();
            double kb = bytes / 1024.0;
            if (mb < 0.01) {
                // 如果小于 0.01 MB，显示为 KB
                return String.format("%s: %.2f KB (数量: %d, %s)", name, kb, count, details);
            } else {
                return String.format("%s: %.2f MB (%.2f KB) (数量: %d, %s)", name, mb, kb, count, details);
            }
        }
    }
}

