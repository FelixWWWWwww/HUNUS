package HUNUS;


/* Copyright (c) 2008-2015 Philippe Fournier-Viger
* 
* This file is part of the SPMF DATA MINING SOFTWARE
* (http://www.philippe-fournier-viger.com/spmf).
* 
* SPMF is free software: you can redistribute it and/or modify it under the
* terms of the GNU General Public License as published by the Free Software
* Foundation, either version 3 of the License, or (at your option) any later
* version.
* 
* SPMF is distributed in the hope that it will be useful, but WITHOUT ANY
* WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
* A PARTICULAR PURPOSE. See the GNU General Public License for more details.
* 
* You should have received a copy of the GNU General Public License along with
* SPMF. If not, see <http://www.gnu.org/licenses/>.
*/

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * This class represents a QMatrix as described in the paper describing USPan
 * @author Philippe Fournier-Viger, 2015
 * @see AlgoUSpan
 */
public class QMatrix {

	/** the qmatrix for items  [item][itemset] -> utility  */
	int[][] matrixItemUtility;
	/** the qmatrix for remaining utility [item][itemset] -> remaining utility*/
	int[][] matrixItemRemainingUtility;
	/** [新] 更紧凑的剩余效用 */
	int[][] matrixItemTighterRemainingUtility;

	/** the item names (the list of all items in a sequence) 序列中包含的所有item的数组*/
	int[] itemNames;
	/** the swu of this sequence **/
	int swu;

	// --- 新增的核心数据结构 ---
	// 1. 存储标准格式的序列数据，避免重复构建
	private ArrayList<ArrayList<Integer>> sequenceData;
	// 2. 存储item到其在矩阵中行号的映射，实现O(1)查询
	private Map<Integer, Integer> itemToIndexMap;
	// --- 新增结束 ---
	
	/**
	 * Constructor
	 * @param nbItem the number of item in the sequence
	 * @param nbItemset the number of itemsets in that sequence
	 * @param sequence 序列的标准格式数据
	 */
	public QMatrix(int nbItem, int nbItemset, int[] itemNames, int itemNamesLength, int swu,
				   ArrayList<ArrayList<Integer>> sequence) {
		matrixItemUtility = new int[nbItem][nbItemset];
		matrixItemRemainingUtility = new int[nbItem][nbItemset];
		this.swu = swu;
		
		this.itemNames = new int[itemNamesLength];
		System.arraycopy(itemNames, 0, this.itemNames, 0, itemNamesLength);

		// --- 初始化新增的数据结构 ---
		this.sequenceData = sequence; // 直接保存，不再重新计算
		this.itemToIndexMap = new HashMap<>();
		for (int i = 0; i < itemNamesLength; i++) {
			this.itemToIndexMap.put(itemNames[i], i);
		}
		// --- 初始化结束 ---
	}

	public int getItemIndex(int item) {
		Integer index = itemToIndexMap.get(item);
		return (index == null) ? -1 : index;
	}

	/**
	 * Register item in the matrix
	 * @param itemPos an item position in "itemNames"
	 * @param itemset the itemset number
	 * @param utility the utility of the item in that itemset
	 * @param remainingUtility the reamining utility of that item at that itemset
	 */
	public void registerItem(int itemPos, int itemset, int utility, int remainingUtility) {
		// we store the utility in the cell for this item/itemset
		matrixItemUtility[itemPos][itemset] = utility;
		// we store the remaining utility in the cell for this item/itemset
		matrixItemRemainingUtility[itemPos][itemset] = remainingUtility;
	}

	
	/**
	 * Get a string representation of this matrix (for debugging purposes)
	 * @return the string representation
	 */
	public String toString() {
		StringBuffer buffer = new StringBuffer();
		buffer.append("QMatrix信息:\n");
		buffer.append("序列效用(SWU): " + swu + "\n");
		buffer.append("项集数量: " + matrixItemUtility[0].length + "\n");
		buffer.append("项数量: " + itemNames.length + "\n");
		buffer.append("项列表: ");
		for(int i = 0; i < itemNames.length; i++) {
			buffer.append(itemNames[i]);
			if(i < itemNames.length - 1) buffer.append(", ");
		}
		buffer.append("\n\n");
		
		// 打印矩阵表头
		buffer.append("效用矩阵 (格式: 效用[剩余效用]):\n");
		buffer.append("项目\\项集\t");
		for(int j = 0; j < matrixItemUtility[0].length; j++) {
			buffer.append("项集" + j + "\t");
		}
		buffer.append("\n");
		
		// 打印矩阵内容
		for(int i=0; i< itemNames.length; i++) {
			buffer.append("项" + itemNames[i] + "\t");
			for(int j=0; j< matrixItemUtility[i].length; j++) {
				buffer.append(matrixItemUtility[i][j] + "[" + 
						matrixItemRemainingUtility[i][j] + "]\t");
			}
			buffer.append("\n");
		}
		buffer.append("\n");
		return buffer.toString();
	}

	/**
	 * Get the number of itemsets in this sequence (Q-Matrix columns).
	 * @return the number of itemsets
	 */
	public int getNbItemsets(){
		return matrixItemUtility[0].length;
	}

	/**
	 * 获取指定项集（itemset）中的所有 item。
	 * @param itemsetIdx 项集索引（列号）
	 * @return 包含该项集所有有效 item 的列表
	 */
	public List<Integer> getItemset(int itemsetIdx) {
		List<Integer> itemset = new ArrayList<>();
		for (int i = 0; i < itemNames.length; i++) {
			itemset.add(itemNames[i]);
		}
		return itemset;
	}

	public int getUtility(int itemsetIdx, int item) {
		Integer index = itemToIndexMap.get(item);
		if (index == null) {
			return 0;
		}
		return this.matrixItemUtility[index][itemsetIdx];
	}

	public int getRemainingUtility(int itemset, int item) {
		Integer index = itemToIndexMap.get(item);
		if (index == null) {
			return 0;
		}
		return matrixItemRemainingUtility[index][itemset];
	}

	/**
	 * [新方法] 获取一个项在给定项集的 "紧凑" 剩余效用。
	 * @param itemset the itemset
	 * @param item the item
	 * @return the tighter remaining utility
	 */
	public int getTighterRemainingUtility(int itemset, int item) {
		// 如果尚未计算，则返回普通剩余效用作为后备
		if (matrixItemTighterRemainingUtility == null) {
			return getRemainingUtility(itemset, item);
		}
		Integer index = itemToIndexMap.get(item);
		if (index == null) {
			return 0;
		}
		return matrixItemTighterRemainingUtility[index][itemset];
	}

	/**
	 * Get the column count
	 * @return the column count
	 */
	public int getNbColumns() {
		return matrixItemUtility[0].length;
	}

	/**
	 * 高效的方法：
	 * 只获取指定项集（itemset）中实际存在（效用>0）的 item。
	 * @param itemsetIdx 项集索引（列号）
	 * @return 包含该项集所有有效 item 的列表
	 */
	 public ArrayList<Integer> getPresentItemsInItemset(int itemsetIdx){
		 ArrayList<Integer> itemset = new ArrayList<>();
		 for(int i = 0; i < this.itemNames.length; i++) {
			 if(matrixItemUtility[i][itemsetIdx] > 0){
				 itemset.add(itemNames[i]);
			 }
		 }
		 return itemset;
	 }

	/**
	 * 将QMatrix内部数据转换为项集序列，每个项集只包含效用>0的item
	 * @return 包含所有项集的列表，每个项集是一个包含有效item的列表
	 * (已优化) 直接返回预存的序列数据
	 */
	public ArrayList<ArrayList<Integer>> getSequenceAsItemsets(){
		return this.sequenceData;
	}

	/**
	 * [新增] 获取从指定项集索引开始的所有后续项。
	 * 用于查找 s-extensions。
	 * @param startColumn 开始搜索的项集索引 (包含此项集)
	 * @return 一个包含所有不重复项的列表
	 */
	public List<Integer> getItemsAfter(int startColumn) {
		List<Integer> items = new ArrayList<>();
		// 使用一个HashSet来保证项的唯一性，避免重复添加
		java.util.Set<Integer> seenItems = new java.util.HashSet<>();

		// 从指定的 startColumn 开始遍历到最后一个项集
		for (int j = startColumn; j < getNbItemsets(); j++) {
			// 遍历该项集中的所有实际存在的项
			for (int item : getPresentItemsInItemset(j)) {
				// 如果这个项还没有被添加过
				if (seenItems.add(item)) {
					items.add(item);
				}
			}
		}
		return items;
	}

	/**
	 * [新增] 查找一个项在指定项集索引之后（或之中）的首次出现位置。
	 * 用于 s-extension 的位置检查。
	 * @param item 要查找的项
	 * @param startColumn 开始搜索的项集索引 (包含此项集)
	 * @return 如果找到，返回项集索引；否则返回 -1。
	 */
	public int getFirstOccurrence(int item, int startColumn) {
		Integer rowIndex = itemToIndexMap.get(item);
		// 如果该项根本不存在于此序列的QMatrix中，直接返回-1
		if (rowIndex == null) {
			return -1;
		}

		// 从指定的 startColumn 开始遍历到最后一个项集
		for (int j = startColumn; j < getNbItemsets(); j++) {
			// 检查该项在当前项集的效用是否大于0
			if (matrixItemUtility[rowIndex][j] > 0) {
				// 如果是，说明项存在于此项集，返回当前项集索引
				return j;
			}
		}

		// 如果遍历完所有后续项集都未找到，返回-1
		return -1;
	}

	public int[] getItemNames() {
		return this.itemNames;
	}

	public Map<Integer, Integer> getItemToIndexMap() {
		return this.itemToIndexMap;
	}

	/**
	 * [新方法] 预计算更紧凑的剩余效用。
	 * 这个方法只累加那些 "有希望的" 项的效用。
	 * @param unpromisingItems 一个包含所有 "无望项" 的集合。
	 */
	public void computeTighterRemainingUtility(Set<Integer> unpromisingItems) {
		int nbRows = this.itemNames.length;
		int nbItemsets = getNbItemsets();
        this.matrixItemTighterRemainingUtility = new int[nbRows][nbItemsets];
        int tighterRemainingUtility = 0;

        // 从序列的最后一个项集向前遍历
        for (int j = nbItemsets - 1; j >= 0; j--) {
            int currentItemsetUtility = 0;
            // 遍历当前项集中的所有项
            for (int i = 0; i < nbRows; i++) {
                int item = this.itemNames[i];
                int utility = this.matrixItemUtility[i][j];

                // 只有当这个项是 "有希望的" (即不在unpromising集合中) 且效用不为0时，才累加
                if (utility > 0 && !unpromisingItems.contains(item)) {
                    currentItemsetUtility += utility;
                }
            }

            // 对于当前项集 j 中的所有项，它们的紧凑剩余效用是相同的，
            // 即在 j 之后所有 "有希望的" 项的效用之和。
            for (int i = 0; i < nbRows; i++) {
                this.matrixItemTighterRemainingUtility[i][j] = tighterRemainingUtility;
            }

            // 更新下一次迭代的剩余效用值
            tighterRemainingUtility += currentItemsetUtility;
        }
    }

}
