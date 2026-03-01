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
	int[][] matrixItemTighterRemainingUtility;

	/** the item names (the list of all items in a sequence) 序列中包含的所有item的数组*/
	int[] itemNames;
	/** the swu of this sequence **/
	int swu;


	private ArrayList<ArrayList<Integer>> sequenceData;
	private Map<Integer, Integer> itemToIndexMap;
	
	/**
	 * Constructor
	 * @param nbItem the number of item in the sequence
	 * @param nbItemset the number of itemsets in that sequence
	 * @param sequence
	 */
	public QMatrix(int nbItem, int nbItemset, int[] itemNames, int itemNamesLength, int swu,
				   ArrayList<ArrayList<Integer>> sequence) {
		matrixItemUtility = new int[nbItem][nbItemset];
		matrixItemRemainingUtility = new int[nbItem][nbItemset];
		this.swu = swu;
		
		this.itemNames = new int[itemNamesLength];
		System.arraycopy(itemNames, 0, this.itemNames, 0, itemNamesLength);

		this.sequenceData = sequence;
		this.itemToIndexMap = new HashMap<>();
		for (int i = 0; i < itemNamesLength; i++) {
			this.itemToIndexMap.put(itemNames[i], i);
		}
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
	 * Get the number of itemsets in this sequence (Q-Matrix columns).
	 * @return the number of itemsets
	 */
	public int getNbItemsets(){
		return matrixItemUtility[0].length;
	}

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



	 public ArrayList<Integer> getPresentItemsInItemset(int itemsetIdx){
		 ArrayList<Integer> itemset = new ArrayList<>();
		 for(int i = 0; i < this.itemNames.length; i++) {
			 if(matrixItemUtility[i][itemsetIdx] > 0){
				 itemset.add(itemNames[i]);
			 }
		 }
		 return itemset;
	 }


	public ArrayList<ArrayList<Integer>> getSequenceAsItemsets(){
		return this.sequenceData;
	}


	public List<Integer> getItemsAfter(int startColumn) {
		List<Integer> items = new ArrayList<>();
		java.util.Set<Integer> seenItems = new java.util.HashSet<>();

		for (int j = startColumn; j < getNbItemsets(); j++) {
			for (int item : getPresentItemsInItemset(j)) {
				if (seenItems.add(item)) {
					items.add(item);
				}
			}
		}
		return items;
	}


	public int getFirstOccurrence(int item, int startColumn) {
		Integer rowIndex = itemToIndexMap.get(item);
		if (rowIndex == null) {
			return -1;
		}


		for (int j = startColumn; j < getNbItemsets(); j++) {
			if (matrixItemUtility[rowIndex][j] > 0) {
				return j;
			}
		}


		return -1;
	}

	public int[] getItemNames() {
		return this.itemNames;
	}

	public Map<Integer, Integer> getItemToIndexMap() {
		return this.itemToIndexMap;
	}


	public void computeTighterRemainingUtility(Set<Integer> unpromisingItems) {
		int nbRows = this.itemNames.length;
		int nbItemsets = getNbItemsets();
        this.matrixItemTighterRemainingUtility = new int[nbRows][nbItemsets];
        int tighterRemainingUtility = 0;


        for (int j = nbItemsets - 1; j >= 0; j--) {
            int currentItemsetUtility = 0;
            for (int i = 0; i < nbRows; i++) {
                int item = this.itemNames[i];
                int utility = this.matrixItemUtility[i][j];

                if (utility > 0 && !unpromisingItems.contains(item)) {
                    currentItemsetUtility += utility;
                }
            }

            for (int i = 0; i < nbRows; i++) {
                this.matrixItemTighterRemainingUtility[i][j] = tighterRemainingUtility;
            }

            tighterRemainingUtility += currentItemsetUtility;
        }
    }

}
