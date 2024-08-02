package misc;

import java.util.ArrayList;
import java.util.HashMap;

public abstract class AbstractClustering {
	public ArrayList<Integer>[] ALL_CLUSTERS;
	final int NUM_CLUSTERS;
	
	public AbstractClustering(int numCluster) {
		NUM_CLUSTERS = numCluster;
		ALL_CLUSTERS = new ArrayList[numCluster];
		for(int i=0;i<NUM_CLUSTERS;i++){
			ALL_CLUSTERS[i] = new ArrayList<Integer>();
		}
	}
	
	void clear() {
		for(int i=0;i<NUM_CLUSTERS;i++){
			ALL_CLUSTERS[i].clear();
		}
	}
	
	public ArrayList<Integer> getCluster(int cluster) {
		return ALL_CLUSTERS[cluster];
	}
	
	/**
	 * https://en.wikipedia.org/wiki/Silhouette_(clustering)
	 * https://github.com/OryxProject/oryx/blob/master/app/oryx-app-mllib/src/main/java/com/cloudera/oryx/app/batch/mllib/kmeans/SilhouetteCoefficient.java
	 * @return
	 */
	public double silhouetteCoefficient() {
		double totalSilhouetteCoefficient = 0.0;
		long sampleCount = 0L;
	
	    for (int i=0;i<NUM_CLUSTERS;i++) {
	    	ArrayList<Integer> objectsInThisCluster = getCluster(i);
	    	int clusterSize = objectsInThisCluster.size();
	    	// Increment the total sample count for computing silhouette coefficient
	    	sampleCount += clusterSize;
	
	        for (int p=0;p<clusterSize;p++) {
	        	int tid = objectsInThisCluster.get(p);
	        	double a_i = clusterDissimilarityForPoint(tid, objectsInThisCluster, true);
	        	double b_i = avgDistToNaerestCluster(tid, i);
	        	totalSilhouetteCoefficient += s_i(a_i, b_i);
	        }
	    }
	    return totalSilhouetteCoefficient / sampleCount;
	}
	
	abstract double clusterDissimilarityForPoint(int tid, ArrayList<Integer> objectsInThisCluster, boolean b);
	abstract double avgDistToNaerestCluster(int tid, int cluster);


	double s_i(double ai, double bi) {
		if (ai < bi) {
			return 1.0 - (ai / bi);
		}else if(ai > bi) {
			return (bi / ai) - 1.0;
		}else{
			return 0.0;
		}
	}
	public HashMap<Integer, Integer> to_hash_map(){
		HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
		for(int cluster=0;cluster<NUM_CLUSTERS;cluster++) {
			for(int id : getCluster(cluster)) {
				map.put(id, cluster);
			}
		}
		return map;
	}
}
