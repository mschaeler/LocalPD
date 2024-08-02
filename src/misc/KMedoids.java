package misc;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

public class KMedoids {
	public static final int MAX_ITERATIONS = 10;//TODO meaningfull value
	
	final int[][] data;
	final int[][] dist_matrix;
	final int k;
	final int num_points;
	
	Clustering myClusters;
	
	public KMedoids(final int[][] data, int k) {
		num_points = data.length;
		this.dist_matrix = new int[num_points][num_points];
		this.data = data;
		this.k = k;
		fill_dist_matrix();
	}

	private void fill_dist_matrix() {
		for(int i=0;i<num_points;i++) {
			for(int j=+1;j<num_points;j++) {//Exploit commutivity
				int dist = l_1_norm(data[i], data[j]);
				dist_matrix[i][j] = dist;
				dist_matrix[j][i] = dist;
			}
		}
	}

	private int l_1_norm(int[] p_1, int[] p_2) {
		int dist = 0;
		for(int i=0;i<p_1.length;i++) {
			dist+=Math.abs(p_1[i]-p_2[i]);
		}
		return dist;
	}

	/**
	 * initial medoids, i.e., their TIDs
	 * @return
	 */
	int[] getMedoids(){
		myClusters = new Clustering(k);
		//initial medoids, i.e., their TIDs
		int[] medoids = myClusters.medoids;

		Random rand = new Random(1234567);
		for(int i=0;i<medoids.length;i++){//XXX wird schon keiner doppelt ausgewählt
			int tid = rand.nextInt(num_points);
			medoids[i] = tid;
		}
		return medoids;
	}
	
	/**
	 * initial medoids, i.e., their TIDs
	 * @return
	 */
	int[] getMedoids2(){
		myClusters = new Clustering(k);
		//initial medoids, i.e., their TIDs
		int[] medoids = myClusters.medoids;

		for(int i=0;i<medoids.length;i++){
			medoids[i] = i;//Very stupid :)
		}
		return medoids;
	}
	
	public Clustering run(){
		getMedoids();
		assignAllPointsToNerarestMedoid();
		
		int iteration = 0;
		boolean changed = true;
		while(changed && (iteration++ < MAX_ITERATIONS)){
			changed = false;
			
			// Try to swap the medoid with a better cluster member
			for(int cluster = 0; cluster<myClusters.NUM_CLUSTERS;cluster++){
				boolean temp = tryFindBetterMedoid(myClusters.getCluster(cluster),cluster); 
				if(!changed){//only consider temp if not at least one medoid alreday changed
					changed = temp;
				}
			}
			//re-assign points
			assignAllPointsToNerarestMedoid();
		}
		return myClusters;
	}

	private boolean tryFindBetterMedoid(ArrayList<Integer> TIDs, int cluster) {
		if(TIDs.isEmpty()) {
			return false;
		}
		double[] cost = new double[TIDs.size()];
		
		//determine all cost, i.e., all pair-wise distances in cluster
		for(int i=0;i<TIDs.size();i++){
			int tid = TIDs.get(i);
			cost[i] = getCost(tid,TIDs);
		}
		// find the point with mininmal cost
		double minCost = cost[0];
		int minTID = TIDs.get(0);
		for(int i=1;i<cost.length;i++){
			if(cost[i]<minCost){
				minCost = cost[i];
				minTID = TIDs.get(i);
			}
		}
		//determine whether we found new medoid
		if(myClusters.medoids[cluster]!=minTID){
			myClusters.medoids[cluster] = minTID;
			return true;
		}else{//no, medoid remins the same
			return false;
		}
	}

	private double getCost(final int tid, final ArrayList<Integer> TIDs) {
		double sumOfCost = 0.0d;
		for(int i=0;i<TIDs.size();i++){
			int toCompare = TIDs.get(i);
			sumOfCost += dist(tid, toCompare);
			//sumOfCost += DIST_MATRIX[tid][toCompare];
		}
		return sumOfCost;
	}

	private void assignAllPointsToNerarestMedoid() {
		myClusters.clear();
		for(int tid=0;tid<num_points;tid++){
			int cluster = 0;
			double minDist = dist(tid,myClusters.medoids[0]);
			
			for(int i=1;i<myClusters.NUM_CLUSTERS;i++){
				double dist = dist(tid,myClusters.medoids[i]);
				if(dist<minDist){
					minDist = dist;
					cluster = i;
				}
			}
			myClusters.ALL_CLUSTERS[cluster].add(tid);
		}
		
	};
	
	final int dist(int tid_1, int tid_2) {
		return this.dist_matrix[tid_1][tid_2];
	}
	
	public class Clustering extends AbstractClustering{
		final int[] medoids = new int[k];
		
		@SuppressWarnings("unchecked")
		Clustering(int numCluster){
			super(numCluster);
		}
		
		public String toString(){
			String ret = "";
			ret += "Clustering result\n";
			for(int cluster=0;cluster<NUM_CLUSTERS;cluster++){
				//Collections.sort(getCluster(cluster));
				ret+="Cluster "+cluster+" of size "+getCluster(cluster).size()+" Medoid "+medoids[cluster]+" "+Arrays.toString(getCluster(cluster).toArray())+"\n";
			}
			return ret;
		}
		
		double clusterDissimilarityForPoint(int tid, ArrayList<Integer> objectsInThisCluster, boolean ownCluster) {
			double totalDissimilarity = 0.0;
			for (int anotherTID : objectsInThisCluster) {
				totalDissimilarity += dist(tid,anotherTID);
			}
			if(ownCluster)
				return totalDissimilarity / (objectsInThisCluster.size() - 1);// I do not count
			else 
				return totalDissimilarity / (objectsInThisCluster.size());
		}
		double avgDistToNaerestCluster(int tid, int cluster) {
			double avgDist = Double.MAX_VALUE;//looking for smallest average distance
			for (int i=0;i<NUM_CLUSTERS;i++) {
				if(i == cluster){continue;}//not me
				avgDist = Math.min(avgDist, clusterDissimilarityForPoint(tid, getCluster(i), false));
			}
			return avgDist;
		}
	}
}
