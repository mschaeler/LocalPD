package misc;

import cern.colt.Arrays;

public final class TopK {
	public final int[] node_ids;
	final double[] node_scores;
	double min_score_in_top_k = Double.POSITIVE_INFINITY;
	int index_min_score = -1;
	final int k;
	
	public TopK(final int k){
		this.k=k;
		this.node_ids = new int[k];
		this.node_scores = new double[k];
	}
	
	/**
	 * update only if distance smaller
	 * @param TID
	 * @return
	 */
	public final boolean tryUpdate(final int node_id, final double node_score){
		if(node_score>this.min_score_in_top_k){
			update(node_score, node_id);		
			return true;
		}else{
			return false;
		}
	}

	/**
	 * No check - safe variant: tryUpdate(final int TID)  
	 * @param distance
	 * @param TID
	 * @return
	 */
	private final void update(final double score, final int id){
		node_scores[index_min_score]	= score;
		node_ids[index_min_score]		= id;
		
		//determine new *maxDistance*
		min_score_in_top_k=Double.POSITIVE_INFINITY;
		for(int j=0;j<k; j++){
			if (node_scores[j]<min_score_in_top_k){
				min_score_in_top_k=node_scores[j];
				index_min_score=j;
			}
		}
	}
	
	public static TopK create_With_First_K_Elements(final int k, final double[] scores){
		TopK col = new TopK(k);		
		col.index_min_score = 0;
		
		//insert data & scores
		for(int tid=0;tid<k;tid++){
			col.node_ids[tid] 	= tid;
			final double score 	= scores[tid];
			col.node_scores[tid]= score;
			
			if(score<col.min_score_in_top_k){
				col.min_score_in_top_k 	= score;
				col.index_min_score 	= tid;
			}
		}
		
		return col;
	}
	
	public String toString() {
		String ret = "min()="+this.min_score_in_top_k+" min_index="+this.index_min_score+"\n";
		ret += Arrays.toString(this.node_ids)+"\n"+Arrays.toString(this.node_scores);
		return ret;
	}
}
