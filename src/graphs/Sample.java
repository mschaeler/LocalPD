package graphs;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Random;

import misc.Util;

public class Sample {
	static double SCORE_EXISTING_EDGE 	= 1;
	static double SCORE_NON_EDGE 		= 0;
	
	final double my_scores[];
	HashSet<Integer> index;
	Random rand = new Random(Util.seed);
	
	public Sample(Graph g, ArrayList<Integer> my_neighbors){
		my_scores = new double[g.num_vertices];
		index = new HashSet<Integer>(my_neighbors.size());
		for(int e : my_neighbors) {
			index.add(e);
		}
	}
	
	double score(int edge){
		if(index.contains(edge)) {
			return SCORE_EXISTING_EDGE;
		}else{
			return SCORE_NON_EDGE;
		}
	}
	
	public int sample(final double epsilon) {
		for(int edge=0;edge<my_scores.length;edge++) {
			my_scores[edge] = Math.pow(Math.E,epsilon*score(edge)/(2*delta_u()));
		}
		final double sum = Util.sum(my_scores);
		final double diced_score = rand.nextDouble()*sum;
		double sum_score = 0.0d;
		for(int edge=0;edge<my_scores.length;edge++) {
			sum_score+=my_scores[edge];
			if(sum_score>=diced_score) {
				return edge;
			}
		}
		System.err.println("sample(double) should never come to here");
		return -1;
	}
	
	public static double delta_u() {
		return SCORE_EXISTING_EDGE-SCORE_NON_EDGE;
	}
}
