package graphs;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Random;

import misc.Util;
import results.Config;

public class RandomSamplePartition {
	ArrayList<Integer> my_edges;
	public RandomSamplePartition(int size){
		my_edges = new ArrayList<Integer>(size);
	}
	public void add_edge(int new_edge) {
		my_edges.add(new_edge);
	}
	public int sample(HashSet<Integer> my_neighbors, double epsilon_q2, Random rand) {
		final double my_scores[] = new double[my_edges.size()];
		
		//(1) Compute the scores
		for(int i=0;i<my_edges.size();i++) {
			int edge = my_edges.get(i);
			if(my_neighbors.contains(edge)) {
				my_scores[i] = Sample.SCORE_EXISTING_EDGE;//for all other edges it remains zero
			}//else : is already zero
		}
		//(2) Compute u()
		for(int i=0;i<my_scores.length;i++) {
			double temp = epsilon_q2*my_scores[i]/(2.0d*Sample.delta_u());
			my_scores[i] = Math.pow(Math.E,temp);
		}
		int edge_to_return = sample(my_scores, rand);
		return edge_to_return;
	}
	
	
	public int sample(final double my_scores[], Random rand) {
		final double sum = Util.sum(my_scores);
		final double diced_score = rand.nextDouble()*sum;
		double sum_score = 0.0d;
		for(int position=0;position<my_scores.length;position++) {
			sum_score+=my_scores[position];
			if(sum_score>=diced_score) {
				int edge = my_edges.get(position);
				return edge;
			}
		}
		System.err.println("sample(double) should never come to here");
		return -1;
	}
	
	public int noisy_max(HashSet<Integer> my_neighbors, double epsilon_q2, Random rand) {
		final double my_scores[] = new double[my_edges.size()];
		
		//(1) Compute the scores
		for(int i=0;i<my_edges.size();i++) {
			int edge = my_edges.get(i);
			if(my_neighbors.contains(edge)) {
				my_scores[i] = Sample.SCORE_EXISTING_EDGE;//for all other edges it remains zero
			}//else : is already zero
		}
		//(2) Sanitize the scores and get noisy max
		double my_max = Double.NEGATIVE_INFINITY;
		int my_max_index = 0;
		double sensitivity = Sample.SCORE_EXISTING_EDGE-Sample.SCORE_NONE_EDGE;
		for(int i=0;i<my_scores.length;i++) {
			double score = my_scores[i];
			double noisy_score;
			if(Config.NON_PRIVATE_NOISY_MAX) {
				noisy_score = score;	
			}else{
				noisy_score = Mechanism.sanitize(score, sensitivity, epsilon_q2);
			}
			 
			if(score>my_max){
				my_max = noisy_score;
				my_max_index = i;
			}
		}
		
		int edge = my_edges.get(my_max_index);
		return edge;
	}
	
}
