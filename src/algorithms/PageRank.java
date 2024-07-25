package algorithms;

import java.util.ArrayList;
import java.util.Arrays;

import graphs.Graph;

//TODO test me

public class PageRank {
	private Graph g;
	private ArrayList<Integer>[] neighbors;
	private double[] result;
	private double[] temp;
	private double d = 0.85;
	private int max_rounds = 20;
	private double convergence_threshold = 0.01d;
	
	public PageRank(Graph g){
		this.g = g;
		this.neighbors = g.get_neighbors();
		this.result = new double[g.num_vertices];
		this.temp = new double[g.num_vertices];
	}
	
	/**
	 * 
	 * @param g
	 * @return array of page rank value per node. Sum over all page rank value = 1.
	 */
	public static double[] run(Graph g) {
		return new PageRank(g).run();
	}
	
	double[] run() {
		//init: equally distribute PR
		Arrays.fill(result, 1.0d/(double)g.num_vertices);
		int round = 0;
		while(round<max_rounds) {
			System.out.print("round "+round+" ");
			double sum_delta = round();
			if(sum_delta<convergence_threshold ) {
				return result;
			}
			round++;
		}
		return result;
	}
	
	double round() {
		Arrays.fill(temp, 0);
		//first compute all the incoming page rank (from last round) 
		for(int node=0;node<g.num_vertices;node++) {
			ArrayList<Integer> my_neighbors = neighbors[node];
			double size = my_neighbors.size();
			double page_rank = result[node];//my page rank of the round before
			double give_away_pr = page_rank / size;
			for(int target_node : my_neighbors){
				temp[target_node] += give_away_pr;
			}
		}
		double sum_delta = sum_delta(temp, result);
		//normalize to 1
		double sum_pr = 0.0d;
		for(int node=0;node<g.num_vertices;node++) {
			result[node] = d*temp[node];
			result[node] += (1.0d-d)/(double)g.num_vertices;
			sum_pr+=result[node];
		}
		for(int node=0;node<g.num_vertices;node++) {
			result[node] /= sum_pr;
		}
		System.out.println("sum_delta="+sum_delta+" sum_pr="+sum_pr);
		return sum_delta;
	}

	double sum_delta(double[] incoming, double[] old_values) {
		double sum_delta = 0;
		for(int node=0;node<g.num_vertices;node++) {
			sum_delta += Math.abs(old_values[node]-incoming[node]);
		}
		return sum_delta;
	}
}
