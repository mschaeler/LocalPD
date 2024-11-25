package graphs;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.HashMap;
import java.util.Map.Entry;

import algorithms.PageRank;
import algorithms.ProximityPrestige;

public class Metrics {
	public static final int correct_edge 		= 0;
	public static final int missing_edge 		= 1;
	public static final int new_fake_edge 		= 2;
	public static final int correct_none_edge 	= 3;
	
	static void out(String algo_name, double[] all_eps, ArrayList<double[]> all_metrics ) {
		System.out.println(algo_name);
		System.out.println("budget\tcorrect\tmissing\tfake edge\tnon edge");
		for(int i=0;i<all_eps.length;i++) {
			System.out.print("eps="+all_eps[i]+"\t");
			System.out.println(to_tsv(all_metrics.get(i)));
		}
	}
	
	private static String to_tsv(double[] arr) {
		if(arr==null) {
			return "null";
		}
		String s = "";
		for(double d : arr) {
			s+=d+"\t";
		};
		return s;
	}

	public static final int[] edge_edit_dist(ArrayList<Integer> l_1, ArrayList<Integer> l_2, int num_vertices) {
		boolean[] l_1_n = new boolean[num_vertices];
		boolean[] l_2_n = new boolean[num_vertices];
		for(int e : l_1) {
			l_1_n[e] = true;
		}
		for(int e : l_2) {
			l_2_n[e] = true;
		}
		int[] error = new int[4];
		for(int e=0;e<num_vertices;e++) {
			if(l_1_n[e] && l_2_n[e]) {
				error[correct_edge]++;
			}else if(l_1_n[e] && !l_2_n[e]){//there should be an edged, but there is none
				error[missing_edge]++;
			}else if(!l_1_n[e] && l_2_n[e]){//we invented a new edge
				error[new_fake_edge]++;
			}else{//in both cases there is no edge
				error[correct_none_edge]++;
			}
		}
		return error;
	}
	
	public static final double[] avg_edge_edit_dist(final BitSet[] ground_truth, final int num_vertices, ArrayList<Graph> sanitized_gs) {
		 
		ArrayList<double[]> all_sums = new ArrayList<double[]>();
		
		for(int i=0;i<sanitized_gs.size();i++) {
			double[] sum = new double[4];
			Graph g = sanitized_gs.get(i);
			final BitSet[] to_ceck = g.get_adjancency_matrix_as_bit_vector();
			final int[][] result = edge_edit_dist(ground_truth, to_ceck);
			
			for(int[] r : result) {
				sum[correct_edge]+= r[correct_edge];
				sum[missing_edge]+= r[missing_edge];
				sum[new_fake_edge]+= r[new_fake_edge];
				sum[correct_none_edge]+= r[correct_none_edge];
			}
			
			sum[correct_edge]/= 	num_vertices;
			sum[missing_edge]/= 	num_vertices;
			sum[new_fake_edge]/= 	num_vertices;
			sum[correct_none_edge]/=num_vertices;
			all_sums.add(sum);
		}

		double[] avg = new double[4];
		for(double[] arr : all_sums) {
			avg[0] += arr[0];
			avg[1] += arr[1];
			avg[2] += arr[2];
			avg[3] += arr[3];
		}
		avg[0] /= all_sums.size();
		avg[1] /= all_sums.size();
		avg[2] /= all_sums.size();
		avg[3] /= all_sums.size();
		
		return avg;
	}
	
	static final double mae(Graph original_g, Graph sanitized_g) {
		final BitSet[] ground_truth = original_g.get_adjancency_matrix_as_bit_vector();
		final BitSet[] to_ceck = sanitized_g.get_adjancency_matrix_as_bit_vector();
		final int[][] result = edge_edit_dist(ground_truth, to_ceck);
		
		double sum = 0;
		for(int[] r : result) {
			sum+= r[missing_edge];
			sum+= r[new_fake_edge];
		}
		
		double mae = sum / original_g.num_vertices; //normalize by Graph size (i.e., |V|)
		
		return mae;
	}
	
	static final double[] mae(Graph original_g, ArrayList<Graph> sanitized_gs) {
		final BitSet[] ground_truth = original_g.get_adjancency_matrix_as_bit_vector();
		double[] all_mae = new double[original_g.num_vertices];
		
		for(int i=0;i<sanitized_gs.size();i++) {
			Graph g = sanitized_gs.get(i);
			final BitSet[] to_ceck = g.get_adjancency_matrix_as_bit_vector();
			final int[][] result = edge_edit_dist(ground_truth, to_ceck);
			
			double sum = 0;
			for(int[] r : result) {
				sum+= r[missing_edge];
				sum+= r[new_fake_edge];
			}
			
			double mae = sum / original_g.num_vertices; //normalize by Graph size (i.e., |V|)
			all_mae[i] = mae;
		}

		return all_mae;
	}
	
	/**
	 * Computes for each vertex the number of correct, missing, fake, and correct-non-edges. I.e., the it return an array[|V|][4]
	 * @param ground_truth
	 * @param to_ceck
	 * @return
	 */
	private static int[][] edge_edit_dist(final BitSet[] ground_truth, final BitSet[] to_ceck) {
		if(ground_truth.length!=to_ceck.length) {
			System.err.println("node_edit_dist() ground_truth.length!=to_ceck.length");
		}
		final int[][] result = new int[ground_truth.length][4];
		for(int node = 0;node<ground_truth.length;node++) {
			for(int e = 0;e<ground_truth.length;e++) {
				if(ground_truth[node].get(e) && to_ceck[node].get(e)) {
					result[node][correct_edge]++;
				}else if(ground_truth[node].get(e) && !to_ceck[node].get(e)){//there should be an edged, but there is none
					result[node][missing_edge]++;
				}else if(!ground_truth[node].get(e) && to_ceck[node].get(e)){//we invented a new edge
					result[node][new_fake_edge]++;
				}else{//in both cases there is no edge
					result[node][correct_none_edge]++;
				}
			}
		}
		return result;
	}

	static final double mre(Graph original_g, Graph sanitized_g) {		
		final BitSet[] ground_truth = original_g.get_adjancency_matrix_as_bit_vector();
		ArrayList<Integer>[] neighbors = original_g.get_neighbors();
		final BitSet[] to_ceck = sanitized_g.get_adjancency_matrix_as_bit_vector();
		final int[][] result = edge_edit_dist(ground_truth, to_ceck);
		
		double sum = 0;
		for(int node=0;node<original_g.num_vertices;node++) {
			int[] r = result[node];
			double mae_node = r[missing_edge] + r[new_fake_edge];
			double mre_node = mae_node / (double)neighbors[node].size();//normalize by original node count
			sum+= mre_node;
		}
		
		double mre = sum / original_g.num_vertices; //normalize by Graph size (i.e., |V|)
		
		return mre;
	}
	
	static final double[] mre(Graph original_g, ArrayList<Graph> sanitized_gs) {
		final BitSet[] ground_truth = original_g.get_adjancency_matrix_as_bit_vector();
		ArrayList<Integer>[] neighbors = original_g.get_neighbors();
		double[] all_mre = new double[original_g.num_vertices];
		
		for(int i=0;i<sanitized_gs.size();i++) {
			Graph g = sanitized_gs.get(i);
			final BitSet[] to_ceck = g.get_adjancency_matrix_as_bit_vector();
			final int[][] result = edge_edit_dist(ground_truth, to_ceck);
			
			double sum = 0;
			for(int node=0;node<original_g.num_vertices;node++) {
				int[] r = result[node];
				double mae_node = r[missing_edge] + r[new_fake_edge];
				double mre_node = mae_node / (double)neighbors[node].size();//normalize by original node count
				sum+= mre_node;
			}
			
			double mae = sum / original_g.num_vertices; //normalize by Graph size (i.e., |V|)
			all_mre[i] = mae;
		}

		return all_mre;
	}
	
	static final double[] page_rank(final Graph original_g, final Graph sanitized_g) {
		ArrayList<Graph> g_s = new ArrayList<Graph>();
		g_s.add(sanitized_g);
		return page_rank(original_g, g_s);
	}
	
	/**
	 * returns average page_rank delta for each sanitized_gs to original_g
	 */
	public static final double[] page_rank(final Graph original_g, final ArrayList<Graph> sanitized_gs) {
		final double[] page_rank_org = PageRank.run(original_g);
		final double[] result = new double[sanitized_gs.size()];
		for(int i=0;i<sanitized_gs.size();i++) {
			final double[] page_rank_san_i = PageRank.run(sanitized_gs.get(i));
			double delta = abs_avg_delta(page_rank_org, page_rank_san_i);
			result[i] = delta;
		}
		return result;
	}
	
	/**
	 * returns average page_rank delta for each sanitized_gs to original_g
	 */
	public static final double page_rank(final Graph g) {
		final double[] page_rank = PageRank.run(g);
		double avg_page_rank = avg(page_rank);
		return avg_page_rank;
	}
	public static final double[] page_rank(final ArrayList<Graph> all_g) {
		final double[] all_page_rank = new double[all_g.size()];
		for(int i=0;i<all_g.size();i++){
			double avg_page_rank = avg(PageRank.run(all_g.get(i)));
			all_page_rank[i] = avg_page_rank;
		}
		return all_page_rank;
	}

	private static double[] avg(double[][] temp_results) {
		int size = temp_results.length;
		int width = temp_results[0].length;
		double[] results = new double[width]; 
		
		for(double[] arr : temp_results) {
			for(int i=0;i<results.length;i++) {
				results[i] += arr[i];
			}
		}
		
		for(int i=0;i<results.length;i++) {
			results[i] /= (double)size;
		}
		
		return results;
	}
	
	private static double[] avg(ArrayList<double[]> temp_results) {
		int size = temp_results.size();
		int width = temp_results.get(0).length;
		double[] results = new double[width]; 
		
		for(double[] arr : temp_results) {
			for(int i=0;i<results.length;i++) {
				results[i] += arr[i];
			}
		}
		
		for(int i=0;i<results.length;i++) {
			results[i] /= (double)size;
		}
		
		return results;
	}

	private static double[] abs_avg_delta(double[][] pp_org, double[][] pp_san_i) {
		if(pp_org.length!=pp_san_i.length) {
			System.err.println("pp_org.length!=pp_san_i.length");
		}
		double[] agg_results = new double[pp_org[0].length];
		for(int node=0;node<pp_org.length;node++) {
			double[] org_r = pp_org[node];
			double[] san_r = pp_san_i[node];
			for(int i=0;i<agg_results.length;i++) {
				agg_results[i] += Math.abs(org_r[i]-san_r[i]);
			}
		}
		double size = pp_org.length;
		for(int i=0;i<agg_results.length;i++) {
			agg_results[i] /= size;
		}
		return agg_results;
	}

	private static double abs_avg_delta(double[] page_rank_org, double[] page_rank_san) {
		if(page_rank_org.length!=page_rank_san.length) {
			System.err.println("page_rank_org.length!=page_rank_san.length");
		}
		double delta = 0;
		for(int i=0;i<page_rank_org.length;i++) {
			delta+=Math.abs(page_rank_org[i]-page_rank_san[i]);
		}
		return delta/(double)page_rank_org.length;
	}

	public static double avg(double[] arr) {
		double sum = sum(arr);
		return sum/(double)arr.length;
	}

	private static double sum(double[] arr) {
		double sum = 0;
		for(double  d: arr) {
			sum+=d;
		}
		return sum;
	}
	
	
	public static double[] count_triangles(Graph org_g, ArrayList<Graph> private_gs) {
		System.out.println("count_triangles(Graph, ArrayList<Graph>)");
		double start = System.currentTimeMillis();
		final double[] count = new double[private_gs.size()];
		final long true_count_triangles = count_triangles(org_g);
		
		for(int i=0;i<private_gs.size();i++) {
			count[i] = count_triangles(private_gs.get(i));
			//deltas[i] = Math.abs(true_count_triangles-my_count);
		}
		double avg = avg(count);
		double[] ret = {true_count_triangles, avg};
		System.out.println("count_triangles(Graph, ArrayList<Graph>) [Done] in "+(System.currentTimeMillis()-start)+" ms");
		return ret;
	}
	
	
	//FIXME
	public static double[] count_triangle_delta(Graph g, ArrayList<Graph> private_gs) {
		//System.out.println("count_triangles(Graph)");
		double start = System.currentTimeMillis();

		final double[] count = new double[private_gs.size()];
		ArrayList<int[]> all_org_triangles = get_triangles(g);
		for(int i=0;i<private_gs.size();i++) {
			count[i] = count_triangle_delta(all_org_triangles, private_gs.get(i));
		}
		double avg = avg(count);
		double[] ret = {all_org_triangles.size(), avg};
		System.out.println("count_triangle_delta(Graph, ArrayList<int[]>) [Done] in " + (System.currentTimeMillis() - start) + " ms");
		return ret;
	}
	
	private static long count_triangle_delta(ArrayList<int[]> all_org_triangles, Graph g) {//XXX Java...
		//System.out.println("count_triangles(Graph)");
		double start = System.currentTimeMillis();
		/**
		 * Adjacency matrix of g
		 */
		final BitSet[] am = g.get_adjancency_matrix_as_bit_vector();
		long count_triangles = 0;
		for(int[] org_triangle : all_org_triangles) {
			final int i = org_triangle[0];
			final int j = org_triangle[1];
			final int k = org_triangle[2];
			if(am[i].get(j) && am[j].get(k) && am[k].get(i)) {
				count_triangles++;
			}
		}
		
		System.out.println("count_triangle_delta(Graph, ArrayList<int[]>) [Done] in " + (System.currentTimeMillis() - start) + " ms");
		return count_triangles;
	}
	
	private static ArrayList<int[]> get_triangles(Graph g) {
		//System.out.println("count_triangles(Graph)");
		double start = System.currentTimeMillis();
		HashMap<Integer, ArrayList<int[]>> all_triangles = new HashMap<Integer, ArrayList<int[]>>();
		/**
		 * Adjacency matrix of g
		 */
		final BitSet[] am = g.get_adjancency_matrix_as_bit_vector();
		final int V = g.num_vertices;
		for (int i = 0; i < V; i++) {
			for (int j = 0; j < V; j++) {
				if(am[i].get(j)) {
					for (int k = 0; k < V; k++) {
						if (am[j].get(k) && am[k].get(i)) {// This is a triangle
							int[] triangle = {i,j,k};
							Arrays.sort(triangle);
							insert_if_not_contained(all_triangles, triangle);
						}
					}
				}
			}
		}

		ArrayList<int[]> result = new ArrayList<int[]>();
		for(Entry<Integer, ArrayList<int[]>> hash_bucket : all_triangles.entrySet()) {
			for(int[] triangle : hash_bucket.getValue()) {
				result.add(triangle);
			}
		}
		
		System.out.println("get_triangles(Graph) [Done] in " + (System.currentTimeMillis() - start) + " ms");
		return result;
	}
	
	/**
	 * The problem is, we find each triangle three times. 
	 * So we must efficiently check, whether we got this one already. We do this by hashing the triangles.
	 * 
	 * @param all_triangles
	 * @param triangle
	 * @return
	 */
	private static boolean insert_if_not_contained(HashMap<Integer, ArrayList<int[]>> all_triangles, int[] triangle) {
		int hash_code = triangle[0] +triangle[1] +triangle[2];
		ArrayList<int[]> my_list = all_triangles.get(hash_code);
		if(my_list == null) {
			my_list = new ArrayList<int[]>();
			all_triangles.put(hash_code, my_list);
		}
		for(int[] other_triangle : my_list) {
			if(other_triangle[0] == triangle[0] && other_triangle[1] == triangle[1] && other_triangle[2] == triangle[2]) {
				return false;
			}
		}
		my_list.add(triangle);
		return true;
	}

	public static long count_triangles(Graph g) {
		//System.out.println("count_triangles(Graph)");
		double start = System.currentTimeMillis();
		/**
		 * Adjacency matrix of g
		 */
		final BitSet[] am = g.get_adjancency_matrix_as_bit_vector();
		long count_triangles = 0;
		final int V = g.num_vertices;
		for (int i = 0; i < V; i++) {
			for (int j = 0; j < V; j++) {
				if(am[i].get(j)) {
					for (int k = 0; k < V; k++) {
						if (am[j].get(k) && am[k].get(i)) {// This is a triangle
							count_triangles++;
						}
					}
				}
			}
		}

		// We find every triangle 3 times in directed graph (6 times in an undirected)
		count_triangles /= 3;
		System.out.println("count_triangles(Graph) [Done] in " + (System.currentTimeMillis() - start) + " ms");
		return count_triangles;
	}
	
	public static double[] get_array_statistics(double expected_error, double[] observed_error) {
		// The mean average
		double mean = 0.0;
		for (int i = 0; i < observed_error.length; i++) {
		        mean += observed_error[i];
		}
		mean /= observed_error.length;

		// The variance
		double variance = 0;
		for (int i = 0; i < observed_error.length; i++) {
		    variance += Math.pow(observed_error[i] - mean, 2);
		}
		variance /= observed_error.length;

		// Standard Deviation
		double std = Math.sqrt(variance);
		
		double[] ret = {mean,variance,std}; 
		return ret; 
	}
}