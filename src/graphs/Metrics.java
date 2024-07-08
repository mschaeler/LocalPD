package graphs;

import java.util.ArrayList;
import java.util.Arrays;

import algorithms.PageRank;
import algorithms.ProximityPrestige;
import algorithms.ProximityPrestige.InfDom;
import results.ResultCollector;

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

	public static final double[] statistics(Graph original_g, ArrayList<Graph> sanitized_gs) {
		final boolean[][] ground_truth = original_g.get_adjancency_matrix_as_bit_vector();
		ArrayList<double[]> all_sums = new ArrayList<double[]>();
		
		for(int i=0;i<sanitized_gs.size();i++) {
			double[] sum = new double[4];
			Graph g = sanitized_gs.get(i);
			final boolean[][] to_ceck = g.get_adjancency_matrix_as_bit_vector();
			final int[][] result = node_edit_dist(ground_truth, to_ceck);
			
			for(int[] r : result) {
				sum[correct_edge]+= r[correct_edge];
				sum[missing_edge]+= r[missing_edge];
				sum[new_fake_edge]+= r[new_fake_edge];
				sum[correct_none_edge]+= r[correct_none_edge];
			}
			
			sum[correct_edge]/= original_g.num_vertices;
			sum[missing_edge]/= original_g.num_vertices;
			sum[new_fake_edge]/= original_g.num_vertices;
			sum[correct_none_edge]/= original_g.num_vertices;
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
		final boolean[][] ground_truth = original_g.get_adjancency_matrix_as_bit_vector();
		final boolean[][] to_ceck = sanitized_g.get_adjancency_matrix_as_bit_vector();
		final int[][] result = node_edit_dist(ground_truth, to_ceck);
		
		double sum = 0;
		for(int[] r : result) {
			sum+= r[missing_edge];
			sum+= r[new_fake_edge];
		}
		
		double mae = sum / original_g.num_vertices; //normalize by Graph size (i.e., |V|)
		
		return mae;
	}
	
	static final double[] mae(Graph original_g, ArrayList<Graph> sanitized_gs) {
		final boolean[][] ground_truth = original_g.get_adjancency_matrix_as_bit_vector();
		double[] all_mae = new double[original_g.num_vertices];
		
		for(int i=0;i<sanitized_gs.size();i++) {
			Graph g = sanitized_gs.get(i);
			final boolean[][] to_ceck = g.get_adjancency_matrix_as_bit_vector();
			final int[][] result = node_edit_dist(ground_truth, to_ceck);
			
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
	
	
	private static int[][] node_edit_dist(final boolean[][] ground_truth, final boolean[][] to_ceck) {
		if(ground_truth.length!=to_ceck.length) {
			System.err.println("node_edit_dist() ground_truth.length!=to_ceck.length");
		}
		final int[][] result = new int[ground_truth.length][4];
		for(int node = 0;node<ground_truth.length;node++) {
			for(int e = 0;e<ground_truth.length;e++) {
				if(ground_truth[node][e] && to_ceck[node][e]) {
					result[node][correct_edge]++;
				}else if(ground_truth[node][e] && !to_ceck[node][e]){//there should be an edged, but there is none
					result[node][missing_edge]++;
				}else if(!ground_truth[node][e] && to_ceck[node][e]){//we invented a new edge
					result[node][new_fake_edge]++;
				}else{//in both cases there is no edge
					result[node][correct_none_edge]++;
				}
			}
		}
		return result;
	}

	static final double mre(Graph original_g, Graph sanitized_g) {		
		final boolean[][] ground_truth = original_g.get_adjancency_matrix_as_bit_vector();
		ArrayList<Integer>[] neighbors = original_g.get_neighbors();
		final boolean[][] to_ceck = sanitized_g.get_adjancency_matrix_as_bit_vector();
		final int[][] result = node_edit_dist(ground_truth, to_ceck);
		
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
		final boolean[][] ground_truth = original_g.get_adjancency_matrix_as_bit_vector();
		ArrayList<Integer>[] neighbors = original_g.get_neighbors();
		double[] all_mre = new double[original_g.num_vertices];
		
		for(int i=0;i<sanitized_gs.size();i++) {
			Graph g = sanitized_gs.get(i);
			final boolean[][] to_ceck = g.get_adjancency_matrix_as_bit_vector();
			final int[][] result = node_edit_dist(ground_truth, to_ceck);
			
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
	 * returns average proximity prestige delta for each sanitized_gs to original_g
	 */
	public static final double[] proximity_prestige(final Graph original_g, final ArrayList<Graph> sanitized_gs) {
		final double[][] pp_org = ProximityPrestige.run(original_g);
		ArrayList<double[]> temp_results = new ArrayList<double[]>(sanitized_gs.size());
		
		for(int i=0;i<sanitized_gs.size();i++) {
			final double[][] pp_san_i = ProximityPrestige.run(sanitized_gs.get(i));
			double[] deltas = abs_avg_delta(pp_org, pp_san_i);//[I,avg dist(),pp]
			temp_results.add(deltas);
		}
		
		final double[] result = avg(temp_results);
		
		return result;
	}
	/**
	 * returns average proximity per sanitized graph
	 */
	public static final double[] proximity_prestige(final ArrayList<Graph> sanitized_gs) {
		ArrayList<double[]> temp_results = new ArrayList<double[]>(sanitized_gs.size());
		
		for(int i=0;i<sanitized_gs.size();i++) {
			final double[][] pp_san_i = ProximityPrestige.run(sanitized_gs.get(i));
			double[] arr = avg(pp_san_i);//[I,avg dist(),pp]
			temp_results.add(arr);
		}
		
		final double[] result = avg(temp_results);
		
		return result;
	}
	
	/**
	 * returns average proximity prestige delta for each sanitized_gs to original_g
	 */
	public static final double[] proximity_prestige(final Graph original_g) {
		ResultCollector.init();
		final double[][] pp_org = ProximityPrestige.run(original_g);
		final double[] result = avg(pp_org);
		ResultCollector.all_proximity_prestige.add(result);
		String name = "org g";
		double[] all_eps = {-1};
		ResultCollector.store(name, all_eps);
		
		return result;
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
}