package graphs;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Random;
import java.util.stream.IntStream;

import algorithms.DK_Series;
import algorithms.LaplaceStream;
import misc.KMedoids;
import misc.KMedoids.Clustering;
import misc.TopK;
import misc.Util;
import results.Config;
import results.TheoreticalAnalysis;

public class Mechanism {
	public static java.util.Random rand = new java.util.Random(Util.seed);
	
	public static boolean use_safe_variant = false;//usually much fast. sufficient for benchmarking, but not suited for productive use
	
	/**
	 * Determines whether counts are rounded to the next integer > 0
	 */
	private static boolean TRUNCATE = true;

	/**
	 * Creates a edge-local DP version of the input Graph g. We assume that g is directed.
	 * @param g - non-private graph as input
	 * @param p - frac{1}{1+exp(\epsilon)}
	 * @return
	 */
	public static Graph radomized_response(Graph g, final double epsilon) {
		if(epsilon<3) {
			return null;//XXX gets to big otherwise
		}
		System.out.println("radomized_response()"+(Config.materialize_graph ? " materialize_graph" : "")+" e="+epsilon);
		double start = System.currentTimeMillis();
		
		final double p = epsilon_to_p(epsilon);
		ArrayList<Integer>[] neighbor_list = g.get_neighbors();
		Graph private_g = new Graph(g.num_vertices, g.name+"_rr");
		
		final boolean[] adjacency_matrix_line = new boolean[g.num_vertices];
		//System.out.println("terror_rr\tmean\tvar\tstd");
		//double error_rr = TheoreticalAnalysis.error_rr(neighbor_list[0].size(), epsilon, g.num_vertices)[0];
		//double[] observed_errors = new double[g.num_vertices];
		
		for(int node=0;node<g.num_vertices;node++) {
			//ArrayList<Integer> my_sanitized_neighbors = new ArrayList<Integer>();//XXX just for debugging
			Arrays.fill(adjacency_matrix_line, false);
			for(int target_node : neighbor_list[node]) {
				adjacency_matrix_line[target_node] = true;
			}
			for(int target_node=0;target_node<g.num_vertices;target_node++) {
				double dice = rand.nextDouble();
				boolean create_edge;
				if(dice <= (1-p)) {//with probability 1-p return the real value
					create_edge = adjacency_matrix_line[target_node];
				}else{//else flip the bit
					create_edge = (adjacency_matrix_line[target_node]) ? false : true;
				}
				if(create_edge){
					//my_sanitized_neighbors.add(target_node);//XXX just for debugging
					if(Config.materialize_graph) {
						Graph.write_edge(node, target_node);
					}else{
						private_g.add_edge(node, target_node);
					}
				}
			}
				
			if(node%1000==0) {
				System.out.print("node="+node+" ");
			}
			
			
			//int[] errors = Metrics.edge_edit_dist(neighbor_list[node],my_sanitized_neighbors, g.num_vertices);
			//double observed_error = errors[1]+errors[2];
			//observed_errors[node] = observed_error;
			//System.out.println(node+"\t"+neighbor_list[node].size()+"\t"+error_rr+"\t"+observed_error);
		}
		//double[] statistics = Metrics.get_array_statistics(error_rr, observed_errors);
		//System.out.println(error_rr+"\t"+statistics[0]+"\t"+statistics[1]+"\t"+statistics[2]);
		Graph.stop(start, "radomized_response()");
		
		return private_g;
	}
	
	static final double expected_error_rr(final double num_vertices, final double sanitized_num_edges, final double p) {
		final double correct_edges_flipped = sanitized_num_edges * p;//Second case in RR mechanism
		final double non_edges_flipped = (num_vertices - sanitized_num_edges) * (1.0d-p);
		return correct_edges_flipped + non_edges_flipped;
	}
	
	static final double count_query_senstivity = 1;
	/**
	 * Version with non-private grouping
	 * @param g
	 * @param p
	 * @return
	 */
	static Graph k_edge_radomized_response_non_private_grouping(Graph g, final double all_epsilon, final double epsilon_q1, boolean use_seq_composition) {
		String name = "k_edge_radomized_response_non_private_grouping() use_seq_composition="+use_seq_composition;
		System.out.println(name);
		double start = System.currentTimeMillis();
		ArrayList<Integer>[] neighbor_list = g.get_neighbors();
		
		//distribute budget
		final double count_query_epsilon = epsilon_q1;
		final double epsilon_q2 = all_epsilon-epsilon_q1;
		
		@SuppressWarnings("unchecked")
		ArrayList<Integer>[] sanitized_neighbor_list = new ArrayList[neighbor_list.length]; 
		
		for(int node=0;node<g.num_vertices;node++) {
			//Execute Q1
			int sanitized_number_of_outgoing_egdes = q_1(neighbor_list[node], count_query_senstivity, count_query_epsilon);
			
			//Partition and Shuffle the data
			ArrayList<Integer> my_neighbors = neighbor_list[node]; 
			Collections.shuffle(my_neighbors);//XXX avoid that the position leaks info whether this is rather true or invented edge, i.e., this is means all submechanisms perform a sequential composition.
			
			boolean[] true_edges = new boolean[g.num_vertices];
			for(int my_edge : my_neighbors) {
				true_edges[my_edge] = true;// Note here I access the true data. Thus the partitioning scheme is not private
			}
			ArrayList<Integer> fake_edges = new ArrayList<Integer>(g.num_vertices);
			for(int target_vertex=0;target_vertex<g.num_vertices;target_vertex++) {
				if(!true_edges[target_vertex]) {//can be used as fake edge
					fake_edges.add(target_vertex);
				}
			}
			Collections.shuffle(fake_edges);//also ensure uniform picks of fake edges
			if(my_neighbors.size()+fake_edges.size()!=g.num_vertices){
				System.err.println("my_neighbors.size()+fake_edges.size()!=g.num_vertices");
			}
			
			//Execute Q2
			ArrayList<Integer> my_sanitized_neighbors;
			if(!use_seq_composition) {
				my_sanitized_neighbors = q_2(my_neighbors, fake_edges, sanitized_number_of_outgoing_egdes, epsilon_q2, rand);	
			}else{
				my_sanitized_neighbors = q_2_seq_composition(my_neighbors, fake_edges, sanitized_number_of_outgoing_egdes, epsilon_q2, rand);
			}
			
			sanitized_neighbor_list[node] = my_sanitized_neighbors;
			
			if(node%10000==0) {
				System.out.println("node="+node);
			}
		}
		Graph.stop(start, "k_edge_radomized_response_non_private_grouping()");
		return new Graph(sanitized_neighbor_list,g.name+" k_edge_radomized_response_non_private_grouping() use_seq_composition="+use_seq_composition);
	}
	
	static Graph k_edge_radomized_response_non_private_grouping_rr_fallback(Graph g, final double all_epsilon, final double epsilon_q1, boolean use_seq_composition) {
		String name = "k_edge_radomized_response_non_private_grouping_rr_fallback() use_seq_composition="+use_seq_composition;
		System.out.println(name);
		double start = System.currentTimeMillis();
		ArrayList<Integer>[] neighbor_list = g.get_neighbors();
		
		//distribute budget
		final double count_query_epsilon = epsilon_q1;
		final double epsilon_q2 = all_epsilon-epsilon_q1;
		final double p = epsilon_to_p(epsilon_q2);
		
		final boolean[] adjacency_matrix_line = new boolean[g.num_vertices];
		final boolean[] sanitized_adjacency_matrix_line = new boolean[g.num_vertices];
		
		@SuppressWarnings("unchecked")
		ArrayList<Integer>[] sanitized_neighbor_list = new ArrayList[neighbor_list.length]; 
		
		for(int node=0;node<g.num_vertices;node++) {
			ArrayList<Integer> my_sanitized_neighbors;
			
			//Execute Q1
			int sanitized_number_of_outgoing_egdes = q_1(neighbor_list[node], count_query_senstivity, count_query_epsilon);
			
			//Decide whether we use k-edge or fall back to randomized response
			double error_rr = TheoreticalAnalysis.error_rr(sanitized_number_of_outgoing_egdes, epsilon_q2, g.num_vertices)[0];
			double error_k_edge;  
			if(use_seq_composition) {//XXX assumes epsilon_q1 = 1
				error_k_edge = TheoreticalAnalysis.error_k_edge_seq_comp(sanitized_number_of_outgoing_egdes, all_epsilon)[2];
			}else {
				error_k_edge = TheoreticalAnalysis.error_k_edge_non_private(sanitized_number_of_outgoing_egdes, all_epsilon)[2];
			}
			
			if(error_rr < error_k_edge) {// Use randomized response
				//System.out.println("Fall back");
				my_sanitized_neighbors = new ArrayList<Integer>();
				Arrays.fill(adjacency_matrix_line, false);
				Arrays.fill(sanitized_adjacency_matrix_line, false);
				for(int target_node : neighbor_list[node]) {
					adjacency_matrix_line[target_node] = true;
				}
				for(int target_node=0;target_node<g.num_vertices;target_node++) {
					double dice = rand.nextDouble();
					boolean create_edge;
					if(dice <= (1-p)) {//with probability 1-p return the real value
						create_edge = adjacency_matrix_line[target_node];
					}else {//else flip the bit
						create_edge = (adjacency_matrix_line[target_node]) ? false : true;
					}
					if(create_edge) {
						my_sanitized_neighbors.add(target_node);
					}
				}
			}else{//Use k-edge-return
				//Partition and Shuffle the data
				ArrayList<Integer> my_neighbors = neighbor_list[node]; 
				Collections.shuffle(my_neighbors);//XXX avoid that the position leaks info whether this is rather a true or invented edge, i.e., this is means all submechanisms perform a sequential composition.
				
				boolean[] true_edges = new boolean[g.num_vertices];
				for(int my_edge : my_neighbors) {
					true_edges[my_edge] = true;// Note here I access the true data. Thus the partitioning scheme is not private
				}
				ArrayList<Integer> fake_edges = new ArrayList<Integer>(g.num_vertices);
				for(int target_vertex=0;target_vertex<g.num_vertices;target_vertex++) {
					if(!true_edges[target_vertex]) {//can be used as fake edge
						fake_edges.add(target_vertex);
					}
				}
				Collections.shuffle(fake_edges);//also ensure uniform picks of fake edges
				if(my_neighbors.size()+fake_edges.size()!=g.num_vertices){
					System.err.println("my_neighbors.size()+fake_edges.size()!=g.num_vertices");
				}
				
				//Execute Q2
				if(!use_seq_composition) {
					my_sanitized_neighbors = q_2(my_neighbors, fake_edges, sanitized_number_of_outgoing_egdes, epsilon_q2, rand);	
				}else{
					my_sanitized_neighbors = q_2_seq_composition(my_neighbors, fake_edges, sanitized_number_of_outgoing_egdes, epsilon_q2, rand);
				}
			}
			sanitized_neighbor_list[node] = my_sanitized_neighbors;
			
			if(node%10000==0) {
				System.out.println("node="+node);
			}
		}
		Graph.stop(start, "k_edge_radomized_response_non_private_grouping_rr_fallback()");
		return new Graph(sanitized_neighbor_list,g.name+" k_edge_radomized_response_non_private_grouping() use_seq_composition="+use_seq_composition);
	}

	static ArrayList<Integer> rr_fall_back(final int node, final double p, ArrayList<Integer>[] neighbor_list, Graph g, Random rand) {
		//System.out.println("Fall back");
		final boolean[] adjacency_matrix_line = new boolean[g.num_vertices];
		final boolean[] sanitized_adjacency_matrix_line = new boolean[g.num_vertices];
		
		ArrayList<Integer> my_sanitized_neighbors = new ArrayList<Integer>();
		Arrays.fill(adjacency_matrix_line, false);
		Arrays.fill(sanitized_adjacency_matrix_line, false);
		for(int target_node : neighbor_list[node]) {
			adjacency_matrix_line[target_node] = true;
		}
		for(int target_node=0;target_node<g.num_vertices;target_node++) {
			double dice = rand.nextDouble();
			boolean create_edge;
			if(dice <= (1-p)) {//with probability 1-p return the real value
				create_edge = adjacency_matrix_line[target_node];
			}else {//else flip the bit
				create_edge = (adjacency_matrix_line[target_node]) ? false : true;
			}
			if(create_edge) {
				my_sanitized_neighbors.add(target_node);
			}
		}
		return my_sanitized_neighbors;
	}
	
	public static Graph m_straw_man(Graph g, final double all_epsilon) {
		String name = "m_straw_man() e="+all_epsilon+" e_1="+all_epsilon;
		System.out.println(name);
		double start = System.currentTimeMillis();
		ArrayList<Integer>[] neighbor_list = g.get_neighbors();
		
		//distribute budget
		final double count_query_epsilon = all_epsilon;
		
		@SuppressWarnings("unchecked")
		ArrayList<Integer>[] sanitized_neighbor_list = new ArrayList[neighbor_list.length]; 
		
		//boolean[] duplicate_check = new boolean[g.num_vertices];
		
		for(int node=0;node<g.num_vertices;node++) {
			//Arrays.fill(duplicate_check, false);

			//Execute Q1
			int sanitized_number_of_outgoing_egdes = q_1(neighbor_list[node], count_query_senstivity, count_query_epsilon);
			//System.out.println(sanitized_number_of_outgoing_egdes+"\t"+neighbor_list[node].size()+"\t"+count_query_epsilon);
			ArrayList<Integer> my_sanitized_neighbors = new ArrayList<Integer>();
			
			if(sanitized_number_of_outgoing_egdes<=0) {
				//do nothing
			}else{
				while(my_sanitized_neighbors.size()<sanitized_number_of_outgoing_egdes) {
					my_sanitized_neighbors.add(rand.nextInt(g.num_vertices));
				}
			}
			if(Config.materialize_graph) {
				Graph.write_edge_list(node, my_sanitized_neighbors);
			}else{
				sanitized_neighbor_list[node] = my_sanitized_neighbors;	
			}
			
			if(node%1000==0) {
				System.out.print("node="+node+" ");
			}
		}
		Graph.stop(start, "m_straw_man()");
		if(Config.materialize_graph) {
			return null;
		}else{
			return new Graph(sanitized_neighbor_list, g.name+" "+name);	
		}
	}
	
	public static Graph m_sample(Graph g, final double all_epsilon, final double epsilon_q1, boolean no_rr_fall_back) {
		String name = "m_sample() e="+all_epsilon+" e_1="+epsilon_q1+(no_rr_fall_back ? "disable_rr_fall_back" : "");
		System.out.println(name);
		double start = System.currentTimeMillis();
			
		ArrayList<Integer>[] neighbor_list = g.get_neighbors();
		
		//distribute budget
		final double count_query_epsilon = epsilon_q1;
		final double epsilon_q2 = all_epsilon-epsilon_q1;
		final double p = epsilon_to_p(epsilon_q2);
		
		if(TheoreticalAnalysis.error_rr(1, all_epsilon, g.num_vertices)[0]< TheoreticalAnalysis.error_m_part(1, epsilon_q1, epsilon_q2, g.num_vertices)) {
			return radomized_response(g, all_epsilon);
		}
		
		@SuppressWarnings("unchecked")
		ArrayList<Integer>[] sanitized_neighbor_list = new ArrayList[neighbor_list.length]; 
		
		boolean[] duplicate_check = new boolean[g.num_vertices];
		
		System.out.println("node\t|N_s|\t|N|\tError rr\terror_m_part\terror_sample\tobserved error");
		
		for(int node=0;node<g.num_vertices;node++) {
			Arrays.fill(duplicate_check, false);

			//Execute Q1
			int sanitized_number_of_outgoing_egdes = q_1(neighbor_list[node], count_query_senstivity, count_query_epsilon);
			sanitized_number_of_outgoing_egdes = neighbor_list[node].size();
			
			//Decide whether we use k-edge or fall back to randomized response
			double error_rr = TheoreticalAnalysis.error_rr(sanitized_number_of_outgoing_egdes, epsilon_q2, g.num_vertices)[0];
			double error_m_part = TheoreticalAnalysis.error_m_part(sanitized_number_of_outgoing_egdes, epsilon_q1, epsilon_q2, g.num_vertices);
			double error_sample = TheoreticalAnalysis.error_m_sample(sanitized_number_of_outgoing_egdes, epsilon_q1, epsilon_q2, g.num_vertices);
			//System.out.println("|N_s|="+sanitized_number_of_outgoing_egdes+"|N|="+neighbor_list[node].size()+" Errors rr="+error_rr+" k_edge="+error_k_edge+" error_m_part="+error_m_part+" error_sample="+error_sample);
			
			if(no_rr_fall_back) {
				error_sample = 0;
			}
			ArrayList<Integer> my_sanitized_neighbors;
			
			
			if(sanitized_number_of_outgoing_egdes<=0) {
				my_sanitized_neighbors = new ArrayList<Integer>();
			}else if(error_rr<error_sample){
				//rr-fallback
				my_sanitized_neighbors = rr_fall_back(node, p, neighbor_list, g, rand);
			}else{
				ArrayList<Integer> my_neighbors = neighbor_list[node];
				my_sanitized_neighbors = new ArrayList<Integer>();
				final double my_scores[] = new double[g.num_vertices];
				for(int edge : my_neighbors) {
					my_scores[edge] = Sample.SCORE_EXISTING_EDGE;//for all other edges it remains zero
				}
				for(int edge=0;edge<my_scores.length;edge++) {
					my_scores[edge] = Math.pow(Math.E,epsilon_q2/(double)sanitized_number_of_outgoing_egdes*my_scores[edge]/(2.0d*Sample.delta_u()));
				}
				HashSet<Integer> edges = new HashSet<Integer>();//We do not want duplicates
				while(edges.size()<sanitized_number_of_outgoing_egdes) {
					edges.add(sample(my_scores, rand));
				}
				for(int edge : edges) {
					my_sanitized_neighbors.add(edge);
				}
			}
			
			int[] erorrs = Metrics.edge_edit_dist(neighbor_list[node],my_sanitized_neighbors, g.num_vertices);
			double observed_error = erorrs[1]+erorrs[2];
			System.out.println(node+"\t"+sanitized_number_of_outgoing_egdes+"\t"+neighbor_list[node].size()+"\t"+error_rr+"\t"+error_m_part+"\t"+error_sample+"\t"+observed_error);
			
			if(Config.materialize_graph) {
				Graph.write_edge_list(node, my_sanitized_neighbors);
			}else {
				sanitized_neighbor_list[node] = my_sanitized_neighbors;
			}
			
			if(node%1000==0) {
				System.out.print("node="+node+" ");
			}
			
			if(node%10000==0) {
				System.out.println("node="+node);
			}
		}
		Graph.stop(start, "m_sample()");
		if(Config.materialize_graph)
			return null;
		return new Graph(sanitized_neighbor_list, g.name+" "+name);
	}
	
	static int counter = 0;
	public static Graph m_sample_weighted(Graph g, final double all_epsilon, final double epsilon_q1, boolean no_rr_fall_back) {
		String name = "m_sample_weighted() e="+all_epsilon+" e_1="+epsilon_q1+(no_rr_fall_back ? " disable_rr_fall_back" : "");
		System.out.println(name);
		double start = System.currentTimeMillis();
		ArrayList<Integer>[] neighbor_list = g.get_neighbors();
		
		//distribute budget
		final double count_query_epsilon = epsilon_q1;
		final double epsilon_q2 = all_epsilon-epsilon_q1;
		final double p = epsilon_to_p(epsilon_q2);
		
		@SuppressWarnings("unchecked")
		ArrayList<Integer>[] sanitized_neighbor_list = new ArrayList[neighbor_list.length]; 
		
		boolean[] duplicate_check = new boolean[g.num_vertices];
		double sum_error_expected = 0;
		double sum_error_observed = 0;
		/**
		 * Frequency of each node as target
		 * XXX make private
		 */
		final double[] histogram   = new double[g.num_vertices];
		
		for(int node=0;node<g.num_vertices;node++) {//TODO locally DP?
			final ArrayList<Integer> my_true_neighbors = neighbor_list[node];
			for(int neighbor : my_true_neighbors) {
				histogram[neighbor]++;
			}
		}
		to_probabilities(histogram);
		
		for(int node=0;node<g.num_vertices;node++) {
			Arrays.fill(duplicate_check, false);

			//Execute Q1
			int sanitized_number_of_outgoing_egdes = q_1(neighbor_list[node], count_query_senstivity, count_query_epsilon);
			
			//Decide whether we use k-edge or fall back to randomized response
			double error_rr = TheoreticalAnalysis.error_rr(sanitized_number_of_outgoing_egdes, epsilon_q2, g.num_vertices)[0];
			//double error_k_edge = TheoreticalAnalysis.error_k_edge_random_part(sanitized_number_of_outgoing_egdes, all_epsilon)[0];//TODO
			double error_sample = TheoreticalAnalysis.error_m_sample(sanitized_number_of_outgoing_egdes, epsilon_q1, epsilon_q2, g.num_vertices);
			if(no_rr_fall_back) {
				error_sample = 0;
			}
			ArrayList<Integer> my_sanitized_neighbors;
			
			
			if(sanitized_number_of_outgoing_egdes<=0) {
				my_sanitized_neighbors = new ArrayList<Integer>();
			}else if(error_rr<error_sample){
				//rr-fallback
				my_sanitized_neighbors = rr_fall_back(node, p, neighbor_list, g, rand);
			}else{
				ArrayList<Integer> my_neighbors = neighbor_list[node];
				my_sanitized_neighbors = new ArrayList<Integer>();
				final double my_scores[] = new double[g.num_vertices];
				for(int edge : my_neighbors) {
					my_scores[edge] = Sample.SCORE_EXISTING_EDGE;//for all other edges it remains zero
				}
				for(int edge=0;edge<my_scores.length;edge++) {//TODO min(c_s,|N|)
					my_scores[edge] = Math.pow(Math.E,epsilon_q2/(double)sanitized_number_of_outgoing_egdes*my_scores[edge]/(2.0d*Sample.delta_u()));
				}
				to_probabilities(my_scores);
				for(int edge=0;edge<my_scores.length;edge++) {
					my_scores[edge] +=  histogram[edge];//this is the trick
				}
				
				HashSet<Integer> edges = new HashSet<Integer>();//We do not want duplicates
				while(edges.size()<sanitized_number_of_outgoing_egdes) {//Simulates removal of sampled edges
					edges.add(sample(my_scores, rand));
				}

				for(int edge : edges) {
					my_sanitized_neighbors.add(edge);
				}
				//int[] local_error = Metrics.edge_edit_dist(my_neighbors, my_sanitized_neighbors, g.num_vertices);
				//double wrong = local_error[Metrics.missing_edge]+local_error[Metrics.new_fake_edge];
				//if(Math.abs(error_m_part_2-wrong)>10.0d) {
				//	System.out.println("id="+node+" Expecting\t"+error_k_edge+"\tgot\t"+(wrong)+"\terror_rr=\t"+error_rr+"\t"+Arrays.toString(local_error)+"\t|N|="+my_neighbors.size());	
				//}
				//sum_error_expected += error_k_edge;
				//sum_error_observed += wrong;
			}
			
			if(Config.materialize_graph) {
				Graph.write_edge_list(node, my_sanitized_neighbors);
			}else {
				sanitized_neighbor_list[node] = my_sanitized_neighbors;
			}
			
			if(node%1000==0) {
				System.out.print("node="+node+" ");
			}
		}
		System.out.println("sum_error_expected="+sum_error_expected+" sum_error_observed="+sum_error_observed+" |G|="+g.num_vertices);
		Graph.stop(start, "m_sample_weighted()");
		if(Config.materialize_graph)
			return null;
		return new Graph(sanitized_neighbor_list, g.name+" "+name);
	}
	
	static int sample(final double my_scores[], Random rand) {
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
	
	public static Graph m_part_2(Graph g, final double all_epsilon, final double epsilon_q1) {
		String name = "m_part_2() e="+all_epsilon+" e_1="+epsilon_q1;
		System.out.println(name);
		double start = System.currentTimeMillis();
		ArrayList<Integer>[] neighbor_list = g.get_neighbors();
		
		//distribute budget
		final double epsilon_q2 = all_epsilon-epsilon_q1;
		
		@SuppressWarnings("unchecked")
		ArrayList<Integer>[] sanitized_neighbor_list = new ArrayList[neighbor_list.length]; 
		
		for(int node=0;node<g.num_vertices;node++) {
			//Execute Q1
			int sanitized_number_of_outgoing_egdes = q_1(neighbor_list[node], count_query_senstivity, epsilon_q1);
			
			//Decide whether we use k-edge or fall back to randomized response
			double error_rr = TheoreticalAnalysis.error_rr(sanitized_number_of_outgoing_egdes, epsilon_q2, g.num_vertices)[0];
			double error_m_part_2 =-1; //TODO

			ArrayList<Integer> my_sanitized_neighbors;
					
			if(sanitized_number_of_outgoing_egdes<=0) {
				my_sanitized_neighbors = new ArrayList<Integer>();
			}else if(error_rr<error_m_part_2){
				//rr-fallback
				my_sanitized_neighbors = rr_fall_back(node, epsilon_to_p(epsilon_q2), neighbor_list, g, rand);
			}else{				
				my_sanitized_neighbors = new ArrayList<Integer>();
				
				//Partition and Shuffle the data
				if(sanitized_number_of_outgoing_egdes>=g.num_vertices) {
					System.err.println("sanitized_number_of_outgoing_egdes>=g.num_vertices");
				}
				ArrayList<Integer> neighbors = neighbor_list[node];
				Collections.shuffle(neighbors);
				
				//Not all real edges necessarily get a partition. They don't get one iff sanitized_number_of_outgoing_egdes < neighbors.size()
				//Therefore, we shuffled the edges above such that each of them a equally like get no partition, but they move to the pool of non-edges
				ArrayList<Integer> my_neighbors = new ArrayList<Integer>();
				for(int i=0;i<Math.min(sanitized_number_of_outgoing_egdes, neighbors.size());i++) {
					my_neighbors.add(neighbors.get(i));
				}

				boolean[] true_edges = new boolean[g.num_vertices];
				for(int my_edge : my_neighbors) {//Trick: Real edges may become fake ones
					true_edges[my_edge] = true;
				}
				ArrayList<Integer> fake_edges = new ArrayList<Integer>(g.num_vertices);
				for(int target_vertex=0;target_vertex<g.num_vertices;target_vertex++) {
					if(!true_edges[target_vertex]) {//can be used as fake edge
						fake_edges.add(target_vertex);
					}
				}
				Collections.shuffle(fake_edges);//also ensure uniform picks of fake edges
				
				final double a = 1.0d-TheoreticalAnalysis.get_probability_of_selecting_a_fake_node_m_part(g.num_vertices, my_neighbors.size()-1, sanitized_number_of_outgoing_egdes);
				final double probability = (Math.pow(Math.E, epsilon_q2)*a)/(Math.pow(Math.E, epsilon_q2)*a+1);
				if(probability>1) {
					System.err.println("probability>1");
				}
				if(probability<0) {
					System.err.println("probability<0");
				}
				
				//Now the sub mechanisms
				int i=0;
				for(;i<my_neighbors.size();i++) {	
					double dice = rand.nextDouble();
					if(dice<=probability) {//Answer truthfully
						my_sanitized_neighbors.add(my_neighbors.get(i));
					}else{
						my_sanitized_neighbors.add(fake_edges.get(i));//Recap. fake edges are shuffled. Thus, we uniformly draw from them.
					}
				}
				
				while(my_sanitized_neighbors.size()<sanitized_number_of_outgoing_egdes) {
					my_sanitized_neighbors.add(fake_edges.get(i++));
				}
			}
			
			if(Config.materialize_graph) {
				Graph.write_edge_list(node, my_sanitized_neighbors);
			}else {
				sanitized_neighbor_list[node] = my_sanitized_neighbors;
			}
			
			if(node%1000==0) {
				System.out.print("node="+node+" ");
			}
		}
		System.out.println();
		Graph.stop(start, "m_part_2()");
		if(Config.materialize_graph) {
			return null;
		}else{
			return new Graph(sanitized_neighbor_list, g.name+" "+name);	
		}
	}
	
	static Graph m_part(Graph g, final double all_epsilon, final double epsilon_q1) {
		String name = "m_part() e="+all_epsilon+" e_1="+epsilon_q1;
		System.out.println(name);
		double start = System.currentTimeMillis();
		ArrayList<Integer>[] neighbor_list = g.get_neighbors();
		
		//distribute budget
		final double epsilon_q2 = all_epsilon-epsilon_q1;
		
		@SuppressWarnings("unchecked")
		ArrayList<Integer>[] sanitized_neighbor_list = new ArrayList[neighbor_list.length]; 
		int num_rr_falbacks = 0;
		
		double sum_error_expected = 0;
		double sum_error_observed = 0;
		
		//boolean[] duplicate_check = new boolean[g.num_vertices];
		
		for(int node=0;node<g.num_vertices;node++) {
			//Arrays.fill(duplicate_check, false);

			//Execute Q1
			int sanitized_number_of_outgoing_egdes = q_1(neighbor_list[node], count_query_senstivity, epsilon_q1);
			
			//Decide whether we use k-edge or fall back to randomized response
			double error_rr = TheoreticalAnalysis.error_rr(sanitized_number_of_outgoing_egdes, epsilon_q2, g.num_vertices)[0];
			double error_m_part_2 =-1;
			//if(node == 322 && all_epsilon == 2) {
				error_m_part_2 = TheoreticalAnalysis.error_m_part(sanitized_number_of_outgoing_egdes, epsilon_q1, epsilon_q2, g.num_vertices);//TODO	
			//}
			//double error_m_part = TheoreticalAnalysis.error_k_edge_random_part(sanitized_number_of_outgoing_egdes, all_epsilon)[0];//TODO
			ArrayList<Integer> my_sanitized_neighbors;
					
			if(sanitized_number_of_outgoing_egdes<=0) {
				my_sanitized_neighbors = new ArrayList<Integer>();
			}else if(error_rr<error_m_part_2){
				//rr-fallback
				my_sanitized_neighbors = rr_fall_back(node, epsilon_to_p(epsilon_q2), neighbor_list, g, rand);
				num_rr_falbacks++;
				//System.out.println(error_rr+" vs. error_m_part="+error_m_part+" for e="+all_epsilon +" |N|="+neighbor_list[node].size()+" at vertex="+node+" num_rr_falbacks="+num_rr_falbacks);
			}else{
				//Partition and Shuffle the data
				if(sanitized_number_of_outgoing_egdes>=g.num_vertices) {
					System.err.println("sanitized_number_of_outgoing_egdes>=g.num_vertices");
				}
				ArrayList<Integer> neighbors = neighbor_list[node];
				Collections.shuffle(neighbors);
				
				//Not all real edges necessarily get a partition. They don't get one iff sanitized_number_of_outgoing_egdes < neighbors.size()
				//Therefore, we shuffled the edges above such that each of them a equally like get no partition, but they move to the pool of non-edges
				ArrayList<Integer> my_neighbors = new ArrayList<Integer>();
				for(int i=0;i<Math.min(sanitized_number_of_outgoing_egdes, neighbors.size());i++) {
					my_neighbors.add(neighbors.get(i));
				}

				boolean[] true_edges = new boolean[g.num_vertices];
				for(int my_edge : my_neighbors) {//Trick. real edges may become fake ones
					true_edges[my_edge] = true;
				}
				ArrayList<Integer> fake_edges = new ArrayList<Integer>(g.num_vertices);
				for(int target_vertex=0;target_vertex<g.num_vertices;target_vertex++) {
					if(!true_edges[target_vertex]) {//can be used as fake edge
						fake_edges.add(target_vertex);
					}
				}
				Collections.shuffle(fake_edges);//also ensure uniform picks of fake edges
				if(my_neighbors.size()+fake_edges.size()!=g.num_vertices){//Sanity check
					System.err.println("my_neighbors.size()+fake_edges.size()!=g.num_vertices");
				}
				M_Part_Partition[] my_partitions = new M_Part_Partition[sanitized_number_of_outgoing_egdes];
				for(int i=0;i<my_partitions.length;i++) {
					if(i<my_neighbors.size()) {//ordinary partition
						my_partitions[i] = new M_Part_Partition(i, my_neighbors);
					}else{//fill-up partition
						my_partitions[i] = new M_Part_Partition(i);
					}
				}
				//Randomly distribute fake edges
				for(int i=0;i<fake_edges.size();i++) {
					int f_edge = fake_edges.get(i);
					my_partitions[i%sanitized_number_of_outgoing_egdes].fake_edges.add(f_edge);
				}
					
				//TODO remove me
				M_Part_Partition.rebalance_partitions(my_partitions, epsilon_q2);
				M_Part_Partition.check_partitions(node, my_partitions, my_neighbors, fake_edges, g, neighbors);
				final double probability = M_Part_Partition.get_p(my_partitions, epsilon_q2);
				my_sanitized_neighbors = new ArrayList<Integer>();
				for(int i=0;i<my_partitions.length;i++) {
					int neighbor_s = my_partitions[i].generalized_rr(probability, rand);
					if(neighbor_s<0||neighbor_s>g.num_vertices-1) {
						System.err.println("neighbor_s<0||neighbor_s>g.num_vertices-1");
					}
					my_sanitized_neighbors.add(neighbor_s);
				}
				
				int[] local_error = Metrics.edge_edit_dist(neighbors, my_sanitized_neighbors, g.num_vertices);
				double wrong = local_error[Metrics.missing_edge]+local_error[Metrics.new_fake_edge];
				/*if(Math.abs(error_m_part_2-wrong)>10.0d) {
					System.out.println("id="+node+" Expecting\t"+error_m_part_2+"\tgot\t"+(wrong)+"\terror_rr=\t"+error_rr+"\t"+Arrays.toString(local_error)+"\t|N|="+neighbors.size());	
				}*/
				sum_error_expected += error_m_part_2;
				sum_error_observed += wrong;
			}
			
			sanitized_neighbor_list[node] = my_sanitized_neighbors;
			
			if(node%10000==0) {
				System.out.println("node="+node);
			}
		}
		System.out.println("sum_error_expected="+sum_error_expected+" sum_error_observed="+sum_error_observed+" |G|="+g.num_vertices);
		Graph.stop(start, "m_part()");
		return new Graph(sanitized_neighbor_list, g.name+" "+name);
	}
	
	static Graph m_part_no_seq(Graph g, final double all_epsilon) {
		String name = "m_part_no_seq() e="+all_epsilon;
		System.out.println(name);
		double start = System.currentTimeMillis();
		ArrayList<Integer>[] neighbor_list = g.get_neighbors();
		
		@SuppressWarnings("unchecked")
		ArrayList<Integer>[] sanitized_neighbor_list = new ArrayList[neighbor_list.length]; 
			
		for(int node=0;node<g.num_vertices;node++) {
			ArrayList<Integer> my_sanitized_neighbors = new ArrayList<Integer>();
			ArrayList<Integer> neighbors = neighbor_list[node];		
			
			if(neighbors.isEmpty()) {
				//do nothing
			}else{
				boolean[] true_edges = new boolean[g.num_vertices];
				for(int my_edge : neighbors) {
					true_edges[my_edge] = true;
				}
				ArrayList<Integer> fake_edges = new ArrayList<Integer>(g.num_vertices);
				for(int target_vertex=0;target_vertex<g.num_vertices;target_vertex++) {
					if(!true_edges[target_vertex]) {//can be used as fake edge
						fake_edges.add(target_vertex);
					}
				}
				Collections.shuffle(fake_edges);//also ensure uniform picks of fake edges

				M_Part_Partition[] my_partitions = new M_Part_Partition[neighbors.size()];
				for(int i=0;i<my_partitions.length;i++) {
					my_partitions[i] = new M_Part_Partition(i, neighbors);
				}
				//Randomly distribute fake edges
				for(int i=0;i<fake_edges.size();i++) {
					int f_edge = fake_edges.get(i);
					my_partitions[i%neighbors.size()].fake_edges.add(f_edge);
				}
				
				double s;
				 
				if(neighbors.size()<=1) {
					s = Math.ceil((g.num_vertices-neighbors.size()+1)/(neighbors.size()+1));//special case only larger neighbor exists
				}else{
					s = Math.ceil((g.num_vertices-neighbors.size()+1)/(neighbors.size()-1));//as in paper
				}
				
				final double probability = Math.pow(Math.E, all_epsilon) / (s+Math.pow(Math.E, all_epsilon));
				my_sanitized_neighbors = new ArrayList<Integer>();
				for(int i=0;i<my_partitions.length;i++) {
					int neighbor_s = my_partitions[i].generalized_rr(probability, rand);
					my_sanitized_neighbors.add(neighbor_s);
				}
			}
			
			sanitized_neighbor_list[node] = my_sanitized_neighbors;
			
			if(node%10000==0) {
				System.out.println("node="+node);
			}
		}
		Graph.stop(start, "m_part_no_seq()");
		return new Graph(sanitized_neighbor_list, g.name+" "+name);
	}
	
	
	@Deprecated
	static Graph k_edge_radomized_response_partitioned(Graph g, final double all_epsilon, final double epsilon_q1, boolean rr_fall_back) {
		String name = "k_edge_radomized_response_partitioned() e="+all_epsilon+" e_1="+epsilon_q1+" rr_fall_back="+rr_fall_back;
		System.out.println(name);
		double start = System.currentTimeMillis();
		ArrayList<Integer>[] neighbor_list = g.get_neighbors();
		
		//distribute budget
		final double count_query_epsilon = epsilon_q1;
		final double epsilon_q2 = all_epsilon-epsilon_q1;
		final double p = epsilon_to_p(epsilon_q2);
		
		final boolean[] adjacency_matrix_line = new boolean[g.num_vertices];
		final boolean[] sanitized_adjacency_matrix_line = new boolean[g.num_vertices];
		
		@SuppressWarnings("unchecked")
		ArrayList<Integer>[] sanitized_neighbor_list = new ArrayList[neighbor_list.length]; 
		
		boolean[] duplicate_check = new boolean[g.num_vertices];
		
		for(int node=0;node<g.num_vertices;node++) {
			Arrays.fill(duplicate_check, false);

			//Execute Q1
			int sanitized_number_of_outgoing_egdes = q_1(neighbor_list[node], count_query_senstivity, count_query_epsilon);
			final int edges_per_group = (int) Math.floor((double)g.num_vertices / sanitized_number_of_outgoing_egdes);
			
			//Decide whether we use k-edge or fall back to randomized response
			double error_rr = TheoreticalAnalysis.error_rr(sanitized_number_of_outgoing_egdes, epsilon_q2, g.num_vertices)[0];
			double error_k_edge = TheoreticalAnalysis.error_k_edge_random_part(sanitized_number_of_outgoing_egdes, all_epsilon)[0];
			ArrayList<Integer> my_sanitized_neighbors;
			
			if(error_rr<error_k_edge){
				//System.out.println("Fall back");
				my_sanitized_neighbors = new ArrayList<Integer>();
				Arrays.fill(adjacency_matrix_line, false);
				Arrays.fill(sanitized_adjacency_matrix_line, false);
				for(int target_node : neighbor_list[node]) {
					adjacency_matrix_line[target_node] = true;
				}
				for(int target_node=0;target_node<g.num_vertices;target_node++) {
					double dice = rand.nextDouble();
					boolean create_edge;
					if(dice <= (1-p)) {//with probability 1-p return the real value
						create_edge = adjacency_matrix_line[target_node];
					}else {//else flip the bit
						create_edge = (adjacency_matrix_line[target_node]) ? false : true;
					}
					if(create_edge) {
						my_sanitized_neighbors.add(target_node);
					}
				}
			}else{
				//Partition and Shuffle the data
				int[] my_neighbors = new int[sanitized_number_of_outgoing_egdes];
				int[] fake_edges = new int[sanitized_number_of_outgoing_egdes];//
				
				for(int group = 0;group<sanitized_number_of_outgoing_egdes;group++) {
					int my_neighbor_in_group = get_group_edge_or_random(group, neighbor_list[node], sanitized_number_of_outgoing_egdes, edges_per_group, rand);
					my_neighbors[group] = my_neighbor_in_group;
					if(duplicate_check[my_neighbor_in_group]) {
						System.err.println("duplicate_check[my_neighbor_in_group] at "+my_neighbor_in_group);
					}
					duplicate_check[my_neighbor_in_group] = true;
					
					while(true){
						//get some random edge form this group
						int i = rand.nextInt(edges_per_group);
						int e=i*sanitized_number_of_outgoing_egdes;
						e+=group;
						if(my_neighbor_in_group!=e) {//return some other edge from this group. It's ok to return an existing edge.
							fake_edges[group] = e;
							if(duplicate_check[e]) {
								System.err.println("duplicate_check[e] at "+e);
							}
							duplicate_check[e] = true;
							break;
						}
					}
				}
				Partition partition = new Partition(my_neighbors, fake_edges);
				
				//Execute Q2
				my_sanitized_neighbors = q_2(partition.my_neighbors, partition.fake_edges, sanitized_number_of_outgoing_egdes, epsilon_q2, rand);	
			
			}
			
			sanitized_neighbor_list[node] = my_sanitized_neighbors;
			
			if(node%10000==0) {
				System.out.println("node="+node);
			}
		}
		Graph.stop(start, "k_edge_radomized_response_partitioned()");
		return new Graph(sanitized_neighbor_list, g.name+" "+name);
	}
	
	static int get_group_edge_or_random(int group, ArrayList<Integer> neighbors, int num_groups, int edges_per_group, Random rand_node) {
		for(int e : neighbors) {
			if(e % num_groups == group) {
				return e;
			}
		}
		//There is no edge in this group: Invent some.
		int e = rand_node.nextInt(edges_per_group);
		e*=num_groups;
		e+=group;
		return e;
	}

	public static final double epsilon_to_p(final double epsilon) {
		if(epsilon<0) {
			System.err.println("epsilon<0");
			return 0.5d;//coin flip
		}
		double p = 1.0d/(1.0d+Math.pow(Math.E, epsilon));
		return p;
	}
	
	static final int q_1(final ArrayList<Integer> my_neighbors, final double senstivity, final double epsilon_1) {
		return sanitize(my_neighbors.size(), senstivity, epsilon_1);
	}
	
	static final ArrayList<Integer> q_2(final ArrayList<Integer> my_neighbors, final ArrayList<Integer> fake_edges,
			final int sanitized_number_of_outgoing_egdes, final double epsilon_q2, final Random rand) {
		//Run the individual mechanisms
		ArrayList<Integer> my_sanitized_neighbors = new ArrayList<Integer>(sanitized_number_of_outgoing_egdes);
		final double p = epsilon_to_p(epsilon_q2);
		
		int position = 0;
		for(;position<Math.min(my_neighbors.size(), sanitized_number_of_outgoing_egdes); position++) {
			double dice = rand.nextDouble();
			if(dice <= (1-p)) {//with probability 1-p return the real value
				my_sanitized_neighbors.add(my_neighbors.get(position));
			}else {//else return some fake neighbor
				my_sanitized_neighbors.add(fake_edges.get(position));
			}
		}
		//In case, we need to invent edges, pick the uniformly form the fake edges.
		for(;position<sanitized_number_of_outgoing_egdes; position++) {
			// uniformly pick some edge (a) not part of the real edges and (b) not picked before
			// This is ensured by shuffling the fake edges
			my_sanitized_neighbors.add(fake_edges.get(position)); 
		}
		Collections.sort(my_sanitized_neighbors);
		return my_sanitized_neighbors;
	}
	
	static final ArrayList<Integer> q_2_seq_composition(final ArrayList<Integer> my_neighbors, final ArrayList<Integer> fake_edges,
			final int sanitized_number_of_outgoing_egdes, final double epsilon_q2, final Random rand) {
		//Run the individual mechanisms
		ArrayList<Integer> my_sanitized_neighbors = new ArrayList<Integer>(sanitized_number_of_outgoing_egdes);
		final double p = epsilon_to_p(epsilon_q2/(double)sanitized_number_of_outgoing_egdes);//Here, we see the seq composition
		
		int position = 0;
		for(;position<Math.min(my_neighbors.size(), sanitized_number_of_outgoing_egdes); position++) {
			double dice = rand.nextDouble();
			if(dice <= (1-p)) {//with probability 1-p return the real value
				my_sanitized_neighbors.add(my_neighbors.get(position));
			}else {//else return some fake neighbor
				my_sanitized_neighbors.add(fake_edges.get(position));
			}
		}
		//In case, we need to invent edges, pick the uniformly form the fake edges.
		for(;position<sanitized_number_of_outgoing_egdes; position++) {
			// uniformly pick some edge (a) not part of the real edges and (b) not picked before
			// This is ensured by shuffling the fake edges
			my_sanitized_neighbors.add(fake_edges.get(position)); 
		}
		Collections.sort(my_sanitized_neighbors);
		return my_sanitized_neighbors;
	}
	
	/**
	 * directed
	 * @param g
	 */
	static int[][] two_k_series(Graph g){
		ArrayList<Integer>[] neighbor_list = g.get_neighbors();
		int[] out_degree = new int[g.num_vertices];
		int[] in_degree = new int[g.num_vertices];
		
		//(1) Compute the degree attributes per node
		HashMap<Integer, ArrayList<Integer>> in_degree_temp = new HashMap<Integer, ArrayList<Integer>>(g.num_vertices);
		for(int id=0;id<neighbor_list.length;id++) { in_degree_temp.put(id, new ArrayList<Integer>()); }
		for(int id=0;id<neighbor_list.length;id++) {
			out_degree[id] = neighbor_list[id].size();
			for(int target : neighbor_list[id]) {
				in_degree_temp.get(target).add(id);//TODO maybe only count??
			}
		}
		//make in degree an array
		for(int id=0;id<neighbor_list.length;id++) {in_degree[id] = in_degree_temp.get(id).size();}
		
		//(2) Grouping - get the keys
		HashSet<ArrayList<Integer>> grouping_keys = new HashSet<ArrayList<Integer>>(neighbor_list.length);//need to use an ArrayList, because int[] does not overwrite equals
		for(int id=0;id<neighbor_list.length;id++) {
			ArrayList<Integer> g_key =  new ArrayList<Integer>();
			g_key.add(out_degree[id]);
			g_key.add(in_degree[id]);
			grouping_keys.add(g_key);
		}
		
		//(2) Grouping - copy grouping keys into an array s.t. I can address them by some offset
		final int num_groups = grouping_keys.size();
		final ArrayList<int[]> all_groups_dense = new  ArrayList<int[]>(num_groups); 
		for(ArrayList<Integer> k : grouping_keys) {
			int[] key = {k.get(0), k.get(1)};
			all_groups_dense.add(key);
		}
		Collections.sort(all_groups_dense, new Comparator<int[]>() {
		    public int compare(int[] a, int[] b) {
		        if(a[0]!=b[0]) {
		        	return Integer.compare(a[0],b[0]);
		        }else {
		        	return Integer.compare(a[1],b[1]);
		        }
		    }
		});
		
		//(2.1) Grouping - compute the grouping - (a) map nodes to groups and (b) groups to list of nodes
		ArrayList<HashSet<Integer>> nodes_in_group 	= new  ArrayList<HashSet<Integer>>(num_groups);
		HashMap<Integer, Integer> node_to_group 	= new HashMap<Integer, Integer>(g.num_vertices);
		for(int group=0;group<num_groups;group++) {
			int[] group_key = all_groups_dense.get(group);
			HashSet<Integer> my_nodes = new HashSet<Integer>();
			nodes_in_group.add(my_nodes);
			//loop over all nodes
			for(int node=0;node<g.num_vertices;node++) {
				int[] node_key = {out_degree[node], in_degree[node]};
				if(group_key[0] == node_key[0] && group_key[1] == node_key[1]) {//arrays in java ...
					my_nodes.add(node);
					node_to_group.put(node, group);
				}
			}
			if(my_nodes.isEmpty()) {
				System.err.println("two_k_series(Graph): my_nodes.isEmpty()");
			}
		}
		if(node_to_group.size()!=g.num_vertices) {
			System.err.println("two_k_series(Graph): my_nodes.isEmpty()");
		}
		//(2.2) Grouping - compute aggregation
		int[][] result = new int[num_groups][num_groups];
		for(int group=0;group<num_groups;group++) {
			//loop over all nodes 
			for(int node : nodes_in_group.get(group)) {
				//loop over all edges per node to find their target group
				for(int edge_target : neighbor_list[node]) {
					//check whether is an edge to the other_group
					Integer edge_target_group_number = node_to_group.get(edge_target);
					if(edge_target_group_number==null) {
						System.err.println("edge_target_group_number==null");
					}
					result[group][edge_target_group_number]++;
				}
			}
			
		}
		
		return result;
	}
	
	static final int sanitize(final double val, final double sensitivity, final double epsilon) {
		final double lambda = sensitivity/epsilon;
		return sanitize(val, lambda);
	}
	
    /**
     * Special method for 1-D Query results.
     * Sanitizes the query release of one time stamp according to the provided noise scale lambda
     * 
     * @param val original 1-dimensional query result at time t
     * @param lambda noise scale to add to each query dimension
     * 
     * @return non-negative integer
     */
    public static final int sanitize(final double val, final double lambda) {
        double noise = LaplaceStream.nextNumber() * lambda;
        double sanVal = val + noise;        
        //System.out.println(sanVal);
        
        if (TRUNCATE) {
        	//Round query result to a non-negative number. Especially, if the stream is sparse (i.e., contains a lot of zero counts), this highly improves utility.
            return (int) Math.max(0, Math.round(sanVal));
        }
        return (int)sanVal;
    }
    
    private static final double sanitize_double(final double val, final double lambda) {
        double noise = LaplaceStream.nextNumber() * lambda;
        double sanVal = val + noise;        
        //System.out.println(sanVal);
        
        return sanVal;
    }


	/**
	 * top_k
	 * @param g
	 * @param p
	 * @return
	 */
	static Graph top_k(Graph g, final double all_epsilon, final double epsilon_q1, boolean use_seq_composition) {
		String name = "top_k() use_seq_composition="+use_seq_composition;
		System.out.println(name);
		double start = System.currentTimeMillis();
		ArrayList<Integer>[] neighbor_list = g.get_neighbors();
		
		//distribute budget
		final double count_query_epsilon = epsilon_q1;
		final double epsilon_q2 = all_epsilon-epsilon_q1;
		
		@SuppressWarnings("unchecked")
		ArrayList<Integer>[] sanitized_neighbor_list = new ArrayList[neighbor_list.length]; 
		
		/**
		 * Buffer for the scoring function
		 */
		double[] scoring_function_u = new double[g.num_vertices];
		
		for(int node=0;node<g.num_vertices;node++) {
			final ArrayList<Integer> my_true_neighbors = neighbor_list[node];
			//Execute Q1
			int sanitized_number_of_outgoing_egdes = q_1(my_true_neighbors, count_query_senstivity, count_query_epsilon);
			
			//Compute the (noisy) scoring function
			Arrays.fill(scoring_function_u, 0);
			compute_scoring_function(scoring_function_u, my_true_neighbors);
			add_noise_to_scoring_function(scoring_function_u, epsilon_q2, sanitized_number_of_outgoing_egdes, use_seq_composition);
			
			TopK col = TopK.create_With_First_K_Elements(sanitized_number_of_outgoing_egdes, scoring_function_u);
			for(int node_id = sanitized_number_of_outgoing_egdes; node_id < g.num_vertices; node_id++) {
				double score = scoring_function_u[node_id];
				col.tryUpdate(node_id, score);
			}
			
			ArrayList<Integer> my_sanitized_neighbors = new ArrayList<Integer>(sanitized_number_of_outgoing_egdes);
			sanitized_neighbor_list[node] = my_sanitized_neighbors;
			for(int target : col.node_ids) {
				my_sanitized_neighbors.add(target);
			}
			Collections.sort(my_sanitized_neighbors);

			if(node%10000==0) {
				System.out.println("node="+node);
			}
		}
		Graph.stop(start, "top_k()");
		return new Graph(sanitized_neighbor_list,g.name+" top_k() use_seq_composition="+use_seq_composition);
	}

	private static void add_noise_to_scoring_function(final double[] scoring_function_u, final double epsilon, int k, final boolean use_seq_composition) {
		add_noise_to_scoring_function(scoring_function_u, epsilon, k, use_seq_composition, 1.0);//Sensitivity is 1
	}
	
	private static void add_noise_to_scoring_function(final double[] scoring_function_u, final double epsilon, int k, final boolean use_seq_composition, final double sensitivity_u) {
		double lambda; 
		if(use_seq_composition) {
			lambda = sensitivity_u*(double)k / epsilon;//distribute epsilon among k queries, i.e., respect seq composition.
		}else{
			lambda = sensitivity_u / epsilon;
		}
		for(int node = 0;node<scoring_function_u.length;node++) {
			double org_score = scoring_function_u[node];
			double noisy_u = sanitize_double(org_score, lambda);
			scoring_function_u[node] = noisy_u;
		}
	}

	private static void compute_scoring_function(double[] scoring_function_u, ArrayList<Integer> my_true_neighbors) {
		for(int neighbor : my_true_neighbors) {
			scoring_function_u[neighbor] = 1;
		}
	}
	
	/**
	 * Implements the private graph release approach based on 2K series statistics. The approach involves two main problems:
	 * (1) Reduction of the noise injection, due to the high sensitivity of the 2kseries. Here, we rely on the approach in "Sharing Graphs Using Differentialy Private Graph Models"
	 * (2) How to reconstruct the graph from the 2kseries. Here we rely on "2K+ Graph Construction Framework: Targeting Joint Degree Matrix and Beyond", https://escholarship.org/uc/item/5tb4f2mn
	 * @param g
	 * @param all_epsilon
	 * @param not_private
	 * @return
	 */
	public static Graph two_k_series(Graph g, final double all_epsilon, final boolean not_private) {
		String name = "two_k_series()";
		System.out.println(name);

		double start = System.currentTimeMillis();
		DK_Series dks = algorithms.DK_Series.two_k_series(g);
		if(!not_private) {
			dks.sanitize(all_epsilon);
		}
		Graph g_san = dks.reconstruct();
		
		Graph.stop(start, name);
		return g_san;
	}
	
	/**
	 * top_k
	 * @param g
	 * @param p
	 * @return
	 */
	static Graph educated_guess(Graph g, final double all_epsilon, final boolean non_private) {
		String name = "educated_guess()";
		System.out.println(name);
		double start = System.currentTimeMillis();
		ArrayList<Integer>[] neighbor_list = g.get_neighbors();
		
		//distribute budget
		final double count_query_epsilon = all_epsilon/2;
		final double epsilon_q2 = all_epsilon-count_query_epsilon;
		
		@SuppressWarnings("unchecked")
		ArrayList<Integer>[] sanitized_neighbor_list = new ArrayList[neighbor_list.length]; 
		
		final int[] out_degress = new int[g.num_vertices];
		/**
		 * Frequency of each node as target
		 */
		final int[] histogram   = new int[g.num_vertices];
		
		for(int node=0;node<g.num_vertices;node++) {//TODO locally DP?
			final ArrayList<Integer> my_true_neighbors = neighbor_list[node];
			for(int neighbor : my_true_neighbors) {
				histogram[neighbor]++;
			}
			//Execute Q1 privately
			out_degress[node] = (non_private) ? my_true_neighbors.size() : q_1(my_true_neighbors, count_query_senstivity, count_query_epsilon);
		}
		
		//make histogram private
		if(!non_private) {
			double lambda = 1.0d / epsilon_q2;
			for(int node=0;node<g.num_vertices;node++) {
				int count = histogram[node];
				count = sanitize(count, lambda);
				histogram[node] = count;
			}
		}
		final double[] probabilites = to_probabilities(histogram);
		
		//dice the edges
		for(int node=0;node<g.num_vertices;node++) {
			int out_degree = out_degress[node];
			
			ArrayList<Integer> my_sanitized_neighbors = new ArrayList<Integer>(out_degree);
			sanitized_neighbor_list[node] = my_sanitized_neighbors;
			int counter = 0;
			HashSet<Integer> edges = new HashSet<Integer>();
			while(counter<g.num_vertices && edges.size()<out_degree) {
				int edge = random_choice(probabilites, rand);
				edges.add(edge);
				if(counter == g.num_vertices-2) {
					System.err.println("Warning educated_guess() counter == g.num_vertices-2");
				}
			}
			for(int egde : edges) {
				my_sanitized_neighbors.add(egde);
			}
			Collections.sort(my_sanitized_neighbors);
		}
		Graph.stop(start, "educated_guess()");
		return new Graph(sanitized_neighbor_list,g.name+" educated_guess()");
	}
	
	static double[] to_probabilities(int[] noissy_u) {
		double sum = 0.0d;
		for(double d: noissy_u) {
			sum += (double)d;
		}
		double[] probabilites = new double[noissy_u.length];
		for(int i=0;i<probabilites.length;i++) {
			probabilites[i] = ((double)noissy_u[i]) / sum;
		}
		return probabilites;
	}
	
	static void to_probabilities(double[] array) {
		double sum = 0.0d;
		for(double d: array) {
			sum += (double)d;
		}
		for(int i=0;i<array.length;i++) {
			array[i] /= sum;
		}
	}
	
	/**
	 * Generating Synthetic Decentralized Social Graphs with Local Differential Privacy, CCS'17
	 * 
	 * @param g - original graph
	 * @param epsilon
	 * @param epsilon_1
	 * @return
	 * 
	 * XXX the problem is it gets worse with increasing, because we are omitting parts of the weight
	 */
	public static Graph LDPGen(final Graph g, final double epsilon, final boolean none_private) {
		String name = "LDPGen";
		if(!none_private) {
			System.out.println(name +" e="+ epsilon);	
		}else{
			System.err.println(name +" e="+ epsilon+ " none private");
		}
		double start = System.currentTimeMillis();
		
		//XXX According to Sec 4.7
		final double epsilon_1 = epsilon/2.0d; 
		final int k_0 = 2;
		
		//Phase 1 - Randomly assign vertices to k partitions
		/**
		 * Xi_0 vertex -> partition
		 */
		HashMap<Integer, Integer> partition_assigenments = new HashMap<Integer, Integer>(g.num_vertices);
		for(int vertex = 0; vertex<g.num_vertices;vertex++) {
			partition_assigenments.put(vertex, vertex%k_0);
		}
		
		//Phase 1 - Compute degree vector
		ArrayList<Integer>[] N = g.get_neighbors();
		double lambda = 1.0d / epsilon_1;
		/**
		 * Noisy out degree estimation for each vertex
		 */
		final int[] all_eta = new int[g.num_vertices];
		/**
		 * Noisy out degree estimation for each vertex to each partition
		 * usage all_delta[vertex][partition]
		 */
		final int[][] all_delta = new int[g.num_vertices][];
		for(int vertex = 0; vertex<g.num_vertices;vertex++) {
			ArrayList<Integer> my_neighbors = N[vertex];
			/**
			 * Noisy degree vector of vertex: \delta^u in the paper. Basically a histogram.
			 */
			final int[] delta_u = new int[k_0];
			for(int target : my_neighbors) {
				int target_partition = partition_assigenments.get(target);
				delta_u[target_partition]++;
			}
			for(int bin=0;bin<k_0;bin++) {
				if(!none_private) {
					int san_bin_value = sanitize(delta_u[bin], lambda);
					delta_u[bin] = san_bin_value;	
				}else{
					delta_u[bin] = delta_u[bin];//The value remains as it is
				}
				
			}
			final int eta = IntStream.of(delta_u).sum();//XXX Equation 1, page 429. Why does that work?
			all_delta[vertex] = delta_u;
			all_eta[vertex] = eta;
		}
		
		//Phase 1 - estimate optimal number of partitions k_1
		final int max_degree = Util.max(all_eta);
		/**
		 * H[degree] -> Pr(degree)
		 * sum over all degrees = 1.0
		 */
		double[] H = new double[max_degree+1];
		for(int degree : all_eta) {
			H[degree]++;
		}
		for(int bin=0;bin<=max_degree;bin++) {//make H the probability histogram
			H[bin] /= (double)g.num_vertices;
		}
		double sum = 0.0d;
		int k_1 = k_0;//XXX Equation 20 wants to estimate k_1 but contains k_1 on the left hand side, guess its an error using k_0 instead.
		for(int bin=0;bin<=max_degree;bin++) {//make H the probability histogram
			sum += H[bin] * (double)k_1*0.5*(double)bin;
		}
		k_1 = (int) Math.ceil(Math.min(sum, 50));//XXX Papers assumes not more than 50 partitions for Eq. 15
		
		Clustering xi_1 = new KMedoids(all_delta, k_1).run();
		HashMap<Integer, Integer> xi_1_hashed = xi_1.to_hash_map();
		
		Graph.stop(start, "Phase 1 k_1="+k_1);
		//System.out.println(xi_1);
		
		//Phase 2 - Compute noisy degree vector to each partition according to xi_1
		double epsilon_2 = epsilon - epsilon_1;
		lambda = 1.0d / epsilon_2;
		/**
		 * \delta^1_u
		 */
		final int[][] all_degree_vector_xi_1 = new int[g.num_vertices][];
		for(int vertex = 0; vertex<g.num_vertices;vertex++) {
			ArrayList<Integer> my_neighbors = N[vertex];
			/**
			 * \delta^u in the paper. Basically a histogram.
			 */
			final int[] my_degree_vector = new int[k_1];
			for(int target : my_neighbors) {
				int target_partition = xi_1_hashed.get(target);//This is the difference to Phase 1
				my_degree_vector[target_partition]++;
			}
			for(int bin=0;bin<k_1;bin++) {
				if(!none_private) {
					int san_bin_value = sanitize(my_degree_vector[bin], lambda);
					my_degree_vector[bin] = san_bin_value;	
				}// else the value remains as it is
			}
			all_degree_vector_xi_1[vertex] = my_degree_vector;
			final int noisy_degree_estimator_eta = IntStream.of(my_degree_vector).sum();//Equation 1, page 429
			all_eta[vertex] = noisy_degree_estimator_eta;
		}
		//Phase 2 - Sec 4.4 Cluster again using the newly computed degree estimation
		Clustering xi_2 = new KMedoids(all_degree_vector_xi_1, k_1).run();
		
		//Phase 2 - Sec 4.4 estimate degree vectors for xi_2 based on noisy degree vector of xi_1
		double[][] final_degree_estimations = new double[g.num_vertices][];
		for(int vertex = 0; vertex<g.num_vertices;vertex++) {
			final int[] delta_xi_1 = all_degree_vector_xi_1[vertex]; 
			double[] my_degree_estimation = new double[k_1];
			for(int i = 0; i < k_1; i++) {//Over all clusters
				ArrayList<Integer> U_2_i = xi_2.getCluster(i);//U^2_i in the paper
				double estimation = estimate(vertex, U_2_i, xi_1.ALL_CLUSTERS, delta_xi_1);
				my_degree_estimation[i] = estimation;
			}
			final_degree_estimations[vertex] = my_degree_estimation;
		}
		
		Graph.stop(start, "Phase 2");
		
		//Phase 3 - generate edges
		Graph san_g = new Graph(g.num_vertices, g.name+"_"+name);
		//ldp_gen_phase_3(san_g, final_degree_estimations, k_1, xi_2);
		create_vertices(san_g, final_degree_estimations, k_1, xi_2);
		
		Graph.stop(start, name);
		return san_g;
	}
	
	
	public static int random_choice(final double[] probabilities, Random rand) {
		final double threshold = rand.nextDouble(); //Some value in [0,1]. We return the
		double prob = 0.0d;
		for(int i=0;i<probabilities.length;i++) {
			double p = probabilities[i];
			prob += p;
			if(prob>=threshold) {
				return i;
			}
		}
		System.err.println("Should never come till here");
		return -1;
	}
	
	static void create_vertices(final Graph san_g, final double[][] final_degree_estimations, final int k_1, final Clustering xi_2){
		for(int i=0;i<k_1;i++) {
			ArrayList<Integer> U_i = xi_2.getCluster(i);
			for(int u : U_i) {
				double[] my_deseired_out_degree = final_degree_estimations[u];
				for(int j = 0; j < k_1; j++) {
					ArrayList<Integer> U_j = xi_2.getCluster(j);
					double[] probabilities = to_probabilities(U_j, i, final_degree_estimations);
					HashSet<Integer> my_neighbors = new HashSet<Integer>();
					for(int neighbor=0;neighbor<Math.round(my_deseired_out_degree[j]);neighbor++) {
						int index = random_choice(probabilities, rand);
						my_neighbors.add(U_j.get(index));//avoids duplicate edges
					}
					for(Integer v : my_neighbors) {
						if(Config.materialize_graph) {
							Graph.write_edge(u, v);
						}else{
							san_g.add_edge(u, v);	
						}
					}
				}
			}
		}
		double sum = Util.sum(final_degree_estimations);
		System.out.println("LDP Gen Phase 3 create_vertices(): Expected "+sum+" vertices and got "+san_g.int_num_edges());
		System.out.println("Done");
	}
	
	private static double[] to_probabilities(ArrayList<Integer> vertices, int target_partition, final double[][] final_degree_estimations) {
		double[] p_s = new double[vertices.size()];
		for(int i=0;i<vertices.size();i++) {
			int vertex = vertices.get(i);
			p_s[i] = 0.01+final_degree_estimations[vertex][target_partition];//avoid all zero case with some outdegree
		}
		double sum = Util.sum(p_s);
		for(int i=0;i<vertices.size();i++) {
			p_s[i] /= sum;
		}
		//System.out.println(sum(p_s));
		return p_s;
	}

	/**
	 * An application/interpretation of the Chun-Lu Model 
	 * There, the input is the desire out degree per node as weight-vector W, with |W| = |V|
	 * The probability of creating an edge (u,v) is W[u] * W[v] / sum(W)
	 * 
	 * Note, it must hold for any W[u], W[v] that W[u] * W[v] < sum(W)
	 * </p>
	 * William Aiello, Fan Chung, and Linyuan Lu. 2000. A random graph model for massive graphs. In Proceedings of the thirty-second annual ACM symposium on Theory of computing. Acm, 171–180.
	 * 
	 * @param san_g - Zhe empty sanitized graph to fill with edges
	 * @param final_degree_estimations - Desired out degree of each vertex to each partition, \delta^u_j = final_degree_estimations[u][j]
	 * @param k_1 - original graph number of partitions
	 * @param xi_2 - original graph the graph partitions
	 * 
	 */
	static void ldp_gen_phase_3(final Graph san_g, final double[][] final_degree_estimations, final int k_1, final Clustering xi_2){
		/*out_tsv(final_degree_estimations);
		for(int partition=0;partition<1;partition++) {
			System.out.println("**********Cluster "+partition);
			for(Integer arr : xi_2.getCluster(partition)) {
				out_tsv(final_degree_estimations[arr]);
			}
		}*/

		ArrayList<double[]> all_partition_weight_sums = new ArrayList<double[]>(k_1);
		
		for(int partition=0;partition<k_1;partition++) {
			double[] partition_weight_sum = new double[k_1];	
			all_partition_weight_sums.add(partition_weight_sum);
			
			final ArrayList<Integer> p = xi_2.getCluster(partition);
			for(int vertex : p) {
				add(partition_weight_sum, final_degree_estimations[vertex]);
			}
		}
		
		for(int i=0;i<k_1;i++) {
			ArrayList<Integer> U_i = xi_2.getCluster(i);
			final double[] sum_weights_U_i = all_partition_weight_sums.get(i);
			for(int u : U_i) {
				for(int j = 0; j < k_1; j++) {
					ArrayList<Integer> U_j = xi_2.getCluster(j); 
					double[] sum_weights_U_j = all_partition_weight_sums.get(j);
					for(int v : U_j) {
						if(Config.USE_FIX) {
							/**
							 * Desired out degree of u to some vertex v \in U_j. This is correct in the paper.
							 */
							final double w_u = (i==j) ? final_degree_estimations[u][j] : final_degree_estimations[u][j]+final_degree_estimations[u][i];	//FIX
							//final double w_u = final_degree_estimations[u][j];
							/**
							 * Average desired out degree of all vertex v \in U_j to a vertex \in U_i. 
							 * Technically, we could also simply use final_degree_estimations[v][i]. I guess the authors wanted to mitigate some noise form the randomization by computing the average.
							 * Note, presuming this is an instance of the Chun-Lu Model. The paper is wrong here. It iterates over all v \in G and consider the out degree to U_j. 
							 * However, canonically w_v is the out degree of v to the rest of the graph, which is here only U_i, because of the partitioning. We corrected that, also in the denominator.
							 * This is backuped by the statement "the denominator is the total number of edges between groups i and j by aggregating all the elements in the nodes’ corresponding degree vectors in two gr"
							 */
							final double w_v = (i==j) ? final_degree_estimations[v][i] :  final_degree_estimations[v][i]+ final_degree_estimations[v][j];//FIX 
							//final double w_v = sum_weights_U_j[i] / U_j.size();
							final double sum_weights = sum_weights_U_i[j] + sum_weights_U_j[i];
							/**
							 * probability of drawing a directed edge (u,v)
							 */
							double p = w_u * w_v / sum_weights;
							if(p>1.0) {
								//System.err.println("p>1.0 p="+p);
								p=1.0d;
							}
							double dice = rand.nextDouble();
							if(dice<=p) {
								san_g.add_edge(u, v);
							}
						}else {
							/**
							 * Desired out degree of u to some vertex v \in U_j. This is correct in the paper.
							 */
							//final double w_u = (i==j) ? final_degree_estimations[u][j] : final_degree_estimations[u][j]+final_degree_estimations[u][i];	//FIX
							final double w_u = final_degree_estimations[u][j];
							/**
							 * Average desired out degree of all vertex v \in U_j to a vertex \in U_i. 
							 * Technically, we could also simply use final_degree_estimations[v][i]. I guess the authors wanted to mitigate some noise form the randomization by computing the average.
							 * Note, presuming this is an instance of the Chun-Lu Model. The paper is wrong here. It iterates over all v \in G and consider the out degree to U_j. 
							 * However, canonically w_v is the out degree of v to the rest of the graph, which is here only U_i, because of the partitioning. We corrected that, also in the denominator.
							 * This is backuped by the statement "the denominator is the total number of edges between groups i and j by aggregating all the elements in the nodes’ corresponding degree vectors in two gr"
							 */
							//final double w_v = (i==j) ? sum_weights_U_j[i] / U_j.size() : sum_weights_U_i[j] / U_i.size()+sum_weights_U_j[i] / U_j.size();//FIX 
							final double w_v = sum_weights_U_j[i] / U_j.size();
							final double sum_weights = sum_weights_U_i[j] + sum_weights_U_j[i];
							/**
							 * probability of drawing a directed edge (u,v)
							 */
							double p = w_u * w_v / sum_weights;
							if(p>1.0) {
								//System.err.println("p>1.0 p="+p);
								p=1.0d;
							}
							double dice = rand.nextDouble();
							if(dice<=p) {
								san_g.add_edge(u, v);
						}
						}
					}
				}	
			}
		}
		
		double sum = Util.sum(final_degree_estimations);
		System.out.println("LDP Gen Phase 3: Expected "+sum+" vertices and got "+san_g.int_num_edges());
		System.out.println("Done");
	}

	private static void add(final double[] sum, final double[] add_me) {
		if(sum.length!=add_me.length) {
			System.err.println("sum.length!=add_me.length");
		}
		
		for(int i=0;i<sum.length;i++) {
			sum[i] +=  add_me[i];
		}
	}
	
	/**
	 * Compute the estimated degree vector of some vertex according to Sec 4.4
	 * @param vertex - u in the paper
	 * @param u_2_i - U^2_i \in \xi_2
	 * @param u_1_all - grouping of \xi_1 in k_1 clusters
	 * @param delta \delta^u from xi_1 - my old degree estimation from \xi_1 
	 * @return
	 */
	private static double estimate(int vertex, ArrayList<Integer> u_2_i, ArrayList<Integer>[] u_1_all, final int[] delta) {
		double sum = 0.0d;
		for(int j = 0;j<u_1_all.length;j++) {//cluster = j in the paper
			ArrayList<Integer> u_1_j = u_1_all[j];//U^1_j in the paper
			if(u_1_j.isEmpty()) {
				//Do nothing
			}else{
				HashSet<Integer> intersection = new HashSet<Integer>();
				intersection.addAll(u_2_i);
				int count = 0;
				for(int e : u_1_j) {
					if(intersection.contains(e)) {
						count++;
					}
				}
				double temp = count;
				temp /= (double) u_1_j.size();
				temp *= (double) delta[j];//should be \delta^u_j, but is \delta_j in the paper
				if (Double.isNaN(temp)) {
				    System.err.println("NaN");
				}
				sum+= temp;
			}
		}
		return sum;
	}
	
	public static Graph Chung_Lu_Model(Graph g, final double epsilon, final boolean none_private, final boolean partioned) {
		String name = "Chung_Lu e="+epsilon;
		System.out.println(name);
		double start = System.currentTimeMillis();
		
		ArrayList<Integer>[] N = g.get_neighbors();
		Graph san_g = new Graph(g.num_vertices, g.name+"_"+name);
		final double[] weights = new double[g.num_vertices];
		
		//Execute Q 1: get desired out degree
		for(int vertex=0;vertex<g.num_vertices;vertex++) {
			if(none_private) {
				weights[vertex] = N[vertex].size();	
			}else{
				weights[vertex] = q_1(N[vertex], count_query_senstivity, epsilon);
			}
		}
		
		if(partioned) {
			int num_clusters = 2;//guess
			HashMap<Integer, Integer> partition_assigenments = new HashMap<Integer, Integer>(g.num_vertices);
			for(int vertex = 0; vertex<g.num_vertices;vertex++) {
				partition_assigenments.put(vertex, vertex%num_clusters);
			}
			
			final int[][] out_degrees_per_partition = new int[g.num_vertices][];
			for(int vertex = 0; vertex<g.num_vertices;vertex++) {
				ArrayList<Integer> my_neighbors = N[vertex];
				final int[] out_degree = new int[num_clusters];
				for(int target : my_neighbors) {
					int target_partition = partition_assigenments.get(target);
					out_degree[target_partition]++;
				}
				out_degrees_per_partition[vertex] = out_degree;
			}
			Clustering clusters = new KMedoids(out_degrees_per_partition, num_clusters).run();
			System.out.println(clusters);
			int k_2 = 3;
			for(int vertex = 0; vertex<g.num_vertices;vertex++) {
				ArrayList<Integer> my_neighbors = N[vertex];
				final int[] out_degree = new int[k_2];
				for(int target : my_neighbors) {
					int target_partition = partition_assigenments.get(target);
					out_degree[target_partition]++;
				}
				out_degrees_per_partition[vertex] = out_degree;
			}
			clusters = new KMedoids(out_degrees_per_partition, k_2).run();
			partition_assigenments = clusters.to_hash_map();
			System.out.println(clusters);
			
			final double[][] temp = new double[g.num_vertices][k_2];
			for(int vertex = 0; vertex<g.num_vertices;vertex++) {
				ArrayList<Integer> my_neighbors = N[vertex];
				final int[] out_degree = new int[k_2];
				for(int target : my_neighbors) {
					int target_partition = partition_assigenments.get(target);
					temp[vertex][target_partition]++;
				}
				out_degrees_per_partition[vertex] = out_degree;
			}
			
			for(int vertex = 0; vertex<g.num_vertices;vertex++) {
				ArrayList<Integer> my_neighbors = N[vertex];
				final int[] out_degree = new int[k_2];
				for(int target : my_neighbors) {
					int target_partition = partition_assigenments.get(target);
					out_degree[target_partition]++;
				}
				out_degrees_per_partition[vertex] = out_degree;
			}
			
			ldp_gen_phase_3(san_g, temp, k_2, clusters);
		}else{
			//Q2 here costs no budget
			final double sum_weights = Util.sum(weights);
			for(int u=0;u<g.num_vertices;u++) {
				for(int v=0;v<g.num_vertices;v++) {
					double p = weights[u] * weights[v] / sum_weights;
					if(p>1.0) {
						//System.err.println("p>1.0 p="+p);
						p=1.0d;
					}
					double dice = rand.nextDouble();
					if(dice<p) {
						if(Config.materialize_graph) {
							Graph.write_edge(u, v);
						}else{
							san_g.add_edge(u, v);	
						}
					}
				}
				if(u%1000==0) {
					System.out.print("node="+u+" ");
				}
			}
			System.out.println();
		}
		
		Graph.stop(start, name);
		return san_g;
	}
}
