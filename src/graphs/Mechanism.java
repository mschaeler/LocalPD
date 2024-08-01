package graphs;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Random;

import algorithms.LaplaceStream;
import misc.TopK;
import results.TheoreticalAnalysis;

public class Mechanism {
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
		System.out.println("radomized_response()");
		double start = System.currentTimeMillis();
		Random rand = new Random(123456);
		final double p = epsilon_to_p(epsilon);
		ArrayList<Integer>[] neighbor_list = g.get_neighbors();
		Graph private_g = new Graph(g.num_vertices, g.name+"_rr");
		
		final boolean[] adjacency_matrix_line = new boolean[g.num_vertices];
		final boolean[] sanitized_adjacency_matrix_line = new boolean[g.num_vertices];
		
		for(int node=0;node<g.num_vertices;node++) {
			Arrays.fill(adjacency_matrix_line, false);
			Arrays.fill(sanitized_adjacency_matrix_line, false);
			for(int target_node : neighbor_list[node]) {
				adjacency_matrix_line[target_node] = true;
			}
			for(int target_node=0;target_node<g.num_vertices;target_node++) {
				double dice = rand.nextDouble();
				if(dice <= (1-p)) {//with probability 1-p return the real value
					sanitized_adjacency_matrix_line[target_node] = adjacency_matrix_line[target_node];
				}else {//else flip the bit
					sanitized_adjacency_matrix_line[target_node] = (adjacency_matrix_line[target_node]) ? false : true;
				}
				if(sanitized_adjacency_matrix_line[target_node]) {
					private_g.add_edge(node, target_node);
				}
			}
			if(node%10000==0) {
				System.out.println("node="+node);
			}
		}
		Graph.stop(start, "radomized_response()");
		
		return private_g;
	}
	
	/**
	 * Implements no shuffle etc. Sufficient for benchmarking utility metrics but not attack resistance. 
	 * @param my_neighbors
	 * @param sanitized_number_of_outgoing_egdes
	 * @return
	 */
	@Deprecated
	static Partition partitioning_fast(final int node_seed, final ArrayList<Integer> my_neighbors, final int sanitized_number_of_outgoing_egdes, final int num_vertices) {
		HashSet<Integer> my_neighbors_sanitized = new HashSet<Integer>();
		my_neighbors_sanitized.addAll(my_neighbors);
		final HashSet<Integer> fake_edges = new HashSet<Integer>(sanitized_number_of_outgoing_egdes);
		Random rand = new Random(node_seed);//This is the trick: the node id is the seed
		int inifinite_loop_counter = 0;
		
		while(fake_edges.size()!=sanitized_number_of_outgoing_egdes) {
			int fake_edge_id = rand.nextInt(num_vertices);
			if(!my_neighbors_sanitized.contains(fake_edge_id)){
				fake_edges.add(fake_edge_id);
			}
			inifinite_loop_counter++;
			if(inifinite_loop_counter>num_vertices) {
				System.err.println("Infinite loop?");
			}
		}
		
		rand = new Random(node_seed+1);//Take a different seed
		while(my_neighbors_sanitized.size()<sanitized_number_of_outgoing_egdes) {
			int fake_edge_id = rand.nextInt(num_vertices);
			if(!my_neighbors_sanitized.contains(fake_edge_id)){
				my_neighbors_sanitized.add(fake_edge_id);//add a fake edge and pretend it is a real one
			}
			inifinite_loop_counter++;
			if(inifinite_loop_counter>num_vertices) {
				System.err.println("Infinite loop?");
			}
		}
		
		return new Partition(my_neighbors_sanitized, fake_edges);
	}
	
	/**
	 * Version with non-private grouping
	 * @param g
	 * @param p
	 * @return
	 */
	@Deprecated
	public static Graph baseline(Graph g, final double all_epsilon, final double epsilon_q1, boolean use_seq_composition) {
		String name = "baseline() use_seq_composition="+use_seq_composition;
		System.out.println(name);
		double start = System.currentTimeMillis();
		Random rand = new Random(123456);
		ArrayList<Integer>[] neighbor_list = g.get_neighbors();
		
		//distribute budget
		final double count_query_epsilon = epsilon_q1;
		final double epsilon_q2 = all_epsilon-epsilon_q1;
		
		@SuppressWarnings("unchecked")
		ArrayList<Integer>[] sanitized_neighbor_list = new ArrayList[neighbor_list.length]; 
		
		for(int node=0;node<g.num_vertices;node++) {
			//Execute Q1
			int sanitized_number_of_outgoing_egdes = q_1(neighbor_list[node], count_query_senstivity, count_query_epsilon);
			
			//Partition the data
			Partition partition = partitioning_fast(node, neighbor_list[node], sanitized_number_of_outgoing_egdes, g.num_vertices);
			if(partition.my_neighbors.size()<sanitized_number_of_outgoing_egdes){
				System.err.println("my_neighbors.size()<sanitized_number_of_outgoing_egdes");
			}
			if(partition.fake_edges.size()!=sanitized_number_of_outgoing_egdes){
				System.err.println("fake_edges.size()!=sanitized_number_of_outgoing_egdes");
			}
			
			//Execute Q2
			ArrayList<Integer> my_sanitized_neighbors = new ArrayList<Integer>(sanitized_number_of_outgoing_egdes);
			double p;
			if(use_seq_composition) {
				p = epsilon_to_p(epsilon_q2/(double)sanitized_number_of_outgoing_egdes);	
			}else{
				p = epsilon_to_p(epsilon_q2);
			}
			
			for(int position=0; position<sanitized_number_of_outgoing_egdes; position++) {
				double dice = rand.nextDouble();
				if(dice <= (1-p)) {//with probability 1-p return the real value
					my_sanitized_neighbors.add(partition.my_neighbors.get(position));
				}else {//else return some fake neighbor
					my_sanitized_neighbors.add(partition.fake_edges.get(position));
				}
			}
			
			sanitized_neighbor_list[node] = my_sanitized_neighbors;
			
			if(node%10000==0) {
				System.out.println("node="+node);
			}
		}
		Graph.stop(start, "baseline()");
		return new Graph(sanitized_neighbor_list, name);
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
	public static Graph k_edge_radomized_response_non_private_grouping(Graph g, final double all_epsilon, final double epsilon_q1, boolean use_seq_composition) {
		String name = "k_edge_radomized_response_non_private_grouping() use_seq_composition="+use_seq_composition;
		System.out.println(name);
		double start = System.currentTimeMillis();
		Random rand = new Random(123456);
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
	
	public static Graph k_edge_radomized_response_non_private_grouping_rr_fallback(Graph g, final double all_epsilon, final double epsilon_q1, boolean use_seq_composition) {
		String name = "k_edge_radomized_response_non_private_grouping_rr_fallback() use_seq_composition="+use_seq_composition;
		System.out.println(name);
		double start = System.currentTimeMillis();
		Random rand = new Random(123456);
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

	public static Graph k_edge_radomized_response_partitioned(Graph g, final double all_epsilon, final double epsilon_q1, boolean rr_fall_back) {
		String name = "k_edge_radomized_response_partitioned() e="+all_epsilon+" e_1="+epsilon_q1+" rr_fall_back="+rr_fall_back;
		System.out.println(name);
		double start = System.currentTimeMillis();
		Random rand = new Random(123456);
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
				Random rand_node = new Random(node);
				int[] my_neighbors = new int[sanitized_number_of_outgoing_egdes];
				int[] fake_edges = new int[sanitized_number_of_outgoing_egdes];//
				
				for(int group = 0;group<sanitized_number_of_outgoing_egdes;group++) {
					int my_neighbor_in_group = get_group_edge_or_random(group, neighbor_list[node], sanitized_number_of_outgoing_egdes, edges_per_group, rand_node);
					my_neighbors[group] = my_neighbor_in_group;
					if(duplicate_check[my_neighbor_in_group]) {
						System.err.println("duplicate_check[my_neighbor_in_group] at "+my_neighbor_in_group);
					}
					duplicate_check[my_neighbor_in_group] = true;
					
					while(true){
						//get some random edge form this group
						int i = rand_node.nextInt(edges_per_group);
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
	 * Creates a edge-local DP version of the input Graph g. We assume that g is directed.
	 * @param g - non-private graph as input
	 * @param p - frac{1}{1+exp(\epsilon)}
	 * @return
	 */
	@Deprecated
	public static Graph k_edge_radomized_response(Graph g, final double p) {
		Random rand = new Random(123456);
		ArrayList<Integer>[] neighbor_list = g.get_neighbors();
		final double count_query_senstivity = 1;
		final double count_query_epsilon = 1;
		@SuppressWarnings("unchecked")
		ArrayList<Integer>[] sanitized_neighbor_list = new ArrayList[neighbor_list.length]; 
		
		for(int node=0;node<g.num_vertices;node++) {
			ArrayList<Integer> my_neighbors = neighbor_list[node];
			int sanitized_number_of_outgoing_egdes = sanitize(my_neighbors.size(), count_query_senstivity, count_query_epsilon);
			sanitized_number_of_outgoing_egdes = Math.min(sanitized_number_of_outgoing_egdes, g.num_vertices);//we do not have more vertices to sell
			ArrayList<Integer> my_sanitized_neighbors = new ArrayList<Integer>(sanitized_number_of_outgoing_egdes);
			
			for(int e=0;e<Math.min(my_neighbors.size(),sanitized_number_of_outgoing_egdes);e++) {
				int edge_target = my_neighbors.get(e);
				double dice = rand.nextDouble();
				if(dice <= (1-p)) {//with probability 1-p return the real value
					my_sanitized_neighbors.add(edge_target);
				}else {//dice a fake edge
					int fake_target;
					/*while(my_neighbors.contains((fake_target=rand.nextInt(g.num_vertices)))) {
						//ensures that the fake edge is not among the real edges
						//TODO maybe also ensure that it is not among the already reported ones
					}*/
					while(true) {
						fake_target=rand.nextInt(g.num_vertices);
						if(!my_neighbors.contains(fake_target)) {//ensures that the fake edge is not among the real edges
							if(!my_sanitized_neighbors.contains(fake_target)) {//or the already reported ones
								break;
							}
						}
					}
					my_sanitized_neighbors.add(fake_target);
				}
			}
			if(sanitized_number_of_outgoing_egdes>my_neighbors.size()){//Invent some fake edges
				int number_of_fake_edges = sanitized_number_of_outgoing_egdes-my_neighbors.size();
				for(int e=0;e<number_of_fake_edges;e++) {
					int fake_target;
					/*while(my_neighbors.contains((fake_target=rand.nextInt(g.num_vertices)))) {
						//ensures that the fake edge is not among the real edges
						//TODO maybe also ensure that it is not among the already reported ones
					}*/
					while(true) {
						fake_target=rand.nextInt(g.num_vertices);
						if(!my_neighbors.contains(fake_target)) {//ensures that the fake edge is not among the real edges
							if(!my_sanitized_neighbors.contains(fake_target)) {//or the already reported ones
								break;
							}
						}
					}
					my_sanitized_neighbors.add(fake_target);
				}
			}
			if(my_sanitized_neighbors.size()!=sanitized_number_of_outgoing_egdes) {
				System.err.println("k_edge_radomized_response() my_sanitized_neighbors.size()!=sanitized_number_of_outgoing_egdes");
			}
			sanitized_neighbor_list[node] = my_sanitized_neighbors;
		}
				
		return new Graph(sanitized_neighbor_list, g.name);
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
     * @return
     */
    static final int sanitize(final double val, final double lambda) {
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
	public static Graph top_k(Graph g, final double all_epsilon, final double epsilon_q1, boolean use_seq_composition) {
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
	 * top_k
	 * @param g
	 * @param p
	 * @return
	 */
	public static Graph educated_guess(Graph g, final double all_epsilon, final boolean non_private) {
		String name = "educated_guess()";
		System.out.println(name);
		Random rand = new Random(123456);
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
	
	static int random_choice(final double[] probabilities, Random rand) {
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
	public static double[] to_probabilities(int[] noissy_u) {
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
}
