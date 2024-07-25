package graphs;

import java.util.ArrayList;

import results.*;


public class TestGraphs {
	static int num_repitions = 3;	
	
	
	static void run_randomized_response(Graph g) {
		ResultCollector.init();
		double[] all_eps = {1,2,3,4,5,6,7,8,9,10};
		System.out.println(g);
		
		for(double epsilon  : all_eps) {
			System.out.print("**************eps="+epsilon+"\t");
			ArrayList<Graph> all_san_g = new ArrayList<Graph>(all_eps.length);
			for(int i=0;i<num_repitions;i++) {
				Graph san_g = Mechanism.radomized_response(g, epsilon);
				all_san_g.add(san_g);
			}
			ResultCollector.collect(g, all_san_g);
		}
		
		String name = "Randomized response num_repitions="+num_repitions;
		ResultCollector.store(name, all_eps);
	}
	
	static void run_sequential_comp(Graph g) {
		ResultCollector.init();	
		double[] all_eps = {1,2,3,4,5,6,7,8,9,10};
		System.out.println(g);
		double e_q1 = 1;
		
		for(double epsilon  : all_eps) {
			System.out.println("**************eps="+epsilon+"\t");
			
			ArrayList<Graph> all_san_g = new ArrayList<Graph>(all_eps.length);
			for(int i=0;i<num_repitions;i++) {
				Graph san_g = Mechanism.k_edge_radomized_response_non_private_grouping_rr_fallback(g, epsilon, e_q1, true);
				all_san_g.add(san_g);
			}
			ResultCollector.collect(g, all_san_g);
		}
		
		String name = "k-edge with seq compostion (with fall_back) e_q1="+e_q1+" num_repitions="+num_repitions;
		ResultCollector.store(name, all_eps);
	}
	
	static void run_k_edge_non_private_grouping(Graph g) {
		ResultCollector.init();
		double[] all_eps = {1,2,3,4,5,6,7,8,9,10};
		System.out.println(g);
		double e_q1 = 1;
		
		for(double epsilon  : all_eps) {
			System.out.println("**************eps="+epsilon+"\t");
			
			ArrayList<Graph> all_san_g = new ArrayList<Graph>(all_eps.length);
			for(int i=0;i<num_repitions;i++) {
				Graph san_g = Mechanism.k_edge_radomized_response_non_private_grouping_rr_fallback(g, epsilon, e_q1, false);
				all_san_g.add(san_g);
			}
			ResultCollector.collect(g, all_san_g);
		}
		
		String name = "k-edge with non-private grouping e_q1="+e_q1+" num_repitions="+num_repitions;
		ResultCollector.store(name, all_eps);
	}
	
	static void run_k_edge_grouping(Graph g) {
		ResultCollector.init();
		double[] all_eps = {1,2,3,4,5,6,7,8,9,10};
		System.out.println(g);
		double e_q1 = 1;
		
		for(double epsilon  : all_eps) {
			System.out.println("**************eps="+epsilon+"\t");
			
			ArrayList<Graph> all_san_g = new ArrayList<Graph>(all_eps.length);
			for(int i=0;i<num_repitions;i++) {
				Graph san_g = Mechanism.k_edge_radomized_response_partitioned(g, epsilon, e_q1, true);
				all_san_g.add(san_g);
			}
			ResultCollector.collect(g, all_san_g);
		}
		
		String name = "k-edge with private grouping e_q1="+e_q1+" num_repitions="+num_repitions;
		ResultCollector.store(name, all_eps);
	}
	
	static final int RANDOM_RESPONSE 	= 0;
	static final int K_EDGE_NON_PRIVATE = 1;
	static final int K_EDGE_SEQ 		= 2;
	static final int K_EDGE_SEQ_RR		= 3;
	static final int K_EDGE_PART 		= 4;
	static final int K_EDGE_PART_RR 	= 5;
	static final int TOP_K			 	= 6;
	static final int GUESS			 	= 7;

	static double e_q1 = 1.0;
	
	static Graph run(final Graph g, final int mechanism, final double epsilon) {
		Graph san_g;
		if(mechanism==RANDOM_RESPONSE) {
			san_g = Mechanism.radomized_response(g, epsilon);
		}else if(mechanism==K_EDGE_NON_PRIVATE){
			san_g = Mechanism.k_edge_radomized_response_non_private_grouping_rr_fallback(g, epsilon, e_q1, false);
		}else if(mechanism==K_EDGE_SEQ){
			san_g = Mechanism.k_edge_radomized_response_non_private_grouping(g, epsilon, e_q1, true);
		}else if(mechanism==K_EDGE_SEQ_RR){
			san_g = Mechanism.k_edge_radomized_response_non_private_grouping_rr_fallback(g, epsilon, e_q1, true);
		}else if(mechanism==K_EDGE_PART){
			san_g = Mechanism.k_edge_radomized_response_partitioned(g, epsilon, e_q1, false);
		}else if(mechanism==K_EDGE_PART_RR){
			san_g = Mechanism.k_edge_radomized_response_partitioned(g, epsilon, e_q1, true);
		}else if(mechanism==TOP_K){
			san_g = Mechanism.top_k(g, epsilon, e_q1, false);
		}else if(mechanism==GUESS){
			san_g = Mechanism.educated_guess(g, epsilon, true);
		}else{
			System.err.println("Unknown mechanism "+mechanism);
			san_g = null;
		}
		return san_g;
	}
	
	static void run(Graph g, int[] all_mechanisms) {
		for(int mechanism : all_mechanisms) {
			run(g, mechanism);
		}
		results.Results.all_out();
	}
	
	static void run(Graph g, int mechanism) {
		ResultCollector.init();
		double[] all_eps = {1,2,3,4,5,6,7,8,9,10};
		//double[] all_eps = {5,6,7,8,9,10};
		System.out.println(g);
		double e_q1 = 1;
		
		for(double epsilon  : all_eps) {
			System.out.println("**************eps="+epsilon+"\t");
			
			ArrayList<Graph> all_san_g = new ArrayList<Graph>(all_eps.length);
			for(int i=0;i<num_repitions;i++) {
				Graph san_g = run(g, mechanism, epsilon);
				all_san_g.add(san_g);
			}
			ResultCollector.collect(g, all_san_g);
		}
		
		//String name = "k-edge with private grouping e_q1="+e_q1+" num_repitions="+num_repitions;
		String name = name(mechanism)+" e_1="+e_q1+" num_repitions="+num_repitions;
		ResultCollector.store(name, all_eps);
	}
	
	static String name(final int mechanism) {
		if(mechanism==RANDOM_RESPONSE) {
			return "RANDOM_RESPONSE";
		}else if(mechanism==K_EDGE_NON_PRIVATE){
			return "K_EDGE_NON_PRIVATE";
		}else if(mechanism==K_EDGE_SEQ){
			return "K_EDGE_SEQ";
		}else if(mechanism==K_EDGE_SEQ_RR){
			return "K_EDGE_SEQ_RR";
		}else if(mechanism==K_EDGE_PART){
			return "K_EDGE_PART";
		}else if(mechanism==K_EDGE_PART_RR){
			return "K_EDGE_PART_RR";
		}else if(mechanism==TOP_K){
			return "TOP_K";
		}else if(mechanism==GUESS){
			return "GUESS";
		}else{
			System.err.println("Unknown mechanism "+mechanism);
			return null;
		}
	}

	public static void main(String[] args) {
		//Config.USE_RESULT_STATISTICS = true;
		Config.USE_RESULT_PROXIMITY_PRESTIGE = true;
		Graph g = DataLoader.get_graph(DataLoader.ADVOGATO);
		//int[] all_mechanisms = {K_EDGE_NON_PRIVATE, K_EDGE_SEQ, K_EDGE_SEQ_RR, K_EDGE_PART, K_EDGE_PART_RR};
		int[] all_mechanisms = {GUESS};
		
		//Metrics.proximity_prestige(g);
		run(g, all_mechanisms);
		//run_randomized_response(g);
		//run_k_edge_non_private_grouping(g);
		//run_sequential_comp(g);
		//run_k_edge_grouping(g);
		results.Results.all_out();
		
		/*Graph example = Graph.get_example();
		System.out.println(example);
		Mechanism.two_k_series(example);
		int k = 0;
		while(k++<20)
			Mechanism.k_edge_radomized_response(example, 0.5);
		
		
		System.out.println(Mechanism.k_edge_radomized_response_non_private_grouping(example, 8, 1, true));
		
		Graph g = DataLoader.get_advogato_graph();
		System.out.println(g);
		
		System.out.println(Mechanism.baseline(g, 8, 1, true));
		System.out.println(Mechanism.baseline(g, 8, 1, true));
		System.out.println(Mechanism.baseline(g, 8, 1, true));
		System.out.println(Mechanism.baseline(g, 8, 1, true));
		
		System.out.println(Mechanism.baseline(g, 8, 1, false));
		System.out.println(Mechanism.baseline(g, 8, 1, false));
		System.out.println(Mechanism.baseline(g, 8, 1, false));
		System.out.println(Mechanism.baseline(g, 8, 1, false));
		
		g = DataLoader.get_enron_graph();
		g = Graph.dedup_edges(g);
		
		System.out.println(Mechanism.baseline(g, 8, 1, true));
		System.out.println(Mechanism.baseline(g, 8, 1, true));
		System.out.println(Mechanism.baseline(g, 8, 1, true));
		System.out.println(Mechanism.baseline(g, 8, 1, true));
		
		System.out.println(Mechanism.baseline(g, 8, 1, false));
		System.out.println(Mechanism.baseline(g, 8, 1, false));
		System.out.println(Mechanism.baseline(g, 8, 1, false));
		System.out.println(Mechanism.baseline(g, 8, 1, false));
		
		/*System.out.println(Mechanism.k_edge_radomized_response_non_private_grouping(g, 8, 1, true));
		System.out.println(Mechanism.k_edge_radomized_response_non_private_grouping(g, 8, 1, true));
		System.out.println(Mechanism.k_edge_radomized_response_non_private_grouping(g, 8, 1, true));
		System.out.println(Mechanism.k_edge_radomized_response_non_private_grouping(g, 8, 1, true));
		
		System.out.println(Mechanism.k_edge_radomized_response_non_private_grouping(g, 8, 1, false));
		System.out.println(Mechanism.k_edge_radomized_response_non_private_grouping(g, 8, 1, false));
		System.out.println(Mechanism.k_edge_radomized_response_non_private_grouping(g, 8, 1, false));
		System.out.println(Mechanism.k_edge_radomized_response_non_private_grouping(g, 8, 1, false));*
		
		System.out.println(Mechanism.radomized_response(g, 8));
		System.out.println(Mechanism.radomized_response(g, 8));
		System.out.println(Mechanism.radomized_response(g, 8));
		System.out.println(Mechanism.radomized_response(g, 8));
		
		
		
		
		//Graph g = DataLoader.get_enron_graph();
		//Graph g_ = Graph.dedup_edges(g);
		
		//System.out.println(Mechanism.radomized_response(g_, 8));
		
		//Graph e_small = DataLoader.get_enron_graph_employee_unique_edges();
		//System.out.println(e_small);
		/*int counter = 0;
		while(counter++<10)
			System.out.println(Mechanism.k_edge_radomized_response_non_private_grouping(g_, 8, 1, true));
		counter = 0;
		while(counter++<10)
			System.out.println(Mechanism.k_edge_radomized_response_non_private_grouping(g_, 8, 1, false));
		while(counter++<10)
			System.out.println(Mechanism.radomized_response(g_, 8));
		*/
	}

}
