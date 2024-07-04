package graphs;

import java.util.ArrayList;

import results.*;


public class TestGraphs {
	static int num_repitions = 3;	
	
	
	static void run_randomized_response() {
		ResultCollector.init();
		int graph_id = DataLoader.ADVOGATO;
		double[] all_eps = {1,2,3,4,5,6,7,8,9,10};
		
		Graph g = DataLoader.get_graph(graph_id);
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
	
	static void run_sequential_comp() {
		ResultCollector.init();
		int graph_id = DataLoader.ADVOGATO;	
		double[] all_eps = {1,2,3,4,5,6,7,8,9,10};
		
		Graph g = DataLoader.get_graph(graph_id);
		System.out.println(g);
		double e_q1 = 1;
		
		for(double epsilon  : all_eps) {
			System.out.println("**************eps="+epsilon+"\t");
			
			ArrayList<Graph> all_san_g = new ArrayList<Graph>(all_eps.length);
			for(int i=0;i<num_repitions;i++) {
				Graph san_g = Mechanism.k_edge_radomized_response_non_private_grouping(g, epsilon, e_q1, true);
				all_san_g.add(san_g);
			}
			ResultCollector.collect(g, all_san_g);
		}
		
		String name = "k-edge with seq compostion e_q1="+e_q1+" num_repitions="+num_repitions;
		ResultCollector.store(name, all_eps);
	}
	
	static void run_k_edge_non_private_grouping() {
		ResultCollector.init();
		int graph_id = DataLoader.ADVOGATO;
		
		double[] all_eps = {1,2,3,4,5,6,7,8,9,10};
		
		Graph g = DataLoader.get_graph(graph_id);
		System.out.println(g);
		double e_q1 = 1;
		
		for(double epsilon  : all_eps) {
			System.out.println("**************eps="+epsilon+"\t");
			
			ArrayList<Graph> all_san_g = new ArrayList<Graph>(all_eps.length);
			for(int i=0;i<num_repitions;i++) {
				Graph san_g = Mechanism.k_edge_radomized_response_non_private_grouping(g, epsilon, e_q1, false);
				all_san_g.add(san_g);
			}
			ResultCollector.collect(g, all_san_g);
		}
		
		String name = "k-edge with non-private grouping e_q1="+e_q1+" num_repitions="+num_repitions;
		ResultCollector.store(name, all_eps);
	}
	
	static void run_k_edge_grouping() {
		ResultCollector.init();
		int graph_id = DataLoader.ADVOGATO;
		
		double[] all_eps = {1,2,3,4,5,6,7,8,9,10};
		
		Graph g = DataLoader.get_graph(graph_id);
		System.out.println(g);
		double e_q1 = 1;
		
		for(double epsilon  : all_eps) {
			System.out.println("**************eps="+epsilon+"\t");
			
			ArrayList<Graph> all_san_g = new ArrayList<Graph>(all_eps.length);
			for(int i=0;i<num_repitions;i++) {
				Graph san_g = Mechanism.k_edge_radomized_response_partitioned(g, epsilon, e_q1);
				all_san_g.add(san_g);
			}
			ResultCollector.collect(g, all_san_g);
		}
		
		String name = "k-edge with private grouping e_q1="+e_q1+" num_repitions="+num_repitions;
		ResultCollector.store(name, all_eps);
	}
	

	public static void main(String[] args) {
		//algorithms.ProximityPrestige.run(DataLoader.get_advogato_graph());
		run_randomized_response();
		//run_k_edge_non_private_grouping();
		//run_sequential_comp();
		//run_k_edge_grouping();
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
