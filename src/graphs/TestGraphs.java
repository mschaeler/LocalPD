package graphs;

import java.util.ArrayList;
import java.util.Arrays;

import results.CentralityResult;
import results.Results;
import results.Statistics;

public class TestGraphs {
	static int num_repitions = 3;	
	
	
	static void run_randomized_response() {

		int graph_id = DataLoader.ADVOGATO;
		double[] all_eps = {1,2,3,4,5,6,7,8,9,10};
		
		Graph g = DataLoader.get_graph(graph_id);
		System.out.println(g);
		ArrayList<double[]> all_metrics = new ArrayList<double[]>();
		ArrayList<Double> all_centralities = new ArrayList<Double>();
		
		for(double epsilon  : all_eps) {
			System.out.print("**************eps="+epsilon+"\t");
			ArrayList<Graph> all_san_g = new ArrayList<Graph>(all_eps.length);
			for(int i=0;i<num_repitions;i++) {
				Graph san_g = Mechanism.radomized_response(g, epsilon);
				all_san_g.add(san_g);
			}
			double[] metrics = Metrics.statistics(g, all_san_g);
			all_metrics.add(metrics);
			double avg_centralities = Metrics.avg(Metrics.page_rank(g, all_san_g));
			all_centralities.add(avg_centralities);
		}
		
		String name = "Randomized response num_repitions="+num_repitions;
		new CentralityResult(name, all_eps, all_centralities);
		new Statistics(name, all_eps, all_metrics);
	}
	
	static void run_sequential_comp() {
		
		int graph_id = DataLoader.ADVOGATO;	
		double[] all_eps = {1,2,3,4,5,6,7,8,9,10};
		
		Graph g = DataLoader.get_graph(graph_id);
		System.out.println(g);
		ArrayList<double[]> all_metrics = new ArrayList<double[]>();
		ArrayList<Double> all_centralities = new ArrayList<Double>();
		double e_q1 = 1;
		
		for(double epsilon  : all_eps) {
			System.out.println("**************eps="+epsilon+"\t");
			
			ArrayList<Graph> all_san_g = new ArrayList<Graph>(all_eps.length);
			for(int i=0;i<num_repitions;i++) {
				Graph san_g = Mechanism.k_edge_radomized_response_non_private_grouping(g, epsilon, e_q1, true);
				all_san_g.add(san_g);
			}
			double[] metrics = Metrics.statistics(g, all_san_g);
			all_metrics.add(metrics);
			double avg_centralities = Metrics.avg(Metrics.page_rank(g, all_san_g));
			all_centralities.add(avg_centralities);
		}
		
		String name = "k-edge with seq compostion e_q1="+e_q1+" num_repitions="+num_repitions;
		new CentralityResult(name, all_eps, all_centralities);
		new Statistics(name, all_eps, all_metrics);
	}
	
	static void run_k_edge_non_private_grouping() {
		int graph_id = DataLoader.ADVOGATO;
		
		double[] all_eps = {1,2,3,4,5,6,7,8,9,10};
		
		Graph g = DataLoader.get_graph(graph_id);
		System.out.println(g);
		ArrayList<double[]> all_metrics = new ArrayList<double[]>();
		ArrayList<Double> all_centralities = new ArrayList<Double>();
		double e_q1 = 1;
		
		for(double epsilon  : all_eps) {
			System.out.println("**************eps="+epsilon+"\t");
			
			ArrayList<Graph> all_san_g = new ArrayList<Graph>(all_eps.length);
			for(int i=0;i<num_repitions;i++) {
				Graph san_g = Mechanism.k_edge_radomized_response_non_private_grouping(g, epsilon, e_q1, false);
				all_san_g.add(san_g);
			}
			double[] metrics = Metrics.statistics(g, all_san_g);
			all_metrics.add(metrics);
			double avg_centralities = Metrics.avg(Metrics.page_rank(g, all_san_g));
			all_centralities.add(avg_centralities);
		}
		
		String name = "k-edge with non-private grouping e_q1="+e_q1+" num_repitions="+num_repitions;
		new CentralityResult(name, all_eps, all_centralities);
		new Statistics(name, all_eps, all_metrics);
	}
	
	static void run_k_edge_grouping() {
		int graph_id = DataLoader.ADVOGATO;
		
		double[] all_eps = {1,2,3,4,5,6,7,8,9,10};
		
		Graph g = DataLoader.get_graph(graph_id);
		System.out.println(g);
		ArrayList<double[]> all_metrics = new ArrayList<double[]>();
		ArrayList<Double> all_centralities = new ArrayList<Double>();
		double e_q1 = 1;
		
		for(double epsilon  : all_eps) {
			System.out.println("**************eps="+epsilon+"\t");
			
			ArrayList<Graph> all_san_g = new ArrayList<Graph>(all_eps.length);
			for(int i=0;i<num_repitions;i++) {
				Graph san_g = Mechanism.k_edge_radomized_response_partitioned(g, epsilon, e_q1);
				all_san_g.add(san_g);
			}
			double[] metrics = Metrics.statistics(g, all_san_g);
			all_metrics.add(metrics);
			double avg_centralities = Metrics.avg(Metrics.page_rank(g, all_san_g));
			all_centralities.add(avg_centralities);
		}
		
		String name = "k-edge with private grouping e_q1="+e_q1+" num_repitions="+num_repitions;
		//Metrics.out(name, all_eps, all_metrics);
		new CentralityResult(name, all_eps, all_centralities);
		new Statistics(name, all_eps, all_metrics);
	}
	

	public static void main(String[] args) {
		run_randomized_response();
		run_k_edge_non_private_grouping();
		run_sequential_comp();
		run_k_edge_grouping();
		Results.all_out();
		
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
