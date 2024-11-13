package graphs;

import java.util.ArrayList;
import java.util.HashSet;

import cern.colt.Arrays;
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
	static final int LDP_GEN			= 8;
	static final int CHUN_LU			= 9;
	static final int TWO_K_SERIES		= 10;
	static final int M_PART				= 11;
	static final int M_SAMPLE			= 12;
	static final int M_SAMPLE_WEIGHTED  = 13;
	static final int M_PART_2			= 14;
	static final int M_NAIVE			= 15;
	static final int M_PART_NO_SEQ		= 16;
	

	static double e_q1 = 1.0;
	
	public static Graph run(final Graph g, final int mechanism, final double epsilon) {
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
			san_g = Mechanism.top_k(g, epsilon, e_q1, true);
		}else if(mechanism==GUESS){
			san_g = Mechanism.educated_guess(g, epsilon, false);
		}else if(mechanism==LDP_GEN){
			boolean non_private = false;
			san_g = Mechanism.LDPGen(g, epsilon, non_private);
		}else if(mechanism==TWO_K_SERIES){
			boolean non_private = true;
			san_g = Mechanism.two_k_series(g, epsilon, non_private);
		}else if(mechanism==CHUN_LU){
			boolean non_private = false;
			boolean grouped = true;
			san_g = Mechanism.Chun_Lu_Model(g, epsilon, non_private, grouped);
		}else if(mechanism==M_PART){
			san_g = Mechanism.m_part(g, epsilon, e_q1);
		}else if(mechanism==M_PART_2){
			san_g = Mechanism.m_part_2(g, epsilon, e_q1);
		}else if(mechanism==M_SAMPLE){
			san_g = Mechanism.m_sample(g, epsilon, e_q1);
		}else if(mechanism==M_SAMPLE_WEIGHTED){
			san_g = Mechanism.m_sample_weighted(g, epsilon, e_q1);
		}else if(mechanism==M_NAIVE){
			san_g = Mechanism.m_straw_man(g, epsilon);
		}else if(mechanism==M_PART_NO_SEQ){
			san_g = Mechanism.m_part_no_seq(g, epsilon);
		}else{
			System.err.println("run() Unknown mechanism "+mechanism);
			san_g = null;
		}
		return san_g;
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
		}else if(mechanism==LDP_GEN){
			return "LDP_GEN";
		}else if(mechanism==CHUN_LU){
			return "CHUN_LU";
		}else if(mechanism==TWO_K_SERIES){
			return "TWO_K_SERIES";
		}else if(mechanism==M_PART){
			return "M_PART";
		}else if(mechanism==M_PART_2){
			return "M_PART_2";
		}else if(mechanism==M_SAMPLE){
			return "M_SAMPLE";
		}else if(mechanism==M_SAMPLE_WEIGHTED){
			return "M_SAMPLE_WEIGHTED";
		}else if(mechanism==M_NAIVE){
			return "M_NAIVE";
		}else if(mechanism==M_PART_NO_SEQ){
			return "M_PART_NO_SEQ";
		}else{
			System.err.println("name() Unknown mechanism "+mechanism);
			return null;
		}
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
	
	static void matrialize_private_graphs(int[] graphs, int[] mechanism, int num_repitions, double[] all_eps) {
		for(int g_id : graphs) {
			Graph g = DataLoader.get_graph(g_id);
			for(int m : mechanism) {
				for(double epsilon : all_eps) {
					for(int i=0;i<num_repitions;i++) {
						Graph san_g = run(g, m, epsilon);
						san_g.to_file(san_g.to_edge_list(), name(m), DataLoader.get_graph_name(g_id), epsilon ,i);	
					}
				}
			}
		}
	}
	
	static void graph_proximity_prestige(int graph_enum, String[] all_m, double[] all_e) {
		System.out.println("***graph_proximity_prestige("+DataLoader.get_graph_name(graph_enum)+") " + Arrays.toString(all_m)+" "+Arrays.toString(all_e));
		
		ArrayList<MaterializedGraphs> mg_s = DataLoader.get_sanitized_graphs(graph_enum);
		ArrayList<ArrayList<String[]>> all_results_i = new ArrayList<ArrayList<String[]>>();
		ArrayList<ArrayList<String[]>> all_results_dist = new ArrayList<ArrayList<String[]>>();
		ArrayList<ArrayList<String[]>> all_results_pp = new ArrayList<ArrayList<String[]>>();
		
		for(String approach : all_m) {
			System.out.println(approach);
			ArrayList<MaterializedGraphs> graphs_of_mechanism = MaterializedGraphs.filter_by_mechanism(mg_s, approach);
			ArrayList<String[]> out_i = new ArrayList<String[]>();
			ArrayList<String[]> out_dist = new ArrayList<String[]>();
			ArrayList<String[]> out_pp = new ArrayList<String[]>();
			
			for(double epsilon : all_e) {
				System.out.println(epsilon);
				ArrayList<MaterializedGraphs> temp = MaterializedGraphs.filter_by_epsilon(graphs_of_mechanism, epsilon);
				ArrayList<Graph> g_s = DataLoader.get_sanitized_graphs(temp);
				double[] res = Metrics.proximity_prestige(g_s);
				System.out.println(Arrays.toString(res));
				String[] line_i = {""+epsilon, ""+(res[0])};
				out_i.add(line_i);
				String[] line_dist = {""+epsilon, ""+(res[1])};
				out_dist.add(line_dist);
				String[] line_pp = {""+epsilon, ""+(res[2])};
				out_pp.add(line_pp);
			}
			System.out.println(approach);
			for(String[] arr : out_i) {
				System.out.println(Arrays.toString(arr));
			}
			for(String[] arr : out_dist) {
				System.out.println(Arrays.toString(arr));
			}
			for(String[] arr : out_pp) {
				System.out.println(Arrays.toString(arr));
			}
			all_results_i.add(out_i);
			all_results_dist.add(out_dist);
			all_results_pp.add(out_pp);
		}
		String header = "eps\t";
		for(String approach : all_m) {
			header+=approach+"\t";
		}
		System.out.println(DataLoader.get_graph_name(graph_enum)+" Avg proximity prestige: |I|");
		System.out.println(header);
		int num_lines = all_results_i.get(0).size();
		for(int line=0;line<num_lines;line++) {
			System.out.print(all_results_i.get(0).get(line)[0]+"\t");//Epsilon
			for(ArrayList<String[]> r : all_results_i) {
				System.out.print(r.get(line)[1]+"\t");
			}
			System.out.println();
		}
		System.out.println(DataLoader.get_graph_name(graph_enum)+" Avg proximity prestige: avg(dist())");
		System.out.println(header);
		num_lines = all_results_dist.get(0).size();
		for(int line=0;line<num_lines;line++) {
			System.out.print(all_results_dist.get(0).get(line)[0]+"\t");//Epsilon
			for(ArrayList<String[]> r : all_results_dist) {
				System.out.print(r.get(line)[1]+"\t");
			}
			System.out.println();
		}
		System.out.println(DataLoader.get_graph_name(graph_enum)+" Avg proximity prestige: pp");
		System.out.println(header);
		num_lines = all_results_pp.get(0).size();
		for(int line=0;line<num_lines;line++) {
			System.out.print(all_results_pp.get(0).get(line)[0]+"\t");//Epsilon
			for(ArrayList<String[]> r : all_results_pp) {
				System.out.print(r.get(line)[1]+"\t");
			}
			System.out.println();
		}
	}
	
	static void graph_statistics(int graph_enum, String[] all_m, double[] all_e) {
		System.out.println("***graph_statistics("+DataLoader.get_graph_name(graph_enum)+") " + Arrays.toString(all_m)+" "+Arrays.toString(all_e));
		Graph g = DataLoader.get_graph(graph_enum);
		ArrayList<MaterializedGraphs> mg_s = DataLoader.get_sanitized_graphs(graph_enum);
		if(mg_s.isEmpty()) {
			System.err.println("graph_statistics() graph enum resutled in no graphs "+graph_enum);
		}
		ArrayList<ArrayList<String[]>> all_results = new ArrayList<ArrayList<String[]>>();
		for(String approach : all_m) {
			System.out.println(approach);
			ArrayList<MaterializedGraphs> graphs_of_mechanism = MaterializedGraphs.filter_by_mechanism(mg_s, approach);
			ArrayList<String[]> out = new ArrayList<String[]>();
			for(double epsilon : all_e) {
				System.out.println(epsilon);
				ArrayList<MaterializedGraphs> temp = MaterializedGraphs.filter_by_epsilon(graphs_of_mechanism, epsilon);
				ArrayList<Graph> g_s = DataLoader.get_sanitized_graphs(temp);
				double[] res = Metrics.avg_edge_edit_dist(g, g_s);
				System.out.println(Arrays.toString(res));
				String[] line = {""+epsilon, ""+(res[Metrics.missing_edge]+res[Metrics.new_fake_edge])};
				out.add(line);
			}
			System.out.println(approach);
			for(String[] arr : out) {
				System.out.println(Arrays.toString(arr));
			}
			all_results.add(out);
		}
		String header = "eps\t";
		for(String approach : all_m) {
			header+=approach+"\t";
		}
		System.out.println(g.name+" Avg edge-edit distance");
		System.out.println(header);
		int num_lines = all_results.get(0).size();
		for(int line=0;line<num_lines;line++) {
			System.out.print(all_results.get(0).get(line)[0]+"\t");//Epsilon
			for(ArrayList<String[]> r : all_results) {
				System.out.print(r.get(line)[1]+"\t");
			}
			System.out.println();
		}
	}
	
	static void count_triangles(int graph_enum, String[] all_m, double[] all_e) {
		System.out.println("***count_tringle_deltas("+DataLoader.get_graph_name(graph_enum)+") " + Arrays.toString(all_m)+" "+Arrays.toString(all_e));
		Graph g = DataLoader.get_graph(graph_enum);
		ArrayList<MaterializedGraphs> mg_s = DataLoader.get_sanitized_graphs(graph_enum);
		ArrayList<ArrayList<String[]>> all_results = new ArrayList<ArrayList<String[]>>();
		for(String approach : all_m) {
			System.out.println(approach);
			ArrayList<MaterializedGraphs> graphs_of_mechanism = MaterializedGraphs.filter_by_mechanism(mg_s, approach);
			ArrayList<String[]> out = new ArrayList<String[]>();
			for(double epsilon : all_e) {
				System.out.println(epsilon);
				ArrayList<MaterializedGraphs> temp = MaterializedGraphs.filter_by_epsilon(graphs_of_mechanism, epsilon);
				ArrayList<Graph> g_s = DataLoader.get_sanitized_graphs(temp);
				double[] res = Metrics.count_triangles(g, g_s);//XXX
				//double[] res = Metrics.count_triangle_delta(g, g_s);
				System.out.println(Arrays.toString(res));
				String[] line = {""+epsilon, ""+res[0], ""+res[1]};
				out.add(line);
			}
			System.out.println(approach);
			for(String[] arr : out) {
				System.out.println(Arrays.toString(arr));
			}
			all_results.add(out);
		}
		String header = "eps\tOriginal\t";
		for(String approach : all_m) {
			header+=approach+"\t";
		}
		System.out.println(g.name);
		System.out.println(header);
		int num_lines = all_results.get(0).size();
		for(int line=0;line<num_lines;line++) {
			System.out.print(all_results.get(0).get(line)[0]+"\t");
			System.out.print(all_results.get(0).get(line)[1]+"\t");
			for(ArrayList<String[]> r : all_results) {
				System.out.print(r.get(line)[2]+"\t");
			}
			System.out.println();
		}
	}
	
	static void count_tringles(int graph_enum, String[] all_m) {
		Graph g = DataLoader.get_graph(graph_enum);
		ArrayList<MaterializedGraphs> mg_s = DataLoader.get_sanitized_graphs(graph_enum);
		double[] all_e = MaterializedGraphs.get_all_epsilons(mg_s);
		count_triangles(graph_enum, all_m, all_e);
	}
	
	static void count_tringles(int graph_enum) {
		Graph g = DataLoader.get_graph(graph_enum);
		ArrayList<MaterializedGraphs> mg_s = DataLoader.get_sanitized_graphs(graph_enum);
		String[] all_m = MaterializedGraphs.get_all_mechanisms(mg_s);
		double[] all_e = MaterializedGraphs.get_all_epsilons(mg_s);
		count_triangles(graph_enum, all_m, all_e);
	}

	public static void out_graph_statistics() {
		int[] graph_enums = {DataLoader.CONGRESS_TWITTER, DataLoader.ADVOGATO, DataLoader.ENRON_SINGLE_EDGE};
		
		for(int grah_enum : graph_enums) {
			Graph g = DataLoader.get_graph(grah_enum);
			System.out.println("Graph enum "+grah_enum);
			System.out.println(g.name+"\t"+g.num_vertices+"\t"+g.int_num_edges()+"\t"+((double)g.int_num_edges()/(double)g.num_vertices)+"\t"+Metrics.count_triangles(g));
		}
		
		
	}
	
	public static void experiment_local_error() {
		int[] all_mechanism = {M_NAIVE, M_PART_2, M_SAMPLE_WEIGHTED, RANDOM_RESPONSE, LDP_GEN, CHUN_LU, TWO_K_SERIES};
		//int[] all_mechanism = {M_PART_2};
		//int[] all_graphs = {DataLoader.CONGRESS_TWITTER, DataLoader.ADVOGATO, DataLoader.ENRON_SINGLE_EDGE};
		int[] all_graphs = {DataLoader.ENRON_SINGLE_EDGE};
		double[] all_epsilon = {1,2,3,4,5,6,7,8,9,10};
		experiment_local_error(all_mechanism, all_epsilon, all_graphs);	
	}
	
	private static void experiment_local_error(int[] all_mechanism, double[] all_epsilon, int[] all_graphs) {
		experiment_local_error(get_mechanism_names(all_mechanism), all_epsilon, all_graphs);
	}

	static void experiment_local_error(String[] all_mechanism, double[] all_epsilon, int[] all_graphs) {
		for(int graph_id  : all_graphs) {
			System.out.println("Graph id="+graph_id+" "+DataLoader.get_graph_name(graph_id));
			graph_statistics(graph_id, all_mechanism, all_epsilon);
		}
	}
	
	static String[] get_mechanism_names(int[] all_mechanism){
		String[] names = new String[all_mechanism.length];
		for(int i=0;i<all_mechanism.length;i++) {
			names[i] = name(all_mechanism[i]);
		}
		return names;
	}
	
	public static void matrialize_private_graphs(){
		int[] graphs = {DataLoader.ENRON_SINGLE_EDGE};
		//int[] mechanism = {K_EDGE_NON_PRIVATE, K_EDGE_SEQ_RR, K_EDGE_PART_RR, TOP_K, GUESS, RANDOM_RESPONSE};
		//int[] mechanism = {M_PART_2};
		int[] mechanism = {M_PART_2, M_SAMPLE_WEIGHTED, RANDOM_RESPONSE, LDP_GEN, CHUN_LU, TWO_K_SERIES};
		int num_repitions = 10;
		double[] all_eps = {1,2,3,4,5,6,7,8,9,10};
		matrialize_private_graphs(graphs, mechanism, num_repitions, all_eps);
	}
	
	public static void main(String[] args) {
		//out_graph_statistics();
		//experiment_local_error();	
		matrialize_private_graphs();

		System.exit(0);
		//String[] all_mechanisms = {name(LDP_GEN), name(CHUN_LU), name(TWO_K_SERIES), name(M_PART)};
		String[] all_mechanisms = {name(M_PART_NO_SEQ), name(M_PART), name(M_SAMPLE_WEIGHTED)};
		//String[] all_mechanisms = {name(K_EDGE_NON_PRIVATE), name(K_EDGE_SEQ_RR), name(K_EDGE_PART_RR), name(TOP_K), name(GUESS), name(RANDOM_RESPONSE)};
		double[] all_eps = {1,2,3,4,5,6,7,8,9,10};
		//graph_statistics(DataLoader.ADVOGATO, all_mechanisms, all_eps);
		count_triangles(DataLoader.CONGRESS_TWITTER, all_mechanisms, all_eps);
		//graph_proximity_prestige(DataLoader.ADVOGATO, all_mechanisms, all_eps);
		System.exit(0);
		//Config.USE_RESULT_STATISTICS = true;
		ArrayList<MaterializedGraphs> mg_s = DataLoader.get_sanitized_graphs(DataLoader.ADVOGATO);
		for(MaterializedGraphs mg : mg_s) {
			DataLoader.load(mg);	
		}
				
		Config.USE_RESULT_PROXIMITY_PRESTIGE = true;
		Graph g = DataLoader.get_graph(DataLoader.ENRON_SINGLE_EDGE);
		//g.to_file();
		//int[] all_mechanisms = {K_EDGE_NON_PRIVATE, K_EDGE_SEQ, K_EDGE_SEQ_RR, K_EDGE_PART, K_EDGE_PART_RR, TOP_K, GUESS};
		//int[] all_mechanisms = {TOP_K};
		
		
		//Metrics.proximity_prestige(g);
		//Metrics.page_rank(g);
		//run(g, all_mechanisms);
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
