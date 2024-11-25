package graphs;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;

import algorithms.ProximityPrestige;
import algorithms.ProximityPrestige.InfDom;
import misc.Histogram;
import misc.Util;
import results.*;


public class TestGraphs {
	static int num_repitions = Config.num_repitions;	
	
	public static final int RANDOM_RESPONSE 	= 0;
	public static final int K_EDGE_NON_PRIVATE 	= 1;
	public static final int K_EDGE_SEQ 			= 2;
	public static final int K_EDGE_SEQ_RR		= 3;
	public static final int K_EDGE_PART 		= 4;
	public static final int K_EDGE_PART_RR 		= 5;
	public static final int TOP_K			 	= 6;
	public static final int GUESS			 	= 7;
	public static final int LDP_GEN				= 8;
	public static final int CHUNG_LU			= 9;
	public static final int TWO_K_SERIES		= 10;
	public static final int M_PART				= 11;
	public static final int M_SAMPLE			= 12;
	public static final int M_SAMPLE_WEIGHTED  	= 13;
	public static final int M_PART_2			= 14;
	public static final int M_NAIVE				= 15;
	public static final int M_PART_NO_SEQ		= 16;
	public static final int M_SAMPLE_NO_RR		= 17;
	public static final int M_SAMPLE_WEIGHTED_NO_RR = 18;
	
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
		}else if(mechanism==CHUNG_LU){
			boolean none_private = false;
			boolean grouped = false;
			san_g = Mechanism.Chung_Lu_Model(g, epsilon, none_private, grouped);
		}else if(mechanism==M_PART){
			san_g = Mechanism.m_part(g, epsilon, e_q1);
		}else if(mechanism==M_PART_2){
			san_g = Mechanism.m_part_2(g, epsilon, e_q1);
		}else if(mechanism==M_SAMPLE){
			boolean no_rr_fall_back = false;
			san_g = Mechanism.m_sample(g, epsilon, e_q1, no_rr_fall_back);
		}else if(mechanism==M_SAMPLE_WEIGHTED){
			boolean no_rr_fall_back = false;
			san_g = Mechanism.m_sample_weighted(g, epsilon, e_q1, no_rr_fall_back);
		}else if(mechanism==M_NAIVE){
			san_g = Mechanism.m_straw_man(g, epsilon);
		}else if(mechanism==M_PART_NO_SEQ){
			san_g = Mechanism.m_part_no_seq(g, epsilon);
		}else if(mechanism==M_SAMPLE_NO_RR){
			boolean no_rr_fall_back = true;
			san_g = Mechanism.m_sample(g, epsilon, e_q1, no_rr_fall_back);
		}else if(mechanism==M_SAMPLE_WEIGHTED_NO_RR){
			boolean no_rr_fall_back = true;
			san_g = Mechanism.m_sample_weighted(g, epsilon, e_q1, no_rr_fall_back);
		}else{
			System.err.println("run() Unknown mechanism "+mechanism);
			san_g = null;
		}
		return san_g;
	}
	
	private static String name(final int mechanism) {
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
		}else if(mechanism==CHUNG_LU){
			return "CHUNG_LU";
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
		}else if(mechanism==M_SAMPLE_NO_RR){
			return "M_SAMPLE_NO_RR";
		}else if(mechanism==M_SAMPLE_WEIGHTED_NO_RR){
			return "M_SAMPLE_WEIGHTED_NO_RR";
		}else{
			System.err.println("name() Unknown mechanism "+mechanism);
			return null;
		}
	}
	
	private static void matrialize_private_graphs(int[] graphs, int[] mechanism, int num_repitions, double[] all_eps) {
		Config.materialize_graph = true;
		for(int g_id : graphs) {
			Graph g = DataLoader.get_graph(g_id);
			String path = "./data/synthetic/"+DataLoader.get_graph_name(g_id);
			File folder = new File(path);
			if(!folder.exists()) {
				folder.mkdirs();
			}
			
			for(int m : mechanism) {
				Mechanism.rand.setSeed(Util.seed);
				for(double epsilon : all_eps) {
					for(int i=0;i<num_repitions;i++) {
						try (PrintWriter out = new PrintWriter(path+"/"+name(m)+"-"+epsilon+"-"+i+".edgelist")) {
							Graph.print_writer = out;
							run(g, m, epsilon);
							Graph.print_writer = null;
						} catch (FileNotFoundException e) {
							e.printStackTrace();
						}	
					}
				}
			}
		}
		Config.materialize_graph = false;
	}
	
	private static void graph_statistics(int graph_enum, String[] all_m, double[] all_e) {
		System.out.println("***graph_statistics("+DataLoader.get_graph_name(graph_enum)+") " + Arrays.toString(all_m)+" "+Arrays.toString(all_e));
		Graph g = DataLoader.get_graph(graph_enum);
		BitSet[] ground_truth = g.get_adjancency_matrix_as_bit_vector();
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
				double[] res = Metrics.avg_edge_edit_dist(ground_truth,g.num_vertices, g_s);
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
	
	private static void count_triangles(int graph_enum, String[] all_m, double[] all_e) {
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

	public static void out_graph_statistics() {
		int[] all_graphs 	= Config.graphs;
				
		for(int grah_enum : all_graphs) {
			Graph g = DataLoader.get_graph(grah_enum);
			System.out.println("Graph enum "+grah_enum);
			System.out.println(g.name+"\t"+g.num_vertices+"\t"+g.int_num_edges()+"\t"+((double)g.int_num_edges()/(double)g.num_vertices)+"\t"+Metrics.count_triangles(g));
		}
		
		
	}
	
	public static void experiment_triangle_error() {
		int[] all_mechanism = Config.mechanism;
		int[] all_graphs 	= Config.graphs;
		double[] all_epsilon= Config.all_eps;
		experiment_triangle_error(all_mechanism, all_epsilon, all_graphs);	
	}
	
	private static void experiment_triangle_error(String[] all_mechanism, double[] all_epsilon, int[] all_graphs) {
		for(int graph_id  : all_graphs) {
			System.out.println("Graph id="+graph_id+" "+DataLoader.get_graph_name(graph_id));
			count_triangles(graph_id, all_mechanism, all_epsilon);
		}
	}
	
	private static void experiment_triangle_error(int[] all_mechanism, double[] all_epsilon, int[] all_graphs) {
		experiment_triangle_error(get_mechanism_names(all_mechanism), all_epsilon, all_graphs);
	}

	public static void experiment_local_error() {
		int[] all_mechanism = Config.mechanism;
		int[] all_graphs 	= Config.graphs;
		double[] all_epsilon= Config.all_eps;
		experiment_local_error(all_mechanism, all_epsilon, all_graphs);	
	}
	
	private static void experiment_local_error(int[] all_mechanism, double[] all_epsilon, int[] all_graphs) {
		experiment_local_error(get_mechanism_names(all_mechanism), all_epsilon, all_graphs);
	}

	private static void experiment_local_error(String[] all_mechanism, double[] all_epsilon, int[] all_graphs) {
		for(int graph_id  : all_graphs) {
			System.out.println("Graph id="+graph_id+" "+DataLoader.get_graph_name(graph_id));
			graph_statistics(graph_id, all_mechanism, all_epsilon);
		}
	}
	
	private static String[] get_mechanism_names(int[] all_mechanism){
		String[] names = new String[all_mechanism.length];
		for(int i=0;i<all_mechanism.length;i++) {
			names[i] = name(all_mechanism[i]);
		}
		return names;
	}
	
	public static void matrialize_private_graphs(){
		//int[] graphs = {DataLoader.CONGRESS_TWITTER};
		int[] graphs = Config.graphs;
		//int[] mechanism = {K_EDGE_NON_PRIVATE, K_EDGE_SEQ_RR, K_EDGE_PART_RR, TOP_K, GUESS, RANDOM_RESPONSE};
		//int[] mechanism = {M_PART_2};
		int[] mechanism = Config.mechanism;
		int num_repitions = TestGraphs.num_repitions;
		double[] all_eps = Config.all_eps;
		matrialize_private_graphs(graphs, mechanism, num_repitions, all_eps);
	}
	
	public static void main(String[] args) {
		matrialize_private_graphs();
		//out_graph_statistics();
		experiment_local_error();	
		//experiment_degree_distribution();
		//experiment_triangle_error();
		//int k= 10; experiment_inf_dom(k);
	}

	public static void experiment_inf_dom(int k) {//TODO k
		int[] all_mechanism = Config.mechanism;
		int[] all_graphs 	= Config.graphs;
		double[] all_epsilon= Config.all_eps;
		
		experiment_inf_dom(all_graphs, get_mechanism_names(all_mechanism), all_epsilon);
	}

	private static void experiment_inf_dom(int[] all_graphs, String[] all_mechanism, double[] all_epsilon) {
		for(int graph_enum : all_graphs) {
			experiment_inf_dom(graph_enum, all_mechanism, all_epsilon);
		}
	}

	private static void experiment_inf_dom(int graph_enum, String[] all_m, double[] all_e) {
		System.out.println("***experiment_inf_dom("+DataLoader.get_graph_name(graph_enum)+") " + Arrays.toString(all_m)+" "+Arrays.toString(all_e));
		Graph g = DataLoader.get_graph(graph_enum);
		String path = "./data/inf_dom/"+DataLoader.get_graph_name(graph_enum);
		
		File folder = new File(path);
		if(!folder.exists()) {
			folder.mkdirs();
		}
		
		try (PrintWriter out = new PrintWriter(path+"/org_g.tsv")) {
			String header = "id\t|inf_dom|\tlength_path\tpp"; 
			//System.out.println(header);
			out.println(header);
			InfDom[] inf_doms = ProximityPrestige.run(g);
			for(InfDom i_d : inf_doms) {
				//System.out.println(i_d);
				out.println(i_d.toString());
			}
			Graph.print_writer = null;
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		
		ArrayList<MaterializedGraphs> mg_s = DataLoader.get_sanitized_graphs(graph_enum);
		
		for(String m : all_m) {
			System.out.println(m);
			ArrayList<MaterializedGraphs> graphs_of_mechanism = MaterializedGraphs.filter_by_mechanism(mg_s, m);

			for(double epsilon : all_e) {
				System.out.println(epsilon);
				ArrayList<MaterializedGraphs> temp = MaterializedGraphs.filter_by_epsilon(graphs_of_mechanism, epsilon);
				ArrayList<Graph> g_s = DataLoader.get_sanitized_graphs(temp);
				
				for(int i=0;i<g_s.size();i++) {
					Graph san_g = g_s.get(i);
					try (PrintWriter out = new PrintWriter(path+"/"+m+"-"+epsilon+"-"+i+".tsv")) {
						String header = "id\t|inf_dom|\tlength_path\tpp"; 
						//System.out.println(header);
						out.println(header);
						InfDom[] inf_doms = ProximityPrestige.run(san_g);
						for(InfDom i_d : inf_doms) {
							//System.out.println(i_d);
							out.println(i_d.toString());
						}
						Graph.print_writer = null;
					} catch (FileNotFoundException e) {
						e.printStackTrace();
					}
				}
				
				
			}
		}
	}

	public static void experiment_degree_distribution() {
		int[] all_mechanism = Config.mechanism;
		int[] all_graphs 	= Config.graphs;
		double[] all_epsilon= Config.all_eps;
		experiment_degree_distribution(all_graphs, get_mechanism_names(all_mechanism), all_epsilon);
	}

	private static void experiment_degree_distribution(int[] all_graphs, String[] all_mechanism, double[] all_epsilon) {
		for(int graph_enum : all_graphs) {
			experiment_degree_distribution(graph_enum, all_mechanism, all_epsilon);
		}
	}

	private static void experiment_degree_distribution(int graph_enum, String[] all_m, double[] all_e) {
		System.out.println("***experiment_degree_distribution("+DataLoader.get_graph_name(graph_enum)+") " + Arrays.toString(all_m)+" "+Arrays.toString(all_e));
		Graph g = DataLoader.get_graph(graph_enum);
		final int[] org_graph_degrees = g.get_out_degree_per_vertex();
		Arrays.sort(org_graph_degrees);
		Histogram hist = new Histogram(org_graph_degrees, 20);
		System.out.println(hist);
		
		ArrayList<MaterializedGraphs> mg_s = DataLoader.get_sanitized_graphs(graph_enum);
		ArrayList<ArrayList<Histogram>> all_results = new ArrayList<ArrayList<Histogram>>();
		for(String approach : all_m) {
			System.out.println(approach);
			ArrayList<MaterializedGraphs> graphs_of_mechanism = MaterializedGraphs.filter_by_mechanism(mg_s, approach);
			ArrayList<Histogram> out = new ArrayList<Histogram>();
			for(double epsilon : all_e) {
				System.out.println(epsilon);
				ArrayList<MaterializedGraphs> temp = MaterializedGraphs.filter_by_epsilon(graphs_of_mechanism, epsilon);
				ArrayList<Graph> g_s = DataLoader.get_sanitized_graphs(temp);
				
				Histogram hist_san = new Histogram(hist);
				for(Graph san_g : g_s) {
					hist_san.add_all(san_g.get_out_degree_per_vertex());
				}
				out.add(hist_san);
				
			}
			System.out.println(approach);
			for(Histogram h : out) {
				System.out.println(h);
			}
			all_results.add(out);
		}
		/*String header = "eps\tOriginal\t";
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
		}*/
	}
}
