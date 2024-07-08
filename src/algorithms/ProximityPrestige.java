package algorithms;

import java.util.ArrayList;
import java.util.HashSet;

import graphs.Graph;

public class ProximityPrestige {
	Graph g;
	ArrayList<Integer>[] incoming_vertices;
	InfDom[] influence_domains;
	
	
	public ProximityPrestige(Graph g){
		this.g = g;
		this.incoming_vertices = g.get_incoming_vertices();
		this.influence_domains = new InfDom[g.num_vertices];  
	}
	
	void get_influence_domain(){
		double start = System.currentTimeMillis();
		System.out.print("get_influence_domain() ");
		for(int node=0;node<g.num_vertices;node++) {
			influence_domains[node] = new InfDom(node, get_influence_domain(node));
		}
		System.out.println("[Done] in "+(System.currentTimeMillis()-start)+" ms");
	}
	
	ArrayList<HashSet<Integer>> get_influence_domain(int node){
		ArrayList<HashSet<Integer>> influence_domain = new ArrayList<HashSet<Integer>>();
		boolean[] reachible = new boolean[g.num_vertices];
		
		reachible[node] = true;
		
		influence_domain.add(new HashSet<Integer>());
		for(int my_direct_neighbor : incoming_vertices[node]) {
			if(!reachible[my_direct_neighbor]) {
				influence_domain.get(0).add(my_direct_neighbor);
				reachible[my_direct_neighbor] = true;
			}
		}
		
		while(true) {
			HashSet<Integer> last_round_new_reachible_nodes = influence_domain.get(influence_domain.size()-1);
			if(last_round_new_reachible_nodes.isEmpty()) {
				return influence_domain;// may be empty directly after start for isolated nodes
			}else{
				HashSet<Integer> new_reachibles = new HashSet<Integer>();
				for(int start_node : last_round_new_reachible_nodes) {
					for(int my_direct_neighbor : incoming_vertices[start_node]) {
						if(!reachible[my_direct_neighbor]) {
							new_reachibles.add(my_direct_neighbor);
							reachible[my_direct_neighbor] = true;
						}
					}
				}
				if(new_reachibles.isEmpty()) {
					return influence_domain;
				}else{
					influence_domain.add(new_reachibles);	
				}
			}
		}
	}
	
	/**
	 * 
	 * @param g
	 * @return array of page rank value per node. Sum over all page rank value = 1.
	 */
	public static double[][] run(Graph g) {
		InfDom[] temp = new ProximityPrestige(g).run();
		double[][] pp_raw = new double[temp.length][3];
		for(int node=0;node<g.num_vertices;node++) {
			pp_raw[node][0] = temp[node].size_influence_domain;
			pp_raw[node][1] = temp[node].average_length_of_shortest_path;
			pp_raw[node][2] = temp[node].proximity_prestige;
		}
		return pp_raw;
	}
	
	InfDom[] run() {
		get_influence_domain();
		
		InfDom[] result = new InfDom[g.num_vertices];
		for(int node=0;node<g.num_vertices;node++) {
			result[node] = influence_domains[node];//TODO why keeping it like this?
		}
		return result;
	}

	public class InfDom{
		int node;
		ArrayList<HashSet<Integer>> inf_dom;
		
		final int size_influence_domain;
		final double average_length_of_shortest_path;
		final double proximity_prestige;
		
		public InfDom(int node, ArrayList<HashSet<Integer>> inf_dom) {
			this.node = node;
			this.inf_dom = inf_dom;
			this.size_influence_domain = size_influence_domain1();
			this.average_length_of_shortest_path = average_length_of_shortest_path1();
			this.proximity_prestige = get_proximity_prestige1();
		}
		private int size_influence_domain1(){
			int size = 0;
			for(HashSet<Integer> i : inf_dom) {
				size +=i.size();
			}
			return size;
		}
		private double average_length_of_shortest_path1() {
			double size_inf_dom = size_influence_domain;
			if(size_inf_dom==0) {
				return 0;
			}
			//We exploit that the incoming nodes are grouped by distance
			double distance_group = 1;
			double sum_distance   = 0;
			for(HashSet<Integer> group : inf_dom) {
				sum_distance += group.size()*distance_group;
				distance_group++;
			}
			return sum_distance / size_inf_dom;
		}
		private double get_proximity_prestige1() {
			double size_i = size_influence_domain;
			if(size_i>0) {
				double nominator = size_i / (double)(g.num_vertices-1);
				double denominator = average_length_of_shortest_path; 
				double pp = nominator / denominator;	
				if(nominator>denominator) {
					System.err.println("nominator>denominator");
				}
				return pp;
			}//else remain zero
			return 0.0d;
		}
	}
}
