package algorithms;

import java.util.ArrayList;
import java.util.HashSet;

import graphs.Graph;

public class InfluenceDomain {
	Graph g;
	ArrayList<Integer>[] neighbors;
	InfDom[] influence_domains;
	
	
	public InfluenceDomain(Graph g){
		this.g = g;
		this.neighbors = g.get_neighbors();
		this.influence_domains = new InfDom[g.num_vertices];  
	}
	
	void get_influence_domain(){
		for(int node=0;node<g.num_vertices;node++) {
			influence_domains[node] = new InfDom(node, get_influence_domain(node));
		}
	}
	
	ArrayList<HashSet<Integer>> get_influence_domain(int node){
		ArrayList<HashSet<Integer>> influence_domain = new ArrayList<HashSet<Integer>>();
		HashSet<Integer> reachible = new HashSet<Integer>();
		
		reachible.add(node);
		
		influence_domain.add(new HashSet<Integer>());
		for(int my_direct_neighbor : neighbors[node]) {
			if(!reachible.contains(my_direct_neighbor)) {
				influence_domain.get(0).add(my_direct_neighbor);
				reachible.add(my_direct_neighbor);
			}
		}
		
		while(true) {
			HashSet<Integer> last_round_new_reachible_nodes = influence_domain.get(influence_domain.size()-1);
			if(last_round_new_reachible_nodes.isEmpty()) {
				return influence_domain;// may be empty directly after start for isolated nodes
			}else{
				HashSet<Integer> new_reachibles = new HashSet<Integer>();
				for(int start_node : last_round_new_reachible_nodes) {
					for(int my_direct_neighbor : neighbors[start_node]) {
						if(!reachible.contains(my_direct_neighbor)) {
							new_reachibles.add(my_direct_neighbor);
							reachible.add(my_direct_neighbor);
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
	public static double[] run(Graph g) {
		return new InfluenceDomain(g).run();
	}
	
	double[] run() {
		// TODO Auto-generated method stub
		return null;
	}

	class InfDom{
		int node;
		ArrayList<HashSet<Integer>> inf_dom;
		public InfDom(int node, ArrayList<HashSet<Integer>> inf_dom) {
			this.node = node;
			this.inf_dom = inf_dom;
		}
	}
}
