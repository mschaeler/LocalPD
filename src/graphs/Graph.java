package graphs;

import java.util.ArrayList;
import java.util.HashSet;

public class Graph {
	int id_e = 0;
	
	public int num_vertices;
	private final ArrayList<Egde> edges = new ArrayList<Graph.Egde>();
	
	public Graph(int num_vertices) {
		this.num_vertices = num_vertices;//They exist only virtually
	}
	
	public Graph(boolean[][] adjacency_matrix) {
		int num_lines 	= adjacency_matrix.length;
		int num_columns = adjacency_matrix[0].length;
		if(num_lines!=num_columns) {
			System.err.println("Graph(bool[][]): num_lines!=num_columns");
		}
		this.num_vertices = num_lines;
		for(int l=0;l<num_lines;l++) {
			for(int c=0;c<num_columns;c++) {
				if(adjacency_matrix[l][c]) {//there is an edge
					new Egde(l, c);
				}
			}
		}
	}

	public Graph(ArrayList<Integer>[] neighbor_list) {
		this(neighbor_list.length);
		for(int node=0;node<neighbor_list.length;node++) {
			for(int target : neighbor_list[node]) {
				new Egde(node,target);
			}
		}
	}

	class Egde{
		final int id;
		final int from;
		final int to;
		
		Egde(int from, int to){
			this.from = from;
			this.to   = to;
			this.id = get_next_edge_id();
			edges.add(this);
		}
		
		public String toString() {
			return id+" ("+from+"->"+to+")";
		}

		public ArrayList<Integer> to_list() {
			ArrayList<Integer> e = new ArrayList<Integer>(2);
			e.add(from);
			e.add(to);
			return e;
		}
	}

	final int get_next_edge_id() {
		return id_e++;
	}
	
	public int[][] get_edges_for_sanitation(){
		int[][] e = new int[edges.size()][2];
		for(int i=0;i<edges.size();i++) {
			Egde my_edge = edges.get(i);
			if(my_edge.from>this.num_vertices) {
				System.err.println("get_edges_for_sanitation");
			}
			if(my_edge.to>this.num_vertices) {
				System.err.println("get_edges_for_sanitation");
			}
			e[i][0] = my_edge.from;
			e[i][1] = my_edge.to;
		}
		return e;
	}
	public ArrayList<Integer>[] get_neighbors(){
		@SuppressWarnings("unchecked")
		ArrayList<Integer>[] edge_list = new ArrayList[this.num_vertices];
		for(int i=0;i<this.num_vertices;i++) {
			edge_list[i] = new ArrayList<Integer>();
		}
		for(Egde e : edges) {
			if(edge_list[e.from].contains(e.to)) {
				System.err.println("get_neighbors(): Duplicate edge");
			}
			edge_list[e.from].add(e.to);
		}
		return edge_list;
	}
	
	public boolean[][] get_adjancency_matrix_as_bit_vector(){
		boolean[][] neighbors = new boolean[this.num_vertices][this.num_vertices];
		for(Egde e : edges) {
			if(neighbors[e.from][e.to] == true) {
				System.err.println("get_adjancency_matrix_as_bit_vector(): Duplicate edge");
			}
			neighbors[e.from][e.to] = true;
		}
		return neighbors;
	}
	
	public static Graph get_example() {
		Graph g = new Graph(60);
		g.add_edge(1-1,2-1);
		g.add_edge(2-1,4-1);
		g.add_edge(3-1,1-1);
		g.add_edge(3-1,2-1);
		g.add_edge(4-1,3-1);
		return g;
	}

	void add_edge(int from, int to) {
		new Egde(from, to);
	}
	
	public String toString() {
		StringBuffer sb = new StringBuffer();
		sb.append("Graph with "+this.num_vertices+" vertices and "+this.edges.size()+" edges");
		
		
		//String ret = "Graph with "+this.num_vertices+" vertices and "+this.edges.size()+" edges";
		int counter = 0;
		for(Egde e : this.edges) {
			if(counter > 10) break;
			sb.append(" "+e);
			counter++;
			//ret+="\n"+e;
		}
		return sb.toString();
	}
	
	static void stop(double start, String text) {
		System.out.println(text+" [DONE] in "+(System.currentTimeMillis()-start));
	}
	
	public static Graph dedup_edges(Graph input) {
		double start = System.currentTimeMillis();
		System.out.println("dedup_edges(Graph g) Input:");
		System.out.println(input);
		
		Graph output = new Graph(input.num_vertices);//creates graph with not edges
		HashSet<ArrayList<Integer>> unique_edges = new HashSet<ArrayList<Integer>>(input.num_vertices);
		int counter = 0;
		for(Egde e : input.edges) {
			unique_edges.add(e.to_list());
			if(counter++%100000==0) {
				System.out.println(counter);
			}
		}
		
		for(ArrayList<Integer> e_as_list : unique_edges) {
			output.add_edge(e_as_list.get(0), e_as_list.get(1));
		}
		//Collections.sort(output.edges);
		
		stop(start, "dedup_edges(Graph g) Ouput:");
		System.out.println(output);
		return output;
	}
}
