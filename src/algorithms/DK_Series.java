package algorithms;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map.Entry;
import java.util.Random;

import graphs.Graph;
import graphs.Graph.Edge;
import graphs.Mechanism;
import misc.Util;
import results.Config;

public class DK_Series {
	static final int NO_MORE_STUBS = -1;
	static Random rand = new Random(12345687);
	/**
	 * Joint degree matrix, needs to be sanitized
	 */
	final int[][] jdm; 
	/**
	 * offsets.get(degree) -> offset of that degree in <code>jdm</code> and <code>degree_histogram</code>
	 */
	final HashMap<Integer, Integer> offsets;
	/**
	 * |V| this must hold
	 */
	final int num_vertices;
	/**
	 * degree_histogram[offsets.get(degree)] -> count of vertices of that degree.
	 * Needs to be sanitized
	 * XXX for now accept it to be not private
	 */
	int[] degree_histogram;
	
	/**
	 * Name of the original graph
	 */
	final String g_name;
	
	/**
	 * 
	 */
	final int[] degrees_sorted;
	
	private DK_Series(int[][] jdm, HashMap<Integer, Integer> offsets, int num_vertices, int[] degree_histogram, String name){
		this.jdm = jdm;
		this.offsets = offsets;
		this.num_vertices = num_vertices;
		this.degree_histogram = degree_histogram;
		this.g_name = name;
		this.degrees_sorted = new int[offsets.size()];
		int i=0;
		for(Entry<Integer, Integer> e : offsets.entrySet()) {
			degrees_sorted[i++] = e.getKey();//key = degree
		}
		Arrays.sort(degrees_sorted);
		
		//sanity checks
		int num_groups = offsets.size();
		if(jdm.length != num_groups) {
			System.err.println("jdm.length != num_groups");
		}
		if(jdm[0].length != num_groups) {
			System.err.println("jdm[0].length != num_groups");
		}
		if(degree_histogram.length != num_groups) {
			System.err.println("degree_histogram.length != num_groups");
		}
		if(Util.sum(degree_histogram) != num_vertices) {//may happen due to sanitizing
			
		}
		//check_jdm(); There is an issue, because we add each edge twice. Since this is done as in the paper, the conditions appear to be not correct
	}
	
	/**
	 * 
	 */
	public void sanitize(final double epsilon) {
		System.out.println("sanitize()");
		Util.out_tsv(this.jdm);
		
		//Global sensitivity is 4 * max(degree) + 1
		//We draw for each pair jdm[k][l] jdm[l][k] the noise only once, this should half the sensitivity
		//We partition the jdm to reduce noise
		System.out.println("Partitioning scheme and respective noise scale lamda");
		for(int i=0;i<degrees_sorted.length;i++) {
			int degree = this.degrees_sorted[i];
			double lambda = (2.0d*(double) degree+1.0d)/epsilon;//Se above
			int offset_degree = offsets.get(degree);
			
			System.out.print("lambda = "+lambda+" for ("+degree+","+degree+") at ["+offset_degree+","+offset_degree+"] ");
			int val_san = Mechanism.sanitize(jdm[offset_degree][offset_degree],lambda);
			jdm[offset_degree][offset_degree] = val_san;
			for(int j=i-1;j>=0;j--) {
				int other_offset = offsets.get(degrees_sorted[j]);
				System.out.print("("+degree+","+degrees_sorted[j]+") at ["+offset_degree+","+other_offset+"] ");
				System.out.print("("+degrees_sorted[j]+","+degree+") at ["+other_offset+","+offset_degree+"] ");
				
				val_san = Mechanism.sanitize(jdm[offset_degree][other_offset],lambda);
				jdm[offset_degree][other_offset] = val_san;
				jdm[other_offset][offset_degree] = val_san;
			}
			System.out.println();
		}
		Util.out_tsv(this.jdm);
	}
	
	public Graph reconstruct(){
		final int num_groups = jdm.length; 
		final int[][] reconstruction_jdm = new int[num_groups][num_groups];
		Graph g = new Graph(num_vertices, g_name+" 2K series");
		boolean[][] adjaency_matrix = new boolean[this.num_vertices][this.num_vertices];
		
		/**
		 * degree class k -> current degree per node in V_k
		 */
		HashMap<Integer, int[]> remaining_degree_histogram = new HashMap<Integer, int[]>();
		for(Entry<Integer, Integer> e : offsets.entrySet()) {
			int degree = e.getKey();
			int offset = e.getValue();
			int num_vertices_of_degree = this.degree_histogram[offset];
			int[] temp = new int[num_vertices_of_degree];
			Arrays.fill(temp, degree);
			remaining_degree_histogram.put(degree, temp);
		}
		
		for(int line=0;line<num_groups;line++) {
			final int k = get_degree(line);
			final int[] remaining_degree_k = remaining_degree_histogram.get(k);
			for(int column=0;column<num_groups;column++) {
				final int j = get_degree(column);
				final int[] remaining_degree_j = remaining_degree_histogram.get(j);
				
				int counter = 0;
				while(reconstruction_jdm[line][column]<jdm[line][column]) {
					int v = random_choice(k, remaining_degree_k);
					if(v==NO_MORE_STUBS) {
						break;
					}
					
					int w = random_choice(j, remaining_degree_j);
					if(w==NO_MORE_STUBS) {
						break;
					}
					
					int v_global = global_vertex_id(k, v);
					int w_global = global_vertex_id(j, w);
					
					if(!adjaency_matrix[v_global][w_global]){
						if(Config.materialize_graph) {
							Graph.write_edge(v_global, w_global);
						}else{
							g.add_edge(v_global, w_global);	
						}
						adjaency_matrix[v_global][w_global] = true;
					}else if(!adjaency_matrix[w_global][v_global]){
						if(Config.materialize_graph) {
							Graph.write_edge(w_global, v_global);
						}else{
							g.add_edge(w_global, v_global);	
						}
						adjaency_matrix[w_global][v_global] = true;
					}else{
						if(counter<100) {
							remaining_degree_k[v]++;
							remaining_degree_j[w]++;
						}else{
							System.err.println("Edge exists already in both directions ignoring this one ("+v_global+","+w_global+")");	
						}
						counter++;
						
					}
					
					reconstruction_jdm[line][column]++;
					reconstruction_jdm[column][line]++;
				}
			}	
			if(line%1000==0) {
				System.out.print("node="+line+" ");
			}
		}
		System.out.println();
		return g;
	}
	
	private int global_vertex_id(final int degree, final int local_id) {
		int global_id_offset = 0;
		int i=0;
		while(this.degrees_sorted[i]<degree) {
			int offset = offsets.get(degrees_sorted[i++]);
			global_id_offset+=this.degree_histogram[offset];
		}
		return global_id_offset+local_id;
	}

	private static int random_choice(final int degree, final int[] remaining_degree) {
		final double sum = Util.sum(remaining_degree);
		if(sum==0) {
			System.err.println("No more stubs available");
			return NO_MORE_STUBS;
		}
		final double[] p = new double[remaining_degree.length];
		for(int i=0;i<remaining_degree.length;i++) {
			p[i] = ((double)remaining_degree[i])/sum;
		}
		final int index = Mechanism.random_choice(p, DK_Series.rand);
		if(remaining_degree[index]<=0) {
			System.err.println("WTF");
		}else{
			remaining_degree[index]--;
		}
		return index;
	}

	private int get_degree(int value) {
		for(Entry<Integer, Integer> entry : offsets.entrySet()) {
			if(entry.getValue()==value) {
				return entry.getKey();
			}
		}
		System.err.println("get_degree(int) Did not find key for offset "+value);
		return -1;
	}

	/**
	 * According to https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=8674808
	 * there are three conditions.
	 */
	private void check_jdm() {
		Util.out_tsv(jdm);
		
		// I - do we need self loops? for all k, JDM(k, k) <= |Vk| · (|Vk| - 1)
		for(int i=0;i<jdm.length;i++) {
			if(jdm[i][i]/2>degree_histogram[i]*(degree_histogram[i]-1)) {//XXX we add each edge twice...
				int k = get_degree(i);
				int jdm_val = jdm[i][i];
				int count = degree_histogram[i];
				System.err.println("for all k, JDM(k, k) <= |Vk| · (|Vk| - 1) does not hold at k="+k);
			}
		}
		// II - describes a complete bipartite graph between a pair	of degree groups: for all k,l, k != l, JDM(k,l) <= |Vk|·|Vl
		for(int i=0;i<jdm.length;i++) {
			for(int j=0;j<jdm.length;j++) {
				if(i==j) continue;
				if(jdm[i][j]>degree_histogram[i]*degree_histogram[j]) {
					System.err.println("for all k,l, k != l, JDM(k,l) <= |Vk|·|Vl| does not hold at k="+i);
				}
			}
		}
		//ensures that size of degree groups are integers and gives the number of nodes with certain degree
		for(int i=0;i<jdm.length;i++) {
			int k=get_degree(i);
			int sum = 0;
			for(int j=0;j<jdm.length;j++) {
				sum += jdm[i][j];
			}
			//is an integer
			if(sum % k != 0) {
				System.err.println("sum % k != 0");
			}
			sum /= k;
			if(sum!=degree_histogram[i]) {
				System.err.println("sum!=degree_histogram[k]");
			}
		}
	}

	public static void main(String[] args) {
		DK_Series dks = two_k_series(get_exmaple());
		System.out.println(dks.reconstruct());
		dks.sanitize(1.0);
	}
	
	public static DK_Series two_k_series(Graph g) {
		ArrayList<Integer>[] outgoing = g.get_neighbors();
		ArrayList<Integer>[] incoming = g.get_incoming_vertices();
		
		final int[] degree_of_vertex = new int[g.num_vertices]; 		
				
		for(int v=0;v<g.num_vertices;v++) {
			degree_of_vertex[v]+=outgoing[v].size();
			degree_of_vertex[v]+=incoming[v].size();
		}
		System.out.println(Arrays.toString(degree_of_vertex));
		
		HashSet<Integer> temp = new HashSet<Integer>();
		for(int val : degree_of_vertex) {
			temp.add(val);
		}
		int[] unique_degrees = new int[temp.size()];
		{int i=0;
			for(int val : temp) {
				unique_degrees[i++]=val;
			}
			Arrays.sort(unique_degrees);
		}
		HashMap<Integer, Integer> offsets = new HashMap<Integer, Integer>(unique_degrees.length);
		for(int i=0;i<unique_degrees.length;i++) {
			offsets.put(unique_degrees[i], i);
		}
		
		/**
		 * Joint degree matrix
		 */
		int[][] jdm = new int[unique_degrees.length][unique_degrees.length];
		
		//int[] keys = new int[];
		for(Edge e : g.get_edges()) {
			int k = degree_of_vertex[e.from];
			int l = degree_of_vertex[e.to];
			int k_offset_jdm = offsets.get(k);
			int l_offset_jdm = offsets.get(l);
			jdm[k_offset_jdm][l_offset_jdm]++;
			jdm[l_offset_jdm][k_offset_jdm]++;
		}
		
		int[] degree_histogram = new int[unique_degrees.length];
		for(int d : degree_of_vertex) {
			int offset = offsets.get(d);
			degree_histogram[offset]++;
		}
				
		System.out.println();
		return new DK_Series(jdm, offsets, g.num_vertices, degree_histogram, g.name);
	}
	
	/**
	 * https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=8674808
	 * @return
	 */
	public static Graph get_exmaple() {
		Graph g = new Graph(7, "2_k_series");
		g.add_edge(0,1);
		g.add_edge(0,2);
		g.add_edge(1,3);
		g.add_edge(1,2);
		g.add_edge(2,3);
		g.add_edge(3,4);
		g.add_edge(2,4);
		g.add_edge(3,5);
		g.add_edge(4,6);
		return g;
	}
}
