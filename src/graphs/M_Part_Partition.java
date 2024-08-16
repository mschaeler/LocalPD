package graphs;

import java.util.ArrayList;
import java.util.Random;

public class M_Part_Partition {
	static final boolean IS_ORDINARY_PARTITION = true;
	static final boolean IS_FILL_UP_PARTITION  = false;
	static final int NO_TRUE_NEIGHBOR = -1;
	
	final boolean type;
	final int partiton;
	final int true_neighbor;
	public ArrayList<Integer> fake_edges;
	
	private M_Part_Partition(int partition, int my_neighbors, boolean type) {
		this.partiton = partition;
		this.true_neighbor = my_neighbors;
		this.type = type;
		this.fake_edges = new ArrayList<Integer>();
	}
	
	public M_Part_Partition(int partition, ArrayList<Integer> my_neighbors) {
		this(partition, my_neighbors.get(partition), IS_ORDINARY_PARTITION);
		
	}
	public M_Part_Partition(int partition) {
		this(partition, NO_TRUE_NEIGHBOR, IS_FILL_UP_PARTITION);
	}
	
	double uniform_sampling_probability(double epsilon) {
		if(type == IS_FILL_UP_PARTITION) {
			return 1.0d/(double)this.fake_edges.size();
		}else{
			double p = Math.pow(Math.E, epsilon)/(((double)this.fake_edges.size())+Math.pow(Math.E, epsilon));
			p = 1.0d-p;
			p = p/(double)this.fake_edges.size();
			return p;
		}
	}
	
	static boolean check_partitions(final int node, M_Part_Partition[] all_partitions, ArrayList<Integer> true_edges_shortened, ArrayList<Integer> all_fake_edges, Graph g, ArrayList<Integer> true_edges) {
		boolean[] found = new boolean[g.num_vertices];
		
		for(int edge : true_edges_shortened) {
			found[edge] = true;
		}
		for(M_Part_Partition part : all_partitions) {
			for(int edge : part.fake_edges) {
				if(!found[edge]) {
					found[edge] = true;		
				}else{
					System.err.println("Duplicate edge "+edge+" for node="+node);
					return false;
				}
			}
		}
		for(int i=0;i<found.length;i++) {
			if(!found[i]) {
				System.err.println("Missing edge "+i+" for node="+node);
				System.out.println("N_c_s"+true_edges_shortened.toString());
				System.out.println("N"+true_edges.toString());
				System.out.println("N_f"+all_fake_edges.toString());
				for(M_Part_Partition part : all_partitions) {
					System.out.println(part.true_neighbor+" "+part.fake_edges.toString());
				}
				return false;
			}
		}
		return true;
	}

	public static double get_p(final M_Part_Partition[] my_partitions, final double epsilon) {
		double p = Double.POSITIVE_INFINITY;
		for(M_Part_Partition part : my_partitions) {
			double my_p = part.get_p(epsilon);
			if(my_p<p) {
				p = my_p;
			}
		}
		if(p==Double.POSITIVE_INFINITY) {//happens if original neighborlist was empty,, but c_s > 0
			p = 0.0d;
		}
		if(p<0.0d || p>1.0d) {
			System.err.println("get_p() p<0.0d || p>1.0d "+p);
		}
		return p;
	}

	private double get_p(final double epsilon) {
		if(type==IS_FILL_UP_PARTITION) {
			return Double.POSITIVE_INFINITY;
		}else{
			double p = Math.pow(Math.E, epsilon)/(((double)this.fake_edges.size())+Math.pow(Math.E, epsilon));
			return p;		
		}
	}

	public int generalized_rr(final double probability, Random rand) {
		if(type==IS_FILL_UP_PARTITION) {
			return this.fake_edges.get(0);//Recap. fake edges are shuffled. Thus, we uniformly draw from them.
		}
		double dice = rand.nextDouble();
		if(dice<=probability) {//Answer truthfully
			return this.true_neighbor;
		}else{
			return this.fake_edges.get(0);//Recap. fake edges are shuffled. Thus, we uniformly draw from them.
		}
	}
	public String toString() {
		String ret = "Part "+partiton+"|N_f|="+fake_edges.size();
		return ret;
	}
}
