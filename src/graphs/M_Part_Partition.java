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
	
	
	public static boolean rebalance_partitions(M_Part_Partition[] all_partitions, double epsilon) {
		int size = all_partitions[0].fake_edges.size();
		boolean has_fill_up_partitions = false;
		boolean[] partition_types = new boolean[all_partitions.length];
		int[] sizes = new int[all_partitions.length];
		int index_first_fill_up_partition = -1;
		ArrayList<M_Part_Partition> fill_up_partitions = new ArrayList<M_Part_Partition>();
		
		for(int i=0;i<all_partitions.length;i++) {
			M_Part_Partition part = all_partitions[i];
			sizes[i] = part.fake_edges.size();
			if(part.type==IS_FILL_UP_PARTITION) {
				partition_types[i] = IS_FILL_UP_PARTITION;
				if(!has_fill_up_partitions) {
					index_first_fill_up_partition = i;
				}
				has_fill_up_partitions = true;
				fill_up_partitions.add(part);
			}else{
				partition_types[i] = IS_ORDINARY_PARTITION;
			}
			
			//check equal size assumption
			if(Math.abs(size-part.fake_edges.size())>1) {
				System.err.println("Violating equal partition size assumption");
			}
		}
		if(!has_fill_up_partitions || partition_types[0]==IS_FILL_UP_PARTITION) {
			return false;
		}
		//Step 1 - Are all ordinary partitions of equal size?
		{
			ArrayList<M_Part_Partition> shorten_me = new ArrayList<M_Part_Partition>();
			for(int i=0;i<index_first_fill_up_partition;i++) {
				if(sizes[i]==size) {
					shorten_me.add(all_partitions[i]);
				}else if(sizes[i]>size) {
					System.err.println("WTF");
				}
			}
			//Does it pay of to distribute |shorten_me| fake edges among fill_up_partitions?
			ArrayList<Integer> fake_edges = new ArrayList<>();
			for(M_Part_Partition part : shorten_me) {
				int edge = remove_last(part);
				fake_edges.add(edge);
			}
			double p = get_p(all_partitions, epsilon);
			double p_uniform_sampling_ordinary_partition = (1.0d-p) / (double) all_partitions[0].fake_edges.size();
			M_Part_Partition shortest = get_shortest(fill_up_partitions);
			double num_edges = ((double) fake_edges.size()) / ((double) fill_up_partitions.size());
			num_edges = Math.ceil(num_edges);
			double p_uniform_sampling_fill_up_partition = 1.0d / (num_edges+(double) shortest.fake_edges.size());
			if(p_uniform_sampling_fill_up_partition>=p_uniform_sampling_ordinary_partition) {
				while(!fake_edges.isEmpty()) {
					shortest = get_shortest(fill_up_partitions);
					int edge = remove_last(fake_edges);
					shortest.fake_edges.add(edge);
				}
			}else{
				//roll back and return
				for(int i=0;i<fake_edges.size();i++) {
					all_partitions[i].fake_edges.add(fake_edges.get(i));//Put them back in again
				}
				return false;
			}
		}
		{
			ArrayList<M_Part_Partition> ordinary_partitions = new ArrayList<M_Part_Partition>();
			for(int i=0;i<index_first_fill_up_partition;i++) {
				ordinary_partitions.add(all_partitions[i]);
			}
			while(true) {
				ArrayList<Integer> fake_edges = new ArrayList<>();
				for(M_Part_Partition part : ordinary_partitions) {
					int edge = remove_last(part);
					fake_edges.add(edge);
				}
				double p = get_p(all_partitions, epsilon);
				double p_uniform_sampling_ordinary_partition = (1.0d-p) / (double) all_partitions[0].fake_edges.size();
				M_Part_Partition shortest = get_shortest(fill_up_partitions);
				double num_edges = ((double) fake_edges.size()) / ((double) fill_up_partitions.size());
				num_edges = Math.ceil(num_edges);
				double p_uniform_sampling_fill_up_partition = 1.0d / (num_edges+(double) shortest.fake_edges.size());
				if(p_uniform_sampling_fill_up_partition>=p_uniform_sampling_ordinary_partition) {
					while(!fake_edges.isEmpty()) {
						shortest = get_shortest(fill_up_partitions);
						int edge = remove_last(fake_edges);
						shortest.fake_edges.add(edge);
					}
				}else{
					//roll back and return
					for(int i=0;i<fake_edges.size();i++) {
						all_partitions[i].fake_edges.add(fake_edges.get(i));//Put them back in again
					}
					return true;
				}
			}
		}
	}
	
	private static M_Part_Partition get_shortest(ArrayList<M_Part_Partition> fill_up_partitions) {
		M_Part_Partition shortest = fill_up_partitions.get(0);
		for(M_Part_Partition part : fill_up_partitions) {
			if(part.fake_edges.size()<shortest.fake_edges.size()){
				shortest = part;
			}
		}
		return shortest;
	}

	private static int remove_last(M_Part_Partition part) {
		return part.fake_edges.remove(part.fake_edges.size()-1);
	}
	private static int remove_last(ArrayList<Integer> l) {
		return l.remove(l.size()-1);
	}


	public static boolean check_partitions(final int node, M_Part_Partition[] all_partitions, ArrayList<Integer> true_edges_shortened, ArrayList<Integer> all_fake_edges, Graph g, ArrayList<Integer> true_edges) {
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
		if(p==Double.POSITIVE_INFINITY) {//happens if original neighbor list was empty,, but c_s > 0
			p = 0.0d;
		}
		if(p<0.0d || p>1.0d) {
			System.err.println("get_p() p<0.0d || p>1.0d "+p);
		}
		return p;
	}

	public static double get_p(final double epsilon, final double partition_size) {
		double p = Math.pow(Math.E, epsilon)/(partition_size+Math.pow(Math.E, epsilon));
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
