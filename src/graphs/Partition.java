package graphs;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;

public class Partition {
	public final ArrayList<Integer> my_neighbors = new ArrayList<Integer>();
	public final ArrayList<Integer> fake_edges = new ArrayList<Integer>();
	
	public Partition(HashSet<Integer> true_edgdes, HashSet<Integer> fake_edges) {
		for(Integer e : true_edgdes) {
			my_neighbors.add(e);
		}
		for(Integer e : fake_edges) {
			this.fake_edges.add(e);
		}
		Collections.sort(this.my_neighbors);
		Collections.sort(this.fake_edges);
	}

	public Partition(int[] true_edgdes, int[] fake_edges) {
		for(Integer e : true_edgdes) {
			my_neighbors.add(e);
		}
		for(Integer e : fake_edges) {
			this.fake_edges.add(e);
		}
		Collections.sort(this.my_neighbors);
		Collections.sort(this.fake_edges);
	}
}

