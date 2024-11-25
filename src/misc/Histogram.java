package misc;

public class Histogram {
	final int num_bins;
	final int[] pop_count;
	final int bin_size;
	
	public Histogram(int[] values, int num_bins){
		this.num_bins = num_bins;
		this.pop_count = new int[num_bins];
		int max_value = misc.Util.max(values);
		this.bin_size = 1+(max_value / num_bins);
		for(int v : values) {
			int my_bin = get_bin(v);
			pop_count[my_bin]++;
		}
	}

	public Histogram(Histogram hist) {
		this.num_bins = hist.num_bins;
		this.pop_count = new int[num_bins];
		this.bin_size = hist.bin_size;
	}

	private int get_bin(int value) {
		int bin = value / bin_size;
		return bin;
	}
	
	public String toString() {
		double num_values = Util.sum(pop_count);
		String s = "";
		for(double count : pop_count) {
			s+=(count/num_values)+"\t";
		}
		return s;
	}

	public void add_all(int[] out_degree_per_vertex) {
		for(int v : out_degree_per_vertex) {
			int my_bin = get_bin(v);
			pop_count[Math.min(my_bin,this.num_bins-1)]++;
		}
	}
}
