package results;

import java.util.ArrayList;

public class Statistics extends Results{
	double[] all_eps;
	ArrayList<double[]> all_metrics;
	
	public Statistics(String name) {
		super(name);
	}
	public Statistics(String name, double[] all_eps, ArrayList<double[]> all_metrics) {
		this(name);
		this.all_eps = all_eps;
		this.all_metrics = all_metrics;
	}
	public void get_header() {
		System.out.println("budget\tcorrect\tmissing\tfake edge\tnon edge");
	}
	public void out() {
		for(int i=0;i<all_eps.length;i++) {
			System.out.print("eps="+all_eps[i]+"\t");
			System.out.println(to_tsv(all_metrics.get(i)));
		}
	}
}
