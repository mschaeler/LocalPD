package results;

import java.util.ArrayList;

public class PageRankResult extends Results{
	double[] all_eps;
	ArrayList<Double> all_metrics;
	
	public PageRankResult(String name) {
		super(name);
	}
	public PageRankResult(String name, double[] all_eps, ArrayList<Double> all_metrics) {
		this(name);
		this.all_eps = all_eps;
		this.all_metrics = all_metrics;
	}

	@Override
	public void get_header() {
		System.out.println("budget\tavg centrality delta");
		
	}

	@Override
	public void out() {
		for(int i=0;i<all_eps.length;i++) {
			System.out.print("eps="+all_eps[i]+"\t");
			System.out.println(all_metrics.get(i));
		}
	}

}
