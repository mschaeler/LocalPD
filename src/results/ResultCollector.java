package results;

import java.util.ArrayList;

import graphs.Graph;
import graphs.Metrics;

public class ResultCollector {
	static ArrayList<double[]> all_metrics = new ArrayList<double[]>();
	static ArrayList<Double> all_page_rank = new ArrayList<Double>();
	static ArrayList<double[]> all_proximity_prestige = new ArrayList<double[]>();
	
	
	public static void init() {
		all_metrics = new ArrayList<double[]>();
		all_page_rank = new ArrayList<Double>();
		all_proximity_prestige = new ArrayList<double[]>();
	}
	public static void collect(Graph g, ArrayList<Graph> all_san_g){
		if(Config.USE_RESULT_STATISTICS) {
			double[] metrics = Metrics.statistics(g, all_san_g);
			all_metrics.add(metrics);	
		}
		if(Config.USE_RESULT_PAGE_RANK) {
			double avg_centralities = Metrics.avg(Metrics.page_rank(g, all_san_g));
			all_page_rank.add(avg_centralities);	
		}
		if(Config.USE_RESULT_PROXIMITY_PRESTIGE) {
			double[] metrics = Metrics.proximity_prestige(g, all_san_g);
			all_metrics.add(metrics);	
		}
		
		
	}
	public static void store(String name, double[] all_eps){
		if(Config.USE_RESULT_STATISTICS) {
			new Statistics(name, all_eps, all_metrics);	
		}
		if(Config.USE_RESULT_PAGE_RANK) {
			new PageRankResult(name, all_eps, all_page_rank);	
		}
		if(Config.USE_RESULT_PROXIMITY_PRESTIGE) {
			new ProximityPrestigeResult(name, all_eps, all_metrics);
		}
	}
}
