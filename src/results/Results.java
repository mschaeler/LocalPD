package results;

import java.util.ArrayList;

public abstract class Results {
	static ArrayList<Results> all_results = new ArrayList<Results>();
	public String name;
	
	public Results(String name) {
		this.name = name;
		all_results.add(this);
	}
	
	public static void all_out() {
		if(Config.USE_RESULT_STATISTICS) {
			System.out.println("*****Results.all_out() Statistics");
			for(Results r : all_results) {
				if(r instanceof Statistics) {
					System.out.println(r.name);
					r.get_header();
					r.out();
				}
			}
		}
		if(Config.USE_RESULT_PAGE_RANK) {
			System.out.println("*****Results.all_out() PageRankResult");
			for(Results r : all_results) {
				if(r instanceof PageRankResult) {
					System.out.println(r.name);
					r.get_header();
					r.out();
				}
			}
		}
		if(Config.USE_RESULT_PROXIMITY_PRESTIGE) {
			System.out.println("*****Results.all_out() ProximityPrestigeResult");
			for(Results r : all_results) {
				if(r instanceof ProximityPrestigeResult) {
					System.out.println(r.name);
					r.get_header();
					r.out();
				}
			}
		}
	}
	
	public abstract void get_header();
	public abstract void out();
	
	static String to_tsv(double[] arr) {
		if(arr==null) {
			return "null";
		}
		String s = "";
		for(double d : arr) {
			s+=d+"\t";
		};
		return s;
	}
}
