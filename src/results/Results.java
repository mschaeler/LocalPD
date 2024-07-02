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
		System.out.println("*****Results.all_out() Statistics");
		for(Results r : all_results) {
			if(r instanceof Statistics) {
				System.out.println(r.name);
				r.get_header();
				r.out();
			}
		}
		System.out.println("*****Results.all_out() CentralityResult");
		for(Results r : all_results) {
			if(r instanceof CentralityResult) {
				System.out.println(r.name);
				r.get_header();
				r.out();
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
