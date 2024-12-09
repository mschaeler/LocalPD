package graphs;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;

import results.Config;

public class MaterializedGraphs implements Comparable<MaterializedGraphs>{
	public final String path;
	public final String graph_name;
	public final String approach;
	public final int repitition;
	public final double epsilon;
	
	public MaterializedGraphs(String path, String graph_name){
		this.path = path;
		this.graph_name = graph_name;
		String[] tokens = path.replace(".edgelist", "").split("-");
		this.approach = tokens[0];
		this.repitition = Integer.parseInt(tokens[2]);
		this.epsilon = Double.parseDouble(tokens[1]);	
	}
	public Graph getGraph(){
		return DataLoader.load(this);
	}
	public String toString(){
		return graph_name+" "+approach+" "+epsilon;
	}
	public static String[] get_all_mechanisms(ArrayList<MaterializedGraphs> mg_s) {
		HashSet<String> set = new HashSet<String>();
		for(MaterializedGraphs mg : mg_s) {
			set.add(mg.approach);
		}
		
		String[] ret = new String[set.size()];
		int i=0;
		for(String s : set) {
			ret[i++] = s;
		}
		Arrays.sort(ret);
		return ret;
	}
	public static double[] get_all_epsilons(ArrayList<MaterializedGraphs> mg_s) {
		HashSet<Double> set = new HashSet<Double>();
		for(MaterializedGraphs mg : mg_s) {
			set.add(mg.epsilon);
		}
		
		double[] ret = new double[set.size()];
		int i=0;
		for(double eps : set) {
			ret[i++] = eps;
		}
		Arrays.sort(ret);
		return ret;
	}
	public static ArrayList<MaterializedGraphs> filter_by_mechanism(ArrayList<MaterializedGraphs> mg_s, String approach) {
		ArrayList<MaterializedGraphs> filtered = new ArrayList<MaterializedGraphs>();
		for(MaterializedGraphs mg : mg_s) {
			if(mg.approach.equals(approach)) {
				filtered.add(mg);
			}
		}
		return filtered;
	}
	public static ArrayList<MaterializedGraphs> filter_by_epsilon(ArrayList<MaterializedGraphs> mg_s, double epsilon) {
		ArrayList<MaterializedGraphs> filtered = new ArrayList<MaterializedGraphs>();
		for(MaterializedGraphs mg : mg_s) {
			if(mg.epsilon==epsilon) {
				filtered.add(mg);
			}
			if(filtered.size()>=Config.num_repitions) {
				return filtered;		
			}
		}
		return filtered;
	}
	@Override
	public int compareTo(MaterializedGraphs arg0) {
		return this.path.compareTo(arg0.path);
	}
}
