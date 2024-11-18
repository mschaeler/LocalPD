package results;

import graphs.DataLoader;
import graphs.TestGraphs;

public class Config {
	/**
	 * Fix for partitioned Chung Lu model: Use fix that does not omit parts of the weigths in case (u,v) belong to different partitions
	 */
	public static final boolean USE_FIX = true;
	public static boolean USE_RESULT_STATISTICS 		= false;
	public static boolean USE_RESULT_PAGE_RANK 		   	= false;
	public static boolean USE_RESULT_PROXIMITY_PRESTIGE = false;
	
	/**
	 * TODO
	 * Defines whether fake edges are returned or no edge at all.
	 */
	public static boolean RETURN_FAKE_EDGE = true;
	
	public static double[] all_eps = {1,2,3,4,5,6,7,8,9,10};
	public static int[] mechanism = {TestGraphs.M_NAIVE, TestGraphs.M_PART_2, TestGraphs.M_SAMPLE_WEIGHTED, TestGraphs.CHUNG_LU, TestGraphs.TWO_K_SERIES};
	//public static int[] mechanism = {TestGraphs.M_NAIVE};
	//public static int[] graphs = {DataLoader.CONGRESS_TWITTER, DataLoader.ADVOGATO, DataLoader.ENRON_SINGLE_EDGE};
	public static int[] graphs = {DataLoader.ENRON_SINGLE_EDGE};
	
	public static boolean materialize_graph = false;
}
