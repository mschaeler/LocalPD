package results;

public class Config {
	/**
	 * Fix for partitioned Chung Lu model: Use fix that does not omit parts of the weigths in case (u,v) belong to differnt partitions
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
}
