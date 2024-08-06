package misc;

public class Util {
	public static double sum(final int[] arr) {
		int sum = 0;
		for(int v : arr) {
			sum += v;
		}
		return sum;
	}
	public static double sum(final double[] arr) {
		double sum = 0;
		for(double v : arr) {
			sum += v;
		}
		return sum;
	}
	
	public static void out_tsv(int[][] matrix) {
		if(matrix==null) {
			System.out.println("Null array");
			return;
		}
		for(int[] arr : matrix) {
			out_tsv(arr);
		}
	}
	
	public static void out_tsv(double[][] matrix) {
		if(matrix==null) {
			System.out.println("Null array");
			return;
		}
		for(double[] arr : matrix) {
			out_tsv(arr);
		}
	}
	public static void out_tsv(int[] arr) {
		if(arr==null) {
			System.out.println("Null array");
			return;
		}
		for(int val : arr) {
			System.out.print(val+"\t");
		}
		System.out.println();
	}
	public static void out_tsv(double[] arr) {
		if(arr==null) {
			System.out.println("Null array");
			return;
		}
		for(double val : arr) {
			System.out.print(val+"\t");
		}
		System.out.println();
	}
	public static double sum(double[][] matrix) {
		double sum = 0.0d;
		for(double[] arr : matrix) {
			sum+=sum(arr);
		}
		return sum;
	}
	public static final int max(int[] arr) {
		int max_value = Integer.MIN_VALUE;
		for(int v : arr) {
			if(v>max_value) {
				max_value = v;  
			}
		}
		return max_value;
	}
}
