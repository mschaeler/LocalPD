package algorithms;

import java.util.ArrayList;
import java.util.Random;

public class UnaryEncoding {
	/**
	 * The classes we survey
	 */
	String[] domain;
	
	/**
	 * Privacy parameters determining the epsilon: p+q=1 must hold;
	 */
	double p,q;
	
	long seed = 1234565;
	Random rand = new Random(seed);
	
	
	public UnaryEncoding(String[] domain) {
		this(domain, 0.75, 0.25);
	}
	
	public UnaryEncoding(String[] domain, final double p, final double q) {
		this.p = p;
		this.q = q;
		this.domain = domain;
	}
	
	boolean[] encode(String my_class){
		boolean[] encoded_response = new boolean[domain.length]; 
		for(int i=0;i<domain.length;i++) {
			String s = domain[i];
			if(s.equals(my_class)){
				encoded_response[i] = true;
				return encoded_response;
			}
		}
		System.err.println("Non of the domain matched");
		return encoded_response;
	}
	
	boolean perturb(final boolean bit) {
		double sample = rand.nextDouble();
		if (bit) {
	        if (sample <= p) return true;
	        else return false;
		}else //bit == false
	        if (sample <= q) return true;
	        else return false;
	}
	
	boolean[] perturb(final boolean[] encoded_response) {
		for(int i=0;i<encoded_response.length;i++) {
			boolean b = perturb(encoded_response[i]);
			encoded_response[i] = b;
		}
		return encoded_response;
	}
	
	double get_epsilon() {
		double eps = Math.log((p*(1-q)) / ((1-p)*q));
		return eps;
	}
	
	double[] aggregate(ArrayList<boolean[]> responses) {
		double[] sums = sum(responses);
		int n = responses.size();
		for(int r=0;r<sums.length;r++) {
			double v = sums[r];
			double sum = (v - n*q) / (p-q);
			sums[r] = sum;
		}
		return sums;
	}

	double[] sum(ArrayList<boolean[]> responses) {
		double[] sums = new double[responses.get(0).length];
		
		for(boolean[] encoded_response : responses) {
			for(int i=0;i<encoded_response.length;i++) {
				if(encoded_response[i]) {
					sums[i]++;
				}
			}
		}
		
		return sums;
	}
}
