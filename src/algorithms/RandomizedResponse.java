package algorithms;

import java.util.Random;

/**
 * Randomized Response to ensure Local Differential privacy
 * Stanley L. Warner. Randomized response: a survey technique for eliminating evasive answer bias. 
 * Journal of the American Statistical Association, 60(309):63–69, 1965. PMID: 12261830. 
 * URL: https://www.tandfonline.com/doi/abs/10.1080/01621459.1965.10480775, doi:10.1080/01621459.1965.10480775.
 * @author b1074672
 *
 */
public class RandomizedResponse {
	long seed = 123456;
	Random rand = new Random(seed);
	final int[] true_data;
	
	RandomizedResponse(int[] true_data){
		this.true_data = true_data;
	}
	
	boolean respond(final int id, final int do_you_belong_to_this_class) {
		if(rand.nextInt(2)==0) {
			//respond truthfully
			if(true_data[id]==do_you_belong_to_this_class) {
				return true;
			}else{
				return false;
			}
		}else{
			//answer randomly (second coin flip)
			return rand.nextInt(2)==0;
		}
	}
	
	boolean[] get_responses(final int[] true_data, final int surveyed_class) {
		boolean[] repsonses = new boolean[true_data.length];
		for(int id=0;id<true_data.length;id++) {
			repsonses[id] = respond(id, surveyed_class);
		}
		return repsonses;
	}
	
	int get_count_true_yesses(final boolean[] responses) {
		
		// we expect 1/4 of the responses to be "yes" based entirely on the coin flip
		// these are "fake" yesses
		int fake_yesses = responses.length/4;
		
		// the total number of yesses recorded
		int num_yesses = sum(responses);

		// the number of "real" yesses is the total number of yesses minus the fake yesses
		int true_yesses = num_yesses - fake_yesses;
		
		return true_yesses;
	}

	private int sum(boolean[] responses) {
		int sum = 0;
		for(boolean b : responses) {
			if(b) sum++;
		}
		return sum;
	}
	
	double pct_error(double true_result, double estimated_result){
		double pct = (estimated_result-true_result)/estimated_result;
		return pct;
	}
	
}
