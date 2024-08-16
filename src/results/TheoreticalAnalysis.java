package results;

import graphs.Mechanism;

public class TheoreticalAnalysis {
	static double[] all_epsilon = {1,2,3,4,5,6,7,8,9,10};
	static final int thousand = 1000;
	
	public static double[] error_rr(double neighbor_count, double epsilon, double num_vertices) {
		double p = Mechanism.epsilon_to_p(epsilon);
		double deleted_edges = neighbor_count * p;
		double new_fake_edges = (num_vertices - neighbor_count) * p;
		double mae = deleted_edges + new_fake_edges;
		double[] result = {mae, deleted_edges, new_fake_edges};
		return result;
	}
	static void error_rr(){
		System.out.println("Expected error Randomized Response");
		System.out.println("eps\tc\t|V|\t\tmae\tdeleted\tfake edges");
		double[] prefix = new double[3];
		for(double epsilon : all_epsilon){
			prefix[0] = epsilon;
			for(int neighbor_count = 1;neighbor_count<1024;neighbor_count*=2) {
				prefix[1] = neighbor_count;
				for(int num_vertices=10*thousand;num_vertices<10*thousand*thousand;num_vertices*=10) {
					prefix[2] = num_vertices;
					double[] expected_error = error_rr(neighbor_count, epsilon, num_vertices);
					System.out.println(Results.to_tsv(prefix)+"\t"+Results.to_tsv(expected_error));
				}
			}
		}
	}
	
	static void error_k_edge_non_private(){
		System.out.println("Expected error k-edge non private grouping");
		System.out.println("eps\tc\t\tmae\tq_1\tq_2");
		double[] prefix = new double[2];
		for(double epsilon : all_epsilon){
			prefix[0] = epsilon;
			for(int neighbor_count = 1;neighbor_count<1024;neighbor_count*=2) {
				prefix[1] = neighbor_count;
				double[] expected_error = error_k_edge_non_private(neighbor_count, epsilon);
				System.out.println(Results.to_tsv(prefix)+"\t"+Results.to_tsv(expected_error));
			}
		}
	}
	
	public static double[] error_k_edge_non_private(double neighbor_count, double epsilon) {
		double e_1 = 1;
		double e_2 = epsilon-e_1;
		double error_q1 = 1/e_1;
		double error_q2 = neighbor_count / (1.0d+Math.pow(Math.E, e_2));
		double error = error_q1 + error_q2;
		double[] result = {error, error_q1, error_q2};
		return result;
	}
	
	static void error_k_edge_seq_comp(){
		System.out.println("Expected error k-edge seq composition");
		System.out.println("eps\tc\t\tmae\tq_1\tq_2");
		double[] prefix = new double[2];
		for(double epsilon : all_epsilon){
			prefix[0] = epsilon;
			for(int neighbor_count = 1;neighbor_count<1024;neighbor_count*=2) {
				prefix[1] = neighbor_count;
				double[] expected_error = error_k_edge_seq_comp(neighbor_count, epsilon);
				System.out.println(Results.to_tsv(prefix)+"\t"+Results.to_tsv(expected_error));
			}
		}
	}
	
	public static double[] error_k_edge_seq_comp(double neighbor_count, double epsilon) {
		double e_1 = 1;
		double e_2 = epsilon-e_1;
		double error_q1 = 1/e_1;
		double error_q2 = neighbor_count / (1.0d+Math.pow(Math.E, e_2/neighbor_count));
		double error = error_q1 + error_q2;
		double[] result = {error, error_q1, error_q2};
		return result;
	}
	
	static void error_k_edge_random_part(){
		System.out.println("Expected error k-edge random partitionig");
		System.out.println("eps\tc\t\tmae\tq_1\tq_2");
		double[] prefix = new double[2];
		for(double epsilon : all_epsilon){
			prefix[0] = epsilon;
			for(int neighbor_count = 1;neighbor_count<1024;neighbor_count*=2) {
				prefix[1] = neighbor_count;
				double[] expected_error = error_k_edge_random_part(neighbor_count, epsilon);
				System.out.println(Results.to_tsv(prefix)+"\t"+Results.to_tsv(expected_error));
			}
		}
	}
	
	public static double[] error_k_edge_random_part(double neighbor_count, double epsilon) {
		double e_1 = 1;
		double e_2 = epsilon-e_1;
		double error_q1 = 1/e_1;
		double error_q2 = 0;
		double p = Mechanism.epsilon_to_p(e_2);
		
		long n=(long) (neighbor_count-1);
		double[] possibilities = new double[(int) neighbor_count];
		
		double nCk = 1.0d;
		for (int k = 0; k < possibilities.length; k++) {
            //System.out.print(nCk + " ");
            possibilities[k] = nCk;
            nCk = nCk * (n-k) / (k+1);
        }
		//System.out.println();
		double all_cases = Math.pow(2, n);
		
		for(int k=0;k<possibilities.length;k++) {
			int i = k+1;
			double error = (1.0d-p)*i;
			error+=(neighbor_count-i);
			error *= possibilities[k];
			error /= all_cases;
			error_q2 += error;
		}
		
		 
		double error = error_q1 + error_q2;
		double[] result = {error, error_q1, error_q2};
		return result;
	}
	
	/**
	 * 
	 * @param epsilon_j
	 * @param num_vertices |V|
	 * @param degree estimation of the degree c_s
	 * @return
	 * 
	 * solve [//math:ln(((1-a)+(1-(a/b)^k))/(1-(a/b)^k)) = x//] for [//math:x//]
	 * solve [//math:ln(((1-p)+(1-(p/v)^k))/(1-(p/v)^k)) = eps//] for [//math:eps//]
	 * solve [//math:ln(((1-(1/(1+exp(a))))+(1-((1+exp(a))/v)^k))/(1-((1+exp(a))/v)^k)) = eps//] for [//math:eps//]
	 * solve [//math: r = (p (1-(p/v))^k-1) / ((1-(p/v))^k-1)//] for [//math:r//]
	 * solve [//math: r = ((1/(1+exp(eps))) (1-((1/(1+exp(eps)))/v))^k-1) / ((1-((1/(1+exp(eps)))/v))^k-1)//] and eps = 1 and v=1 and k=2 for [//math:r//]
	 */
	static double to_delta_pr(final double epsilon_j, final double num_vertices, final double degree) {
		if(degree>num_vertices) {
			System.err.println(degree>num_vertices);
			return -1;
		}
		if(epsilon_j<0.0d) {
			System.err.println("epsilon_j<0.0d");
			return -1;
		}
		double p = Mechanism.epsilon_to_p(epsilon_j);
		
		double p_by_v = p / num_vertices;
		/**
		 * 
		 */
		double p_by_uniform_sampling = 1.0d - (Math.pow(1.0d-p_by_v, degree));
		double delta_p = ((1.0d-p)+p_by_uniform_sampling)/p_by_uniform_sampling;
		System.out.println("e_j ="+epsilon_j+" |V|="+num_vertices+" c_s="+degree+" p="+p+" p_by_v="+p_by_v+" p_by_uniform_sampling="+p_by_uniform_sampling+" delta_p="+delta_p);
		return delta_p;
	}
	
	
	public static void main(String[] args) {
		//error_rr();
		//error_k_edge_non_private();
		//error_k_edge_seq_comp();
		//error_k_edge_random_part();
		to_delta_pr(1, 1, 1);
		to_delta_pr(1, 2, 1);
		to_delta_pr(1, 2, 2);
	}
}

