/**
Boosting HMM
main algorithm: HMM solution 3
Train a model: Given O, N, and M, find lambda that maximizes probability of O
iteratively adjust lambda = (A, B, pi) to better fit the given observations of O
need to train an HMM per each of the three malware families
using samples from two families on one model to analyze the accuracy of the model 
boosting strategy; randomly assign the initial conditons
*/
import java.io.File;
import java.util.*;
import java.util.ArrayList;
import java.io.FileNotFoundException;
import java.io.IOException;
import javafx.util.Pair;

public class BoostedHMM {
	int N; // rows
	int M; // cols
	int[] observation; //O observation sequence
	int T;
	double [][] a;//state transition probability NxN
	double [][] b; //observation probability matrix; NxM
	double [] Pi;
	double [][] alpha;
	double [][] beta;
	double [][] gamma;
	double [][][] gammaIJ; // used this as gamma[t][i][j]
	double [] c; // scalars at any time t
	double maxStateVal = 0;
	int maxIters = 1000; // max number of re-estiamtion steps
	int iters = 0;

	/**
	* construct a model lambda with matrices A,B and pi
	*/
	public BoostedHMM(int N, int M) {
		this.N = N;
		this.M = M;

		this.a = new double[N][N];
		this.b = new double[N][M];
		this.Pi = new double[N];

		Random generator = new Random(0);

		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				a[i][j] = 1.0 / N * (generator.nextDouble() / 5 + 0.9);
			}
			for (int j = 0; j < M; j++) {
				b[i][j] = 1.0 / M * (generator.nextDouble() / 5 + 0.9);
			}
			Pi[i] = 1.0 / N  * (generator.nextDouble() / 5 + 0.9);
		}
	}

	/**
	* Alpha pass for implementing the forward algorithm with scaling
	* alpha[t][i]: measures the relevant probability up to time t
	*/
	public double[][] alpha_pass() {
		System.out.println("Running alpha pass ...");
		//double old_c0 = c[0];
		c[0] = 0;

		// alpha at time t = 0
		for (int i = 0; i <= N-1; i++) {
			alpha[0][i] = Pi[i]*b[i][observation[0]];
			// scale alpha[i] at time 0
			c[0] += alpha[0][i];
		}
		
		//scale the alpha fo i at time 0
		c[0] = 1/c[0];
		
		for (int i = 0; i <= N-1; i++) {
			alpha[0][i] *= c[0];
		}

		// compute alpha[t][i] at the rest time after time 0
		for (int t = 1; t <= T-1; t++) {
			// c at time t initialize as 0
			c[t] = 0;
			
			// updates each probability of alpha after t=0 based on the pervious state's probability 
			for (int i = 0; i <= N-1; i++ ) {
				alpha[t][i] = 0;
				for (int j = 0; j <= N-1; j++) {
					alpha[t][i] += alpha[t-1][j] * a[j][i];
				}
				// updatesd alpha for one t
				alpha[t][i] *= b[i][observation[t]];
				c[t] += alpha[t][i];
			}

			c[t] = 1/c[t];

			for (int i = 0; i <= N-1; i++) {
				alpha[t][i] *= c[t];
			}
		}
		return alpha;
	}


	/**
	* Beta pass for implementing the forward algorithm with scaling
	* beta[t][i]: measures the relevant probability after time t
	*/
	public double[][] beta_pass() {
		System.out.println("Running beta pass ...");
		
		// let beta of i at time t-1 scaled by scalar at time t-1
		for (int i = 0; i <= N-1; i++) {
			beta[T-1][i] = 1 * c[T-1];
		}

		for (int t = T-2; t >= 0; t --) {
			// updates each probabilty of beta matrix based on the probability after current state
			for (int i = 0; i <= N-1; i++) {
				beta[t][i] = 0;
				for (int j = 0; j <= N-1; j++) {
					beta[t][i] += a[i][j] * b[j][observation[t+1]] * beta[t+1][j];
				}
				
				// scale beta[t][i]
				beta[t][i] = c[t] * beta[t][i];
			}
		}
		return beta;
	}

	/**
	* Gammas using scaled aplhas and betas; solution 3
	* Most likely state at time t: 
	* at time t, gamma[i][j] = (alpha[t][i]*aij*b[j][O[t+1]]beta[t+1][j])/p(O|lambda model)
	* p(O|lambda model) = sum of alpha[T-1][i] form i to N-1
	* at timt t, gamma[i] = sum of gamma[i][j]
	*/
	public double[][] gammas() {
		System.out.println("Running gamma pass ...");
		for (int t = 0; t <= T-2; t++) {
			double denom = 0;
			// rememebr previous gamma[i][j] at time t, this will be used to compute gamma[t][i]
			for (int i = 0; i <= N-1; i++) {
				for (int j = 0; j <= N-1; j++) {
					denom += alpha[t][i]* a[i][j] * b[j][observation[t+1]] * beta[t+1][j];
				}
			}
			for (int i = 0; i <= N-1; i++) {
				gamma[t][i] = 0;
				for (int j = 0; j <= N-1; j++) {
					gammaIJ[t][i][j] = (alpha[t][i] * a[i][j] * b[j][observation[t+1]] * beta[t+1][j])/denom;
					gamma[t][i] += gammaIJ[t][i][j];
				}
			}
		}
		// special case for gamma[T-1][i]
		double denom = 0;
		for (int i = 0; i <= N-1; i++) {
			denom += alpha[T-1][i];
		}
		for (int i = 0; i <= N-1; i++) {
			gamma[T-1][i] = alpha[T-1][i]/denom;
		}
        
		return gamma;
	}

	/**
	* Reestimate A, B, and pi
	*/
	public void reEstimation() {
		System.out.println("Running re-estimate pass ...");
		// re-estimate pi
        for (int i = 0; i <= N-1; i++) {
        	Pi[i] = gamma[0][i];
        }
        
        // re-estimate A
        for (int i = 0; i <= N-1; i++) {
        	for (int j = 0; j <= N-1; j++) {
        		double numer = 0;
        		double denom = 0;
        		for (int t = 0; t <= T-2; t++) {
        			numer += gammaIJ[t][i][j];
        			denom += gamma[t][i];
        		}
        		a[i][j] = numer/denom;
        	}
        }

        // re-estimate B
        for (int i = 0; i <= N-1; i++) {
        	for (int j = 0; j <= M-1; j++) {
        		double numer = 0;
        		double denom = 0;
        		for (int t = 0; t <= T-2; t++) {
        			if (observation[t] == j) {
        				numer += gamma[t][i];
        			}
        			denom += gamma[t][i];
        		}
        		b[i][j] = numer/denom;
        	}
        }
	}

	/**
	* Run all the passes
	*/
	public void run() {
		alpha_pass();
    	beta_pass();
	    gammas();
	    reEstimation();		
	}


	/**
	* Run one iteration for all the training samples-files
	*/
	public void runOneIteration2(ArrayList<ArrayList<Integer>> obsList) {
		int totalLength = 0;
		for (int i = 0; i < obsList.size(); i++) {
			totalLength += obsList.get(i).size();
		}		
		int[] obs = new int[totalLength];
		int idx = 0;
		for (int i = 0; i < obsList.size(); i++) {
			for (int j = 0; j < obsList.get(i).size(); j++) {
				obs[idx] = obsList.get(i).get(j);
				idx += 1;
			}
		}
		System.out.println("Total obs seq length: " + obs.length);

		// set model data
		this.observation = obs;
		this.T = this.observation.length;
		this.alpha = new double[this.T][this.N];
	    this.beta = new double[this.T][this.N];
	    this.gamma = new double[this.T][this.N];
	    this.gammaIJ = new double[this.T][this.N][this.N]; // used this for gamma[t][i][j]
	    this.c = new double[this.T]; // scalars at any time t	

	    run();
		//showResults();
	}

	/**
	* Train model with multiple iterations with stopping criteria
	*/
	public void train(ArrayList<ArrayList<Integer>> obsList) {
		double oldLogProb = Double.NEGATIVE_INFINITY;

		//showResults();
		while (iters < maxIters) {
			System.out.println("======================= Iteration " + iters);
			iters += 1;

			runOneIteration2(obsList); // for all observations

			// compute log[p(O|lambda)]
			double logProb = 0;
			for (int i = 0; i <= T-1; i++) {
				logProb += Math.log(c[i]);
			}
			logProb = -logProb;

			System.out.println("Iter " + iters + ", logP(O|lambda) = " + logProb + ", prev-log = " + oldLogProb);

			// stopping criteria
			if (logProb > oldLogProb) {
				oldLogProb = logProb;
				//showResults();

			} else {
				System.out.println("Stopped at Iter " + iters);
				break;
			}
		}
	}

	/**
	* Print out some inermediate results
	* including pi, A, and B matrices
	*/
	public void showResults() {
	    System.out.println("========");
	    System.out.println("A: " + Arrays.deepToString(this.a));
	    System.out.println("Pi: " + Arrays.toString(this.Pi));

	    System.out.println("B Matrix:");
	    for (int i = 0; i < this.M; i++) {
	    	System.out.print(" " + i + "\t");
	    	for (int j = 0; j < this.N; j++) {
	    		System.out.printf(" %20g ", this.b[j][i]);
	    	}
	    	System.out.println();
	    }
	}

	/**
	* Find the maximum of the gammas matrix
	* This was used for solution 2
	*/
	public double maxGammas(double [][] g) {
		double max = Integer.MIN_VALUE;
		for (int t = 0; t <= T-1; t++) {
			for (int i =0; i <= N-1; i++) {
				if (g[t][i] > max) {
					max = g[t][i];
					maxStateVal = max;
				}
			}
		}
		return maxStateVal;
	}


	/**
	* API for getting A matrix
	*/
	public double[][] getA() {
		System.out.println("The A matrix converges to : " + Arrays.deepToString(this.a));
		return this.a;
	}

	/**
	* API for getting B matrix
	*/
	public double[][] getB() {
		System.out.println("The B Matrix converges to: \n-----------------");
	    for (int i = 0; i < this.M; i++) {
	    	System.out.print(" " + i + "\t");
	    	for (int j = 0; j < this.N; j++) {
	    		System.out.printf(" %20g ", this.b[j][i]);
	    	}
	    	System.out.println();
	    }
		return this.b;
	}

	/**
	* API for getting pi matrix
	*/
	public double[] getPi() {
		System.out.println("The Pi converges to: " + Arrays.toString(this.Pi));
		return this.Pi;
	}

	/**
	* Compute scores for testing samples
	*/
	public ArrayList<Double> computeScore(ArrayList<ArrayList<Integer>> obsList) {
		ArrayList<Double> logProbList = new ArrayList<Double>();
		
		for (int i = 0; i < obsList.size(); i++) {
			int[] obs = new int[obsList.get(i).size()];
			for (int j = 0; j < obsList.get(i).size(); j++) {
				obs[j] = obsList.get(i).get(j);
			}
			System.out.println("Total obs seq length: " + obs.length);
			// set model data
		    this.observation = obs;
			this.T = this.observation.length;
			this.alpha = new double[this.T][this.N];
	    	this.beta = new double[this.T][this.N];
	    	this.gamma = new double[this.T][this.N];
	    	this.gammaIJ = new double[this.T][this.N][this.N]; // used this for gamma[t][i][j]
	    	this.c = new double[this.T]; // scalars at any time t
	 
	    	// only run alpha pass
	    	alpha_pass();
	    	// compute log[p(O|lambda)]
			double logProb = 0;
			for (int k = 0; k <= T-1; k++) {
				logProb += Math.log(c[k]);
			}
		    logProb = -logProb;
		    logProbList.add(logProb);
			//System.out.println("Test one sample logP(O|lambda) = " + logProb);
		}
		return logProbList;
	}

	/**
	* Save score and label into a Tuple
	*/
	public ArrayList<Pair<Double, String>> getTuple(ArrayList<Double> logProbList, String label) {
		
		ArrayList<Pair<Double, String>> tupleList = new ArrayList<Pair<Double, String>>();
		for (int i = 0; i < logProbList.size(); i++) {
			Pair<Double, String> scoreLabelTuple = new Pair<Double, String>(logProbList.get(i), label);
			tupleList.add(scoreLabelTuple);
		}
		System.out.println("Tuple list of " + label + ": " + "\n-----------------------\n" + tupleList);

		return tupleList;
	}
 
    /*
    *  Main function
    */
	public static void main(String[] args) throws FileNotFoundException, IOException {
		int N = 2;
		SampleHandling object = new SampleHandling();
		//ArrayList<ArrayList<Integer>> obsList = object.getObservations("./Malicia (.txt Opcodes)/zbot", 0, 900);
		//ArrayList<ArrayList<Integer>> obsList = object.getObservations("./Malicia (.txt Opcodes)/winwebsec", 0, 900);
		ArrayList<ArrayList<Integer>> obsList = object.getObservations("./Malicia (.txt Opcodes)/zeroaccess", 0, 900);
		//ArrayList<ArrayList<Integer>> obsList = object.getObservations("./Malicia (.txt Opcodes)/zeroaccess", 0, 10);

		int M = object.getNumUniqueOpcode();

		System.out.println(" M = " + M + " N = " + N);
		System.out.println("Number of the observations:" + obsList.size());
		
	    BoostedHMM model = new BoostedHMM(N, M);
	    model.train(obsList);
	    model.getA();
	    model.getB();
	    model.getPi();

	    // test samples
		ArrayList<ArrayList<Integer>> testList = object.getTestObservation("./Malicia (.txt Opcodes)/winwebsec", 900, 1000);
		//ArrayList<ArrayList<Integer>> testList = object.getTestObservation("./Malicia (.txt Opcodes)/zbot", 0, 100);
		//ArrayList<ArrayList<Integer>> testList = object.getTestObservation("./Malicia (.txt Opcodes)/zeroaccess", 900, 1001);
		System.out.println("Number of the test observations:" + testList.size());
		model.computeScore(testList);
		//model.getTuple(model.computeScore(testList), "zbot");
		model.getTuple(model.computeScore(testList), "winwebsec");
		//model.getTuple(model.computeScore(testList), "zeroaccess");
	}
}