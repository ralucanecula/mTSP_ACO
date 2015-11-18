package multiobjectiveACO;

import java.util.ArrayList;
import java.util.Map;

import multiobjectiveACO.Ants.Ant;
import multiobjectiveACO.ParetoFront.Solution;

/* this Eclipse project represents the ACO algorithm adapted from the mTSP_ACO global 
 * project in Eclipse, but in which there are considered more than one objective, such that 
 * a multiobjective ACO algorithm called MoACO/D-ACS (multi-objective ant colony optimization 
 * based on decomposition combined with ACS) - see the article "Multi-objective ant colony optimization 
 * based on decomposition for bi-objective traveling salesman problems", is applied 
 * to solve the mTSP problem and find a good and diversified set of non-dominated solutions 
 * belonging to the Pareto front; the 2 considered objectives that must be optimized simultaneously 
 * are: 1. minimizing the total distance of tours traveled by the salesmen and 2. balancing the
 *  tours traveled by the salesmen, i.e. minimizing the length of the longest tour or minimizing the
 *  amplitude, which is the difference between the length of the longest tour and the shortest tour 
 * 
 * it's an improved ACO algorithm that actually add the depot city in the tour that each salesman travel,
 * and pheromone is deposited on the edges incident (that enters or leave) with the depot city 
 * this approach should be a better one, since in nature this is how things actually happen: the ant 
 * leaves/deposit pheromone anywhere it travels in the search for food
 */

/* contains the main entry point for the implementation of ACO for solving the mTSP problem on mTSP 
 * instances provided by Carter
 */
public class TSP_ACO {
	
	/** File name to be used for input data set */
	private static String fileName = "rat99.tsp";
	//private static String fileName = "mtsp20.txt";
	
	/** the number of salesmen from the mTSP instance */
    private static int m = 7;
    		
    /** the number of objectives to be considered in the multiobjective (bi-criteria) mTSP problem */
    public static int k = 2; 
    
    /** the number of subproblems or subcolonies or ants in colony */
    //public static int N = 30;
    public static int N = 5;
    
    /** the number of ants in each subcolony or how many weight vectors to be considered in the 
     * neighborhood of each weight vector or by how many ants to be solved each subproblem */
    //public static int T = 5;
    public static int T = 3;
    
    /** name of the TSP input instance from TSPLIB **/
    public static String instanceName = new String();
    
    /** for each subproblem keeps weight values for each objective */
    static double weightVectors[][];

	
	//checks whether termination condition is met
    static boolean termination_condition() {   
    	//return (((InOut.n_tours >= InOut.max_tours) && (Timer.elapsed_time() >= InOut.max_time)) || (Ants.best_so_far_ant.tour_length <= InOut.optimal));
    	//return ((InOut.n_tours >= InOut.max_tours) && (Timer.elapsed_time() >= InOut.max_time));
    	return (InOut.iteration >= InOut.max_iterations);
    }
    
    static int[][] divideProblem() {
    	double nr;
    	double[][] distance = new double[N][N];
    	int[][] nearestWeightVectors = new int[N][T]; 
    	
    	//generate N random weight vectors
    	for (int i = 0; i < N; i++) {
    		//nr should be in the range of (0, 1)
    		nr = Math.random();
    		while (nr == 0.0) {
    			nr = Math.random();
    		}
    		weightVectors[i][0] = nr;
    		weightVectors[i][1] = 1 - nr;		
    	}
    	
    	//compute distances between each pair of weight vectors
    	distance = Tsp.computeWeightVectorsDistances(weightVectors);
    	
    	//determine T closest weight vectors for each weight vector
    	nearestWeightVectors = Tsp.computeNearestWeightVectors(distance);
    	
    	return nearestWeightVectors;
    }
    
    //check if there is still an ant with left cities to visit
    static boolean isDone(int indexSubproblem) {
    	boolean done = true;
    	
    	for (int j = 0; j < T; j++) {
    		Ant a = Ants.subcolonies[indexSubproblem][j];	
    		if (a.nextSubproblemId < a.idSolvedSubproblems.size()) {
	    		if (a.toVisit.get(a.nextSubproblemId) > 0) {
	    			return false;
	    		}	
    		}
    		
    	}
    	
    	return done;
    }

    //manage the solution construction phase (construct a set of complete and closed tours, 
    //each tour for one of the m salesmen)
    static void construct_solutions() {
		int step; /* counter of the number of construction steps */
		int salesman;
		
		/* Mark all cities as unvisited */
		for (int j = 0; j < Ants.n_ants; j++) {
		    Ants.ant_empty_memory(Ants.ants[j]);
		}
	
		step = 0;
		
		//for every subproblem/subcolony
		for (int i = 0; i < N; i++) {
			for (int l = 0; l < T; l++) {
				/* Place the ants on same initial city */
				Ant a = Ants.subcolonies[i][l];
				for (int j = 0; j < MTsp.m; j++) {
					//place each ant on the depot city, which is the start city of each tour
					// -1 is a special marker for the deport city, so that it's not be confused with the rest of the cities
					// all the rest of the cities are represented by integer values > 0
					a.tours.get(a.nextSubproblemId)[j].add(-1); 				
				}
				a.nextSubproblemId++;
			}
		}
		
		//reset value of nextSubproblemId field of all the ants
		for (int l = 0; l < Ants.n_ants; l++) {
			Ants.ants[l].nextSubproblemId = 0;
		}
		
		//for every subproblem/subcolony
		for (int i = 0; i < N; i++) {
			while (!isDone(i)) {
				for (int j = 0; j < T; j++) {				
					Ant a = Ants.subcolonies[i][j];
					if (a.nextSubproblemId < a.idSolvedSubproblems.size() && a.toVisit.get(a.nextSubproblemId) > 0) {
						//choose for each ant in a probabilistic way by some type of roullette wheel selection 
						//which salesman to consider next, that will visit a city
			    		salesman = (int)(Math.random() * MTsp.m);
						Ants.neighbour_choose_and_move_to_next(a, salesman, i);
						if (Ants.acs_flag)
						    Ants.local_acs_pheromone_update(a, salesman, i);						
					}
				}
			}
			for (int j = 0; j < T; j++) {	
				Ant a = Ants.subcolonies[i][j];
				a.nextSubproblemId++;
			}
		}

		//reset value of nextSubproblemId field of all ants
		for (int j = 0; j < Ants.n_ants; j++) {
			Ants.ants[j].nextSubproblemId = 0;
		}
		
		//for every subproblem/subcolony
		for (int i = 0; i < N; i++) {
		  for (int l = 0; l < T; l++) {
			  Ant a = Ants.subcolonies[i][l];
			  for (int j = 0; j < MTsp.m; j++) {
				step = a.tours.get(a.nextSubproblemId)[j].size();
				a.tours.get(a.nextSubproblemId)[j].add(step, -1);
				
				a.tour_lengths.get(a.nextSubproblemId)[j] = Tsp.compute_tour_length_(a.tours.get(a.nextSubproblemId)[j]);
				double value = a.total_tour_lengths.get(a.nextSubproblemId);
				a.total_tour_lengths.set(a.nextSubproblemId, value + a.tour_lengths.get(a.nextSubproblemId)[j]);
				
			    if (Ants.acs_flag)
			    	Ants.local_acs_pheromone_update(a, j, i);
			  }
			  //Ants.ants[k].total_tour_length = Tsp.compute_tour_lengths(Ants.ants[k].tours);
			  a.costObjectives.get(a.nextSubproblemId)[0] = a.total_tour_lengths.get(a.nextSubproblemId);
			  a.costObjectives.get(a.nextSubproblemId)[1] = Ants.computeToursAmplitude(a);
			  a.nextSubproblemId++;
		  }
		}
		InOut.n_tours += (Ants.n_ants * MTsp.m); //each ant constructs a complete and closed tour
		
		//reset value of nextSubproblemId field of all ants
		for (int j = 0; j < Ants.n_ants; j++) {
			Ants.ants[j].nextSubproblemId = 0;
		}
		
    }

    //initialize variables appropriately when starting a trial
    static void init_try() {

		Timer.start_timers();
		InOut.time_used = Timer.elapsed_time();
		InOut.time_passed = InOut.time_used;
	
		/* Initialize variables concerning statistics etc. */
		InOut.n_tours = 1 * MTsp.m;
		InOut.iteration = 1;
		//Ants.best_so_far_ant.total_tour_length = Integer.MAX_VALUE;
		InOut.found_best = 0;
		double sum;
		
		/*
		 * Initialize the Pheromone trails, only if ACS is used, Ants.pheromones
		 * have to be initialized differently
		 */
		double[] objectiveValues = Ants.nn_tour();
		if (!(Ants.acs_flag)) {
			for (int i = 0; i < N; i++) {
				sum = 0.0;
				for (int nrObj = 0; nrObj < TSP_ACO.k; nrObj++) {
					sum +=  weightVectors[i][nrObj] * objectiveValues[nrObj];
				}
				Ants.trail_0[i] = (double)T / sum;		
				
				/*
			     * in the original papers on Ant System it is not exactly defined what the
			     * initial value of the Ants.pheromones is. Here we set it to some
			     * small constant, analogously as done in MAX-MIN Ant System.
			     */
			    Ants.init_pheromone_trails(i, Ants.trail_0[i]);			
			}
		}
		
		if (Ants.acs_flag) {
			for (int i = 0; i < N; i++) {
				sum = 0.0;
				for (int nrObj = 0; nrObj < TSP_ACO.k; nrObj++) {
					sum += weightVectors[i][nrObj] * objectiveValues[nrObj];
				}
				Ants.trail_0[i] = 1. / ((double)MTsp.n * sum);

			    Ants.init_pheromone_trails(i, Ants.trail_0[i]);
			}
		}
	
		//this is no longer possible since the total information depend on the weights chosen random
		//at each iteration by each ant, which is a dynamic information
		/* Calculate combined information Ants.pheromone times heuristic information */
		//Ants.compute_total_information();
		
		//there is no sense in doing this since the heuristic value for the second objective(amplitude)
		//is not a static value and its value can be determined only during the tour construction process
		//for each subproblem, initialize the matrix with heuristic information
		/*for (int i = 0; i < N; i++) {
			Ants.init_heuristic_matrix(i);
		}*/
		
		System.out.print("Initial values for z: (");
		//compute the initial value for the reference point
		for (int nrObj = 0; nrObj < TSP_ACO.k; nrObj++) {
			Ants.z[nrObj] = objectiveValues[nrObj]; 
			System.out.print(Ants.z[nrObj] + " ");
		}	
		System.out.print(")" + "\n");
    }
    
    //manage some statistical information about the trial, especially if a new best solution
    //(best-so-far) is found and adjust some parameters if a new best solution is found
    static void update_statistics(int trial, int iterationNr) {
    	int[] idBestAnts = new int[N];
    	
    	/*if (iterationNr == 0 || iterationNr == InOut.max_iterations - 1) {
    		System.out.println("Trail #" + trial +  " iteration " + iterationNr + " size bestSoFarPareto=" + ParetoFront.bestSoFarPareto.size());
    	}*/
    	
    	//the initial pheromone trail is reinitialized by applying the weighted average objective function values 	
    	for (int i = 0; i < N; i++) {
    	  double[] sum1 = new double[TSP_ACO.k];
    	  double[] averageValues = new double[TSP_ACO.k];
  		  for (int k1 = 0; k1 < T; k1++) {
  			  Ant a = Ants.subcolonies[i][k1];
  			  for (int nrObj = 0; nrObj < TSP_ACO.k; nrObj++) {
  				//System.out.println("Colony:" + i + " ant id:" + a.id + " nextSubproblemId:" + a.nextSubproblemId + " total solved problems:" + a.idSolvedSubproblems.size());
  				sum1[nrObj] += a.costObjectives.get(a.nextSubproblemId)[nrObj];
  			  }		  			  
  			  a.nextSubproblemId++;
  		  }
  		  for (int nrObj = 0; nrObj < TSP_ACO.k; nrObj++) {
  			  averageValues[nrObj] = sum1[nrObj] / (double) T;
		  }
  		  double sum2 = 0.0; 
		  for (int nrObj = 0; nrObj < TSP_ACO.k; nrObj++) {
			  sum2 += weightVectors[i][nrObj] * averageValues[nrObj]; 
		  } 
		  Ants.trail_0[i] = 1.0 / sum2;
    	}
    	
    	//reset value of nextSubproblemId field of all the ants
		for (int k1 = 0; k1 < Ants.n_ants; k1++) {
			Ants.ants[k1].nextSubproblemId = 0;
		}
		
    	//update reference point
    	double min[] = new double[TSP_ACO.k];
    	for (int i = 0; i < TSP_ACO.k; i++) {
    		min[i] = Double.MAX_VALUE;
    	}
    	
    	for (int i = 0; i < N; i++) {
    		 for (int k1 = 0; k1 < T; k1++) {
    			 Ant a = Ants.subcolonies[i][k1];
     			 for (int nrObj = 0; nrObj < TSP_ACO.k; nrObj++) {
     				 if (a.costObjectives.get(a.nextSubproblemId)[nrObj] < min[nrObj]) {
     					min[nrObj] = a.costObjectives.get(a.nextSubproblemId)[nrObj];
     				 }
     			 }
     			 a.nextSubproblemId++;
    		 }
    	}
    	for (int nrObj = 0; nrObj < TSP_ACO.k; nrObj++) {
    		Ants.z[nrObj] = min[nrObj];
    	}
    	
    	//reset value of nextSubproblemId field of all the ants
		for (int j = 0; j < Ants.n_ants; j++) {
			Ants.ants[j].nextSubproblemId = 0;
		}
    	
    	//determine best ant for each subproblem/subcolony
    	for (int i = 0; i < N; i++) {
    		idBestAnts[i] = Ants.find_best(i);
    		Ants.copy_from_to(Ants.subcolonies[i][idBestAnts[i]], Ants.best_so_far_ants[i]);	
    	}
    	
    	//reset value of nextSubproblemId field of all the ants
		for (int j = 0; j < Ants.n_ants; j++) {
			Ants.ants[j].nextSubproblemId = 0;
		}
	
		//compute non-dominated set of solutions (iteration non-dominated front)
		ParetoFront.iterationPareto.clear();
		Solution copyAnt;
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < T; j++) {
				Ant a = Ants.subcolonies[i][j];
				copyAnt = Ants.copyAnt(a);
				ParetoFront.paretoUpdateWithSolution(ParetoFront.iterationPareto, copyAnt);
				a.nextSubproblemId++;
			}		
		}
		
		//reset value of nextSubproblemId field of all the ants
		for (int l = 0; l < Ants.n_ants; l++) {
			Ants.ants[l].nextSubproblemId = 0;
		}
		
		//update BestSoFarPareto external set
		ParetoFront.paretoUpdate(ParetoFront.bestSoFarPareto, ParetoFront.iterationPareto);
		
    }

    //occasionally compute some statistics
   //at every 5 iterations save the value of the best cost of the tour given by the best so far ant
    static void search_control_and_statistics(ArrayList<Double> iterTotalCost, boolean saveDetailedOutput)
    {
    	double totalCost;
    	
    	if (saveDetailedOutput) {
    		if ((InOut.iteration % 5) == 0) {
			    //System.out.println("TSP(" + tspIndex + "): best tour length so far " + Ants.best_so_far_ant[tspIndex].tour_length + ", iteration: " + InOut.iteration);
	    		/*totalCost = Ants.best_so_far_ant.total_tour_length;
	    		iterTotalCost.add(totalCost);*/
    		}
    	}
    }

    //manage global Ants.pheromone deposit for Ant System
    static void as_update() {	
		//for every subproblem/subcolony
		for (int i = 0; i < N; i++) {
			Ants.global_update_pheromone(i);		
		}
		
		/*for (j = 0; j < Ants.n_ants; j++)
		    Ants.global_update_pheromone(Ants.ants[j]);*/
    }

    //manage global Ants.pheromone deposit for Ant Colony System
    static void acs_global_update() {
    	//for every subproblem/subcolony
    	for (int i = 0; i < N; i++) {
		    Ants.global_acs_pheromone_update(i);	
    	}
    	
    }
    
    static void acs_share_pheromone_trails() {
    	//for every subproblem/subcolony
    	for (int i = 0; i < N; i++) {
		    Ants.global_acs_share_pheromone_trails(i);	
    	}
    }

    //manage global Ants.pheromone trail update for the ACO algorithms
	static void pheromone_trail_update()
	{
		//Simulate the Ants.pheromone evaporation of all Ants.pheromones; this is not necessary for ACS
		if (Ants.as_flag) {
			/* evaporate all Ants.pheromone trails */
			Ants.evaporation();
		}
	
		/* Next, apply the Ants.pheromone deposit for the various ACO algorithms */
		if (Ants.as_flag)
		    as_update();
		else if (Ants.acs_flag)
		    acs_global_update();
	
	
		/*
		 * Compute combined information Ants.pheromone times heuristic info after
		 * the Ants.pheromone update for all ACO algorithms except ACS; in the ACS case
		 * this is already done in the Ants.pheromone update procedures of ACS
		 */
		/*if (Ants.as_flag) {
			Ants.compute_total_information();
	    }*/
	}
	
	//after an iteration, a pheromone trail share procedure is evoked to implement the 
	//information share of those subproblems solved by common ants
	static void share_pheromone_trails() {
		if (Ants.acs_flag) {
			acs_share_pheromone_trails();
		}
	}

	public static void main(String[] args) {
		long startTime = System.currentTimeMillis();
		for (int trial = 0; trial < 1; trial++) {
	    	//read the data from the input file
			DataReader reader = new DataReader(fileName);
	        //read the data from the file
	        Data d = reader.read();
	        
	        String[] str = fileName.split("\\.");
			//we are dealing with a input file in TSP format 
			if (str[1].equals("tsp")) {
				instanceName = str[0];
			}
	        
	    	int n = d.getRecords().size();  //nr of cities except the depot city
	    	ArrayList<Record> records = d.getRecords();
	    	Record depotCity = d.getDepotCity();
	    	
	        MTsp mtsp = new MTsp(records, depotCity, n, m);
	        mtsp.printDepotCity();
	        //mtsp.printRecords(records);        
	        
	        Ants.subcolonies = new Ant[N][T];
	        weightVectors = new double[N][TSP_ACO.k];
	        Ants.z = new double[TSP_ACO.k];
	        
			//decompose mTSP problem into N scalar optimization subproblems
	        //divide ant colony
	        int[][] closestWeightVectors = divideProblem();
	        //store indexes of subproblems solved by each ant; it corresponds to A(h) set
	        ArrayList<Integer>[] indexSubproblems = (ArrayList<Integer>[])new ArrayList[N];
	        for (int i = 0; i < N; i++) {
	        	indexSubproblems[i] = new ArrayList<Integer>();
	        }
	        for (int i = 0; i < N; i++) {
	        	for (int j = 0; j < T; j++) {
	        		for (int indexAnt = 0; indexAnt < N; indexAnt++) {
	        			if (closestWeightVectors[i][j] == indexAnt) {
	        				indexSubproblems[indexAnt].add(i);
	        			}
	        		}
	        	}
	        }      
	           
	        System.out.println("Indexes of subproblems solved by each ant");
	        for (int i = 0; i < N; i++) {
	        	System.out.print("Ant " + i + " will solve: ");
	        	for (int j = 0; j < indexSubproblems[i].size(); j++) {
	        		System.out.print(indexSubproblems[i].get(j) + " ");
	        	}
	        	System.out.println();
	        }        
	        
			InOut.init_program(args, closestWeightVectors, indexSubproblems);	
			
		    MTsp.instance.nn_list = Tsp.compute_nn_lists();
			
		    Ants.pheromone = new double[N][][];
			Ants.heuristic = new double[N][][];
			
			for (int i = 0; i < N; i++) {
				Ants.pheromone[i] = new double[MTsp.n + 1][MTsp.n + 1];
				Ants.heuristic[i] = new double[MTsp.n + 1][MTsp.n + 1];
			}
			
			init_try(); 
	
			boolean saveIterCosts = false;
			ArrayList<Double> iterTotalCost = null;
			if (saveIterCosts) {
				//for saving detailed cost values at each 5 iteration
				iterTotalCost = new ArrayList<Double>();
			}
			
			InOut.iteration = 0;
		    while (!termination_condition()) {
				construct_solutions();
				update_statistics(trial, InOut.iteration);
				pheromone_trail_update();
				share_pheromone_trails();
				search_control_and_statistics(iterTotalCost, saveIterCosts);
				InOut.iteration++;
		    }
		    
		    //Utilities.setIterTotalCost(iterTotalCost);
		    
		    InOut.exit_try(trial);
		    ////Utilities.writeParetoSet(ParetoFront.bestSoFarPareto, trial);
		    //System.out.println("Reached " + InOut.iteration + " iterations");
		    //Utilities.writeExcel(MTsp.n, MTsp.m, totalLength);
		    //Utilities.writeResultsExcel(trial, saveIterCosts);
		    ////Utilities.writeParetoSolutions(ParetoFront.bestSoFarPareto);
		    ParetoFront.bestSoFarPareto.clear();
		}
		long endTime = System.currentTimeMillis();
		double difference = (endTime - startTime)/1000.0;
		System.out.println("\nElapsed seconds: " + difference);
    }
}
