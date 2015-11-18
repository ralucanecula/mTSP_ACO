package multiobjectiveACO;

import java.util.ArrayList;
import java.util.Map;

import multiobjectiveACO.Ants.Ant;

/* this Eclipse project represents the ACO algorithm adapted from the mTSP_ACO global 
 * project in Eclipse, but in which there are considered more than one objective, such that 
 * a multiobjective ACO algorithm called P-ACO (Pareto ant colony optimization) is applied 
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
    
    /** the number of objectives to be considered in the multiobjective (bi-criteria) mTSP problem*/
    public static int k = 2; 
    
    /** name of the TSP input instance from TSPLIB **/
    public static String instanceName = new String();

	
	//checks whether termination condition is met
    static boolean termination_condition() {   
    	//return (((InOut.n_tours >= InOut.max_tours) && (Timer.elapsed_time() >= InOut.max_time)) || (Ants.best_so_far_ant.tour_length <= InOut.optimal));
    	//return ((InOut.n_tours >= InOut.max_tours) && (Timer.elapsed_time() >= InOut.max_time));
    	return (InOut.iteration >= InOut.max_iterations);
    }
    
    /*static int chooseSalesman(double[] probabilities) {
    	int i;
    	double rnd, partialSum = 0.;
    	
    	rnd = Utilities.random01();
    	i = 0;
    	partialSum = probabilities[i];
    	
    	 This loop always stops because last value from probabilities vector is a huge value 
    	while (partialSum <= rnd) {
			i++;
			partialSum += probabilities[i];
	    }
    	return i;
    }*/
    
    //check if there is still an ant with left cities to visit
    static boolean isDone() {
    	boolean done = true;
    	
    	for (int k = 0; k < Ants.n_ants; k++) {
    		if (Ants.ants[k].toVisit > 0) {
    			return false;
    		}
    	}
    	
    	return done;
    }

    //manage the solution construction phase (construct a set of complete and closed tours, 
    //each tour for one of the m salesmen)
    static void construct_solutions() {
		int k; /* counter variable */
		int step; /* counter of the number of construction steps */
		int salesman;
	    double nr;
		
		/* Mark all cities as unvisited */
		for (k = 0; k < Ants.n_ants; k++) {
		    Ants.ant_empty_memory(Ants.ants[k]);
		}
	
		step = 0;
		
		/* Place the ants on same initial city */
		for (k = 0; k < Ants.n_ants; k++) {
			//generate random weights (objective weights) for each ant for each objective
			//objective weights reflects the ant’s individual preferences
			nr = Math.random();
			//nr = Math.random() * 100;
			/*double nr1 = Math.random();
			double nr2 = Math.random();*/
			Ants.ants[k].weights[0] = 1 - nr;
			Ants.ants[k].weights[1] = nr;
			/*Ants.ants[k].weights[0] = nr1;
			Ants.ants[k].weights[1] = nr2;*/
			
			/* Place the ants on same initial city */
			for (int i = 0; i < MTsp.m; i++) {
				//place each ant on the depot city, which is the start city of each tour
				// -1 is a special marker for the deport city, so that it's not be confused with the rest of the cities
				// all the rest of the cities are represented by integer values > 0
				Ants.ants[k].tours[i].add(-1);   
			}
		}
	
		while (!isDone()) {
		    for (k = 0; k < Ants.n_ants; k++) {
		    	if (Ants.ants[k].toVisit > 0) {
		    		//choose for each ant in a probabilistic way by some type of roullette wheel selection 
					//which salesman to consider next, that will visit a city
		    		salesman = (int)(Math.random() * MTsp.m);
					Ants.neighbour_choose_and_move_to_next(Ants.ants[k], salesman);
					if (Ants.acs_flag)
					    Ants.local_acs_pheromone_update(Ants.ants[k], salesman);
		    	}
		    	
		    }
		}
	
		for (k = 0; k < Ants.n_ants; k++) {
			for (int i = 0; i < MTsp.m; i++) {
				step = Ants.ants[k].tours[i].size();
				Ants.ants[k].tours[i].add(step, -1);
				
				Ants.ants[k].tour_lengths[i] = Tsp.compute_tour_length_(Ants.ants[k].tours[i]);
				Ants.ants[k].total_tour_length += Ants.ants[k].tour_lengths[i];
				
			    if (Ants.acs_flag)
			    	Ants.local_acs_pheromone_update(Ants.ants[k], i);
			}
			//Ants.ants[k].total_tour_length = Tsp.compute_tour_lengths(Ants.ants[k].tours);
			Ants.ants[k].costObjectives[0] = Ants.ants[k].total_tour_length;
			Ants.ants[k].costObjectives[1] = Ants.computeToursAmplitude(Ants.ants[k]);
		}
		InOut.n_tours += (Ants.n_ants * MTsp.m); //each ant constructs a complete and closed tour
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
		
		/*
		 * Initialize the Pheromone trails, only if ACS is used, Ants.pheromones
		 * have to be initialized differently
		 */
		if (!(Ants.acs_flag)) {
			for (int nrObj = 0; nrObj < k; nrObj++) {
				Ants.trail_0[nrObj] = 1. / ((Ants.rho) * Ants.nn_tour(nrObj));
			}
		    /*
		     * in the original papers on Ant System it is not exactly defined what the
		     * initial value of the Ants.pheromones is. Here we set it to some
		     * small constant, analogously as done in MAX-MIN Ant System.
		     */
		    Ants.init_pheromone_trails(Ants.trail_0);
		}
		
		if (Ants.acs_flag) {
			for (int nrObj = 0; nrObj < k; nrObj++) {
				double length = (double) Ants.nn_tour(nrObj);
				Ants.trail_0[nrObj] = 1. / ((double) MTsp.n * length);
				/*System.out.println("Length nn" + (nrObj + 1) + ": "+ length);
				System.out.println("Initial trail " + (nrObj + 1) + ": "+ Ants.trail_0[nrObj]);	*/	
			}
			
			//bring the 2 values of initial trial to the same order of magnitude so that the
			//one wouldn't be much bigger than the other; the idea is to have close values one to the other for the initial trails
			double max = Double.MIN_VALUE, min = Double.MAX_VALUE;
			int indexMax = 0, indexMin = 0;
			for (int nrObj = 0; nrObj < k; nrObj++) {
				if (Ants.trail_0[nrObj] < min) {
					min = Ants.trail_0[nrObj];
					indexMin = nrObj;
				}
				if (Ants.trail_0[nrObj] > max) {
					max = Ants.trail_0[nrObj];
					indexMax = nrObj;
				}
			}
			double raport = max/min;
			//System.out.println("Before Raportul intre cele 2 este: " + max/min);
			int nrDigits = (int) Math.log10(raport) + 1;
			double additionalValue =  Math.pow(10, nrDigits);
			Ants.trail_0[indexMin] = Ants.trail_0[indexMin] * additionalValue;
			System.out.println("After Raportul intre cele 2 este: " + Ants.trail_0[indexMax]/Ants.trail_0[indexMin]);
			System.out.println("Min:" + Ants.trail_0[indexMin] + " max:" + Ants.trail_0[indexMax]);		
			
			Ants.indexMin = indexMin;
			Ants.addedWeight = additionalValue;
			
			Ants.init_pheromone_trails(Ants.trail_0);
		}
	
		//this is no longer possible since the total information depend on the weights chosen random
		//at each iteration by each ant, which is a dynamic information
		/* Calculate combined information Ants.pheromone times heuristic information */
		//Ants.compute_total_information();
		
    }
    
    //manage some statistical information about the trial, especially if a new best solution
    //(best-so-far) is found and adjust some parameters if a new best solution is found
    static void update_statistics(int trial, int iterationNr) {
		//int iteration_best_ant;
    	int[] antsIds = new int[TSP_ACO.k];
    	
    	/*if (iterationNr == 0 || iterationNr == InOut.max_iterations - 1) {
    		System.out.println("Trail #" + trial +  " iteration " + iterationNr + " size bestSoFarPareto=" + ParetoFront.bestSoFarPareto.size());
    	}*/
    	
		for (int indexObj = 0; indexObj < TSP_ACO.k; indexObj++) {
			antsIds = Ants.findBestSecondBest(indexObj);  	
			Ants.copy_from_to(Ants.ants[antsIds[0]], Ants.bestIterationAnts[indexObj]);
			Ants.copy_from_to(Ants.ants[antsIds[1]], Ants.secondBestIterationAnts[indexObj]);
    	}		
		
		//compute non-dominated set of solutions (iteration non-dominated front)
		ParetoFront.iterationPareto.clear();
		Ant copyAnt;
		for (int i = 0; i < Ants.n_ants; i++) {
			copyAnt = Ants.copyAnt(Ants.ants[i]);
			ParetoFront.paretoUpdateWithSolution(ParetoFront.iterationPareto, copyAnt);
		}
		
		//update BestSoFarPareto external set
		ParetoFront.paretoUpdate(ParetoFront.bestSoFarPareto, ParetoFront.iterationPareto);
		
	
		//iteration_best_ant = Ants.find_best(); /* iteration_best_ant is a global variable */
	
		/*if (Ants.ants[iteration_best_ant].total_tour_length < Ants.best_so_far_ant.total_tour_length) {
	
		    InOut.time_used = Timer.elapsed_time();  best solution found after time_used 
		    Ants.copy_from_to(Ants.ants[iteration_best_ant], Ants.best_so_far_ant);
	
		    InOut.found_best = InOut.iteration;	 
		}*/
		
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
		int k;
	
		for (k = 0; k < Ants.n_ants; k++)
		    Ants.global_update_pheromone(Ants.ants[k]);
    }

    //manage global Ants.pheromone deposit for Ant Colony System
    static void acs_global_update() {
    	Map<Integer, Integer>[] edgesMapArray;
    	
    	//global pheromone update is performed by the best iteration and second best iteration ants
    	//in each pheromone matrix corresponding to the objectives
    	for (int indexObj = 0; indexObj < TSP_ACO.k; indexObj++) {
    		edgesMapArray = Ants.global_acs_pheromone_update1(Ants.bestIterationAnts[indexObj], indexObj);	
		    Ants.global_acs_pheromone_update2(Ants.secondBestIterationAnts[indexObj], indexObj, edgesMapArray);	
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
	        
			InOut.init_program(args);	
			
		    MTsp.instance.nn_list = Tsp.compute_nn_lists();
			
		    Ants.pheromone = new double[k][][];
			//Ants.total = new double[k][][];
			
			for (int i = 0; i < k; i++) {
				Ants.pheromone[i] = new double[MTsp.n + 1][MTsp.n + 1];
				//Ants.total[i] = new double[MTsp.n + 1][MTsp.n + 1];
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
