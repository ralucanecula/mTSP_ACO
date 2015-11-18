package multiobjectiveACO;

import java.util.ArrayList;
import java.util.Map;

import multiobjectiveACO.Ants.Ant;

/* this Eclipse project represents the ACO algorithm adapted from the mTSP_ACO global 
 * project in Eclipse, but in which there are considered more than one objective, such that 
 * a multiobjective ACO algorithm called MACS (Multiple ant colony system) is applied 
 * to solve the mTSP problem and find a good and diversified set of non-dominated solutions 
 * belonging to the Pareto front; the 2 considered objectives that must be optimized simultaneously 
 * are: 1. minimizing the total distance of tours traveled by the salesmen and 2. balancing the
 * tours traveled by the salesmen, i.e. minimizing the length of the longest tour or minimizing the
 * amplitude, which is the difference between the length of the longest tour and the shortest tour 
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
    
    public static int foundBestParetoSet = 0;

	
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
		
		/* Mark all cities as unvisited */
		for (k = 0; k < Ants.n_ants; k++) {
		    Ants.ant_empty_memory(Ants.ants[k]);
		}
	
		step = 0;
		
		/* Place the ants on same initial city and compute weight for each ant */
		for (k = 0; k < Ants.n_ants; k++) {
			Ants.ants[k].weight = (double)k / (double)Ants.n_ants;
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
				Ants.trail_0 = 1. / ((Ants.rho) * Ants.nn_tour());	
		    /*
		     * in the original papers on Ant System it is not exactly defined what the
		     * initial value of the Ants.pheromones is. Here we set it to some
		     * small constant, analogously as done in MAX-MIN Ant System.
		     */
		    Ants.init_pheromone_trails(Ants.trail_0);
		}
		
		if (Ants.acs_flag) {
			Ants.trail_0 = 1. / ((double) MTsp.n * Ants.nn_tour());
			
			Ants.init_pheromone_trails(Ants.trail_0);
		}
	
		//this is no longer possible since the total information depend on the weights chosen random
		//at each iteration by each ant, which is a dynamic information
		/* Calculate combined information Ants.pheromone times heuristic information */
		//Ants.compute_total_information();
		
    }
    
    //manage some statistical information about the trial, especially if a new best solution
    //(best-so-far) is found and adjust some parameters if a new best solution is found
    static boolean update_statistics(int trial, int iterationNr) {
    	boolean doPheromoneGlobalUpdate = false;
    	double averageValues[] = new double[TSP_ACO.k];  
    		
    	/*if (iterationNr == 0 || iterationNr == InOut.max_iterations - 1) {
    		System.out.println("Trail #" + trial +  " iteration " + iterationNr + " size bestSoFarPareto=" + ParetoFront.bestSoFarPareto.size());
    	}*/	
		
		//compute non-dominated set of solutions (iteration non-dominated front)
		ParetoFront.iterationPareto.clear();
		Ant copyAnt;
		for (int i = 0; i < Ants.n_ants; i++) {
			copyAnt = Ants.copyAnt(Ants.ants[i]);
			ParetoFront.paretoUpdateWithSolution(ParetoFront.iterationPareto, copyAnt);
		}
		
		//update BestSoFarPareto external set
		ParetoFront.paretoUpdate(ParetoFront.bestSoFarPareto, ParetoFront.iterationPareto);
		
		//compute average value on each objective function
		for (int indexObj = 0; indexObj < TSP_ACO.k; indexObj++) {
			for (int i = 0; i < ParetoFront.bestSoFarPareto.size(); i++) {
				averageValues[indexObj] += ParetoFront.bestSoFarPareto.get(i).costObjectives[indexObj];
			}
			averageValues[indexObj] = averageValues[indexObj] / (double)ParetoFront.bestSoFarPareto.size();			
		}
		
	    //compute new value of initial trail information
		Ants.trail_0_prim = 1.0 / ((double) MTsp.n * averageValues[0] * averageValues[1]);
		//a better Pareto set P was found -> improve exploration
		if (Ants.trail_0_prim > Ants.trail_0) {
			doPheromoneGlobalUpdate = false;
			//System.out.println("A better Pareto set was found");
			foundBestParetoSet++;
		}
		//exploitation is favored by globally updating the pheromone trail
		else {
			doPheromoneGlobalUpdate = true;
		}	
		
		return doPheromoneGlobalUpdate;
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
		for (int i = 0; i < ParetoFront.bestSoFarPareto.size(); i++) {
			Ants.global_update_pheromone(ParetoFront.bestSoFarPareto.get(i));
		}
		    
    }

    //manage global Ants.pheromone deposit for Ant Colony System
    static void acs_global_update() {	
    	//global pheromone update is performed with each solution of the current Pareto optimal set P 
    	for (int i = 0; i < ParetoFront.bestSoFarPareto.size(); i++) {
    		 Ants.global_acs_pheromone_update(ParetoFront.bestSoFarPareto.get(i));	
    	}
    	
    }

    //manage global Ants.pheromone trail update for the ACO algorithms
	static void pheromone_trail_update(boolean performGlobalUpdate)
	{
		//a better Pareto set P was found 
		if (!performGlobalUpdate) {
			Ants.trail_0 = Ants.trail_0_prim;
			//the pheromone trail matrix is reinitialized with the new value
			Ants.init_pheromone_trails(Ants.trail_0);
		}
		else {
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
	        
			InOut.init_program(args);	
			
		    MTsp.instance.nn_list = Tsp.compute_nn_lists();
			
		    Ants.pheromone = new double[MTsp.n + 1][MTsp.n + 1];
			
			init_try(); 
	
			boolean saveIterCosts = false;
			ArrayList<Double> iterTotalCost = null;
			if (saveIterCosts) {
				//for saving detailed cost values at each 5 iteration
				iterTotalCost = new ArrayList<Double>();
			}
			
			InOut.iteration = 0;
			boolean performGlobalUpdate;
		    while (!termination_condition()) {
				construct_solutions();
				performGlobalUpdate = update_statistics(trial, InOut.iteration);
				pheromone_trail_update(performGlobalUpdate);
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
