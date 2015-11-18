package multiobjectiveACO;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

public class Ants {

    static class Ant {
    	//for each of the m salesmen an ant will construct a tour, so that a candidate solution constructed by
    	//an ant will be represented by a list of tours, one for each salesman
    	ArrayList<Integer>[] tours;
		boolean[] visited;
		double[] tour_lengths;
		double total_tour_length;
		//cities left to be visited by an ant (initially toVisit = n, which is the number of cities from the mTSP instance)
		int toVisit;
		//stores the cost of each solution according to the considered objectives (2 in this case)
		double costObjectives[];
		//stores the vector of weights to be used in the transition formula when computing the
		//weighted sum of the pheromone matrices for each objective
		double weights[];
    }

    public static final int MAX_ANTS = 1024;
    public static final int MAX_NEIGHBOURS = 512;

    static Ant ants[];
    //static Ant best_so_far_ant;
    
    //used in the global pheromone update of the pheromone matrices for each objective
    //therefore in case of 2 objectives there will be in total 2 bestIterationAnts and 2 secondBestIterationAnts
    static Ant bestIterationAnts[];
    static Ant secondBestIterationAnts[];

    static double pheromone[][][];
    //static double total[][][];  //keeps heuristic information times pheromone for each arc

    static double prob_of_selection[];

    static int n_ants; /* number of ants */
 
    static int nn_ants; /* length of nearest neighbor lists for the ants' solution construction */

    static double rho; /* parameter for evaporation */
    static double alpha; /* importance of trail */
    static double beta; /* importance of heuristic evaluate */
    static double q_0; /* probability of best choice in tour construction */

    static boolean as_flag; /* ant system */
    static boolean acs_flag; /* ant colony system (ACS) */

    static int u_gb; /* every u_gb iterations update with best-so-far ant */

    //according to the considered objectives there will be 2 different values for the initial trail in case of 2 objectives
    static double trail_0[] = new double[TSP_ACO.k]; /* initial pheromone level in ACS */
    
    //index of the initial trail which has the minimum value
    static int indexMin = 0;
    
    //addtional weight to be added to the value of the pheromones from the matrix where the value of
    //1/f(cost according to objective) is smaller than the other one; the intent is to reduce the
    //difference between their values
    static double addedWeight = 0;

/*    static double HEURISTIC(int m, int n) {
    	return (1.0 / (double) MTsp.instance.distance[m][n]);
    }
*/
    static double HEURISTIC(int m, int n) {
    	return (double) MTsp.instance.distance[m][n];
    }
    
    //allocate the memory for the ant colony, the best-so-far ant
    static void allocate_ants() {
		int i, j;
	
		ants = new Ant[n_ants];
	
		for (i = 0; i < n_ants; i++) {
		    ants[i] = new Ant();
		    ants[i].tours = (ArrayList<Integer>[])new ArrayList[MTsp.m];
		    ants[i].tour_lengths = new double[MTsp.m];
		    for (j = 0; j < MTsp.m; j++) {
		    	ants[i].tours[j] = new ArrayList<Integer>();
		    	ants[i].tour_lengths[j] = 0;
		    }
		    ants[i].visited = new boolean[MTsp.n];
		    ants[i].toVisit = MTsp.n;
		    
		    ants[i].costObjectives = new double[TSP_ACO.k];
		    ants[i].weights = new double[TSP_ACO.k];
		    for (int indexObj = 0; indexObj < TSP_ACO.k; indexObj++) {
		    	ants[i].costObjectives[indexObj] = 0;
		    	ants[i].weights[indexObj] = 0;
	    	}

		}
		
		bestIterationAnts = new Ant[TSP_ACO.k];
		secondBestIterationAnts = new Ant[TSP_ACO.k];
		
		for (int indexObj = 0; indexObj < TSP_ACO.k; indexObj++) {
			bestIterationAnts[indexObj] = new Ant();
			bestIterationAnts[indexObj].tours = (ArrayList<Integer>[])new ArrayList[MTsp.m];
			bestIterationAnts[indexObj].tour_lengths = new double[MTsp.m];
		    for (j = 0; j < MTsp.m; j++) {
		    	bestIterationAnts[indexObj].tours[j] = new ArrayList<Integer>();
		    	bestIterationAnts[indexObj].tour_lengths[j] = 0;
		    }
		    bestIterationAnts[indexObj].visited = new boolean[MTsp.n];
		    bestIterationAnts[indexObj].toVisit = MTsp.n;
		    bestIterationAnts[indexObj].costObjectives = new double[TSP_ACO.k];
		    bestIterationAnts[indexObj].weights = new double[TSP_ACO.k];
		    bestIterationAnts[indexObj].costObjectives[indexObj] = 0;
		    bestIterationAnts[indexObj].weights[indexObj] = 0;
	    	    
		    secondBestIterationAnts[indexObj] = new Ant();
		    secondBestIterationAnts[indexObj].tours = (ArrayList<Integer>[])new ArrayList[MTsp.m];
		    secondBestIterationAnts[indexObj].tour_lengths = new double[MTsp.m];
		    for (j = 0; j < MTsp.m; j++) {
		    	secondBestIterationAnts[indexObj].tours[j] = new ArrayList<Integer>();
		    	secondBestIterationAnts[indexObj].tour_lengths[j] = 0;
		    }
		    secondBestIterationAnts[indexObj].visited = new boolean[MTsp.n];
		    secondBestIterationAnts[indexObj].toVisit = MTsp.n;
		    secondBestIterationAnts[indexObj].costObjectives = new double[TSP_ACO.k];
		    secondBestIterationAnts[indexObj].weights = new double[TSP_ACO.k];
		    secondBestIterationAnts[indexObj].costObjectives[indexObj] = 0;
		    secondBestIterationAnts[indexObj].weights[indexObj] = 0;
			
		}
		/*best_so_far_ant = new Ant();
		best_so_far_ant.tours = (ArrayList<Integer>[])new ArrayList[MTsp.m];
		best_so_far_ant.tour_lengths = new double[MTsp.m];
	    for (j = 0; j < MTsp.m; j++) {
	    	best_so_far_ant.tours[j] = new ArrayList<Integer>();
	    	best_so_far_ant.tour_lengths[j] = 0;
	    }
		best_so_far_ant.visited = new boolean[MTsp.n];
		best_so_far_ant.toVisit = MTsp.n;*/
	
		prob_of_selection = new double[nn_ants + 1];
		for (i = 0; i < nn_ants + 1; i++) {
		    prob_of_selection[i] = Double.POSITIVE_INFINITY;
		}
    }

    // find the best ant of the current iteration (the one with the minimum total tour length)
    static int find_best()
    {
    	double min;
		int k, k_min;
	
		min = ants[0].total_tour_length;
		k_min = 0;
		for (k = 1; k < n_ants; k++) {
		    if (ants[k].total_tour_length < min) {
				min = ants[k].total_tour_length;
				k_min = k;
		    }
		}
		return k_min;
    }
    
    //find the best and the second best ants in the current iteration according to a specified objective
    static int[] findBestSecondBest(int indexObj) {
    	double min1, min2;
		int k, k_min1, k_min2;
		int result[] = new int[2];
	
		min1 = min2 = Integer.MAX_VALUE;
		k_min1 = 0; k_min2 = 0;
		for (k = 0; k < n_ants; k++) {
			/* If current element is smaller than first then update both first and second */
			if (ants[k].costObjectives[indexObj] < min1) {
				min2 = min1;
				min1 = ants[k].costObjectives[indexObj];
				k_min2 = k_min1;
				k_min1 = k;
			}
			
			/* If is in between first and second then update second  */
			else if (ants[k].costObjectives[indexObj] < min2 && ants[k].costObjectives[indexObj] != min1) {
				min2 = ants[k].costObjectives[indexObj];
				k_min2 = k;
			}
		}		
		result[0] =  k_min1; result[1] =  k_min2;
		
		return result;		
    }

    //initialize pheromone trails for each considered objective
    //matricea cu urmele de feromoni trebuie sa se faca relativ la toate cele n orase
    static void init_pheromone_trails(double[] initial_trail)
    {
		int i, j;
	
		//for each objective 
		for (int indexObj = 0; indexObj < TSP_ACO.k; indexObj++) {
			/* Initialize pheromone trails */
			for (i = 0; i < MTsp.n + 1; i++) {
			    for (j = 0; j <= i; j++) {
					pheromone[indexObj][i][j] = initial_trail[indexObj];
					pheromone[indexObj][j][i] = initial_trail[indexObj];
					/*total[indexObj][i][j] = initial_trail[indexObj];
					total[indexObj][j][i] = initial_trail[indexObj];*/
			    }
			}
		}
		
    }

    //implements the pheromone trail evaporation
    static void evaporation()
    {
		int i, j;
		
		for (int indexObj = 0; indexObj < TSP_ACO.k; indexObj++) {
		    for (i = 0; i < MTsp.n + 1; i++) {
			    for (j = 0; j <= i; j++) {
					pheromone[indexObj][i][j] = (1 - rho) * pheromone[indexObj][i][j];
					pheromone[indexObj][j][i] = pheromone[indexObj][i][j];
			    }
		    }		
    	}
	
    }

    //reinforces edges used in ant k's solution
    static void global_update_pheromone(Ant a)
    {
		int i, j, h, k, size;
		double d_tau;

		for (int indexObj = 0; indexObj < TSP_ACO.k; indexObj++) {
		    d_tau = 1.0 / (double) a.costObjectives[indexObj];
			for (i = 0; i < MTsp.m; i++) {
				size = a.tours[i].size();
				for (k = 0; k < size - 1; k++) {
				    j = a.tours[i].get(k); 
				    h = a.tours[i].get(k + 1);
				    
				    j++;
		            h++;
		            
				    pheromone[indexObj][j][h] += d_tau;
				    pheromone[indexObj][h][j] = pheromone[indexObj][j][h];
				}
			} 		
    	}
		
    }

    //calculates heuristic info times pheromone for each arc
    /*static void compute_total_information()
    {
		int i, j;
	
		for (i = 0; i < MTsp.n + 1; i++) {
		    for (j = 0; j < i; j++) {
				total[i][j] = Math.pow(pheromone[i][j], alpha) * Math.pow(HEURISTIC(i, j), beta);
				total[j][i] = total[i][j];
		    }
		}
    }*/

    //empty the ants's memory regarding visited cities
    static void ant_empty_memory(Ant a)
    {
        int i, j;

		a.total_tour_length = 0;
	    for (int indexObj = 0; indexObj < TSP_ACO.k; indexObj++) {
	    	a.costObjectives[indexObj] = 0;
	    	a.weights[indexObj] = 0.0;
	    }
		
		//clear all the elements (cities) from the m tours of an ant
		for (i = 0; i < MTsp.m; i++) {
		    a.tours[i].clear();
		    a.tour_lengths[i] = 0;
		}
		
		for (j = 0; j < MTsp.n; j++) {
		    a.visited[j] = false;
		}
		a.toVisit = MTsp.n;
		
    }

    //choose for an ant as the next city the one with maximal value of heuristic information times pheromone
    static void choose_best_next(Ant a, int indexSalesman)
    {
		int city, current_city, next_city;
		double value_best;
		double help, weightedSum = 0, heuristicValue = 0, amplitude;
	
		next_city = MTsp.n;
		int lastPos = a.tours[indexSalesman].size() - 1;
		current_city = a.tours[indexSalesman].get(lastPos);
    	current_city++;
    	
		value_best = -1.; /* values in total matrix are always >= 0.0 */
		for (city = 0; city < MTsp.n; city++) {
			weightedSum = 0;
			heuristicValue = 0;
			amplitude = 0;
		    if (a.visited[city])
			; /* city already visited, do nothing */
		    else {
		    	for (int indexObj = 0; indexObj < TSP_ACO.k; indexObj++) {
		    		weightedSum += a.weights[indexObj] * pheromone[indexObj][current_city][city + 1];		    		
		    	}
		    	//add temporary the node under verification to its corresponding subtour to take it 
		    	//into account when computing the amplitude of the partial constructed subtours
		    	a.tours[indexSalesman].add(city);
	
		    	//compute lengths for all partial subtours constructed so far so as to compute 
		    	//their amplitude
		    	for (int i1 = 0; i1 < MTsp.m; i1++) {
		    		a.tour_lengths[i1] = Tsp.compute_tour_length_(a.tours[i1]);
		    	}
		    	amplitude = computeToursAmplitude(a);
		    	
		    	//since the value for amplitude is much greater than the value of HEURISTIC(current_city, help_city)
		    	//we must bring the value of amplitude on the same scale with value of HEURISTIC(current_city, help_city)
		    	if (amplitude > 0) {
		    		double raport = amplitude / HEURISTIC(current_city, city + 1);
			    	int nrDigits = (int) Math.log10(raport) + 1;
					double additionalValue =  Math.pow(10, nrDigits);
					amplitude = amplitude / additionalValue;
		    	}
		    	
		    	//remove the node that we've previously added temporary
		    	int lastPos1 = a.tours[indexSalesman].size() - 1;
		    	a.tours[indexSalesman].remove(lastPos1);
		    	
		    	heuristicValue = ((double)TSP_ACO.k) / (amplitude + HEURISTIC(current_city, city + 1));
		    	
		    	help = Math.pow(weightedSum, alpha) * Math.pow(heuristicValue, beta);
				if (help > value_best) {
				    next_city = city;
				    value_best = help;
				}
		    }
		}
		a.tours[indexSalesman].add(next_city);
		a.visited[next_city] = true;
		a.toVisit--;
    }

    //chooses for an ant as the next city the one with maximal value of heuristic information times pheromone
    static void neighbour_choose_best_next(Ant a, int indexSalesman) {
		int i, current_city, next_city, help_city;
		double value_best, help;
	    double weightedSum = 0, heuristicValue = 0, amplitude;
		
		next_city = MTsp.n;   //next_city = Integer.MAX_VALUE;
		int lastPos = a.tours[indexSalesman].size() - 1;
		current_city = a.tours[indexSalesman].get(lastPos);
        current_city++;
 
		value_best = -1.;  //values in total matrix are always >= 0.0 
		for (i = 0; i < nn_ants; i++) {
			weightedSum = 0;
			heuristicValue = 0;
			amplitude = 0;
		    help_city = MTsp.instance.nn_list[current_city][i];
		    if (a.visited[help_city - 1])
		    	; // city already visited, do nothing 
		    else {
		    	for (int indexObj = 0; indexObj < TSP_ACO.k; indexObj++) {
		    		weightedSum += a.weights[indexObj] * pheromone[indexObj][current_city][help_city];		    		
		    	}
		    	//add temporary the node under verification to its corresponding subtour to take it 
		    	//into account when computing the amplitude of the partial constructed subtours
		    	a.tours[indexSalesman].add(help_city - 1);
	
		    	//compute lengths for all partial subtours constructed so far so as to compute 
		    	//their amplitude
		    	for (int i1 = 0; i1 < MTsp.m; i1++) {
		    		a.tour_lengths[i1] = Tsp.compute_tour_length_(a.tours[i1]);
		    	}
		    	amplitude = computeToursAmplitude(a);
		    	
		    	/*if (amplitude >= HEURISTIC(current_city, help_city)) {
		    		System.out.println("Amplitude is greater");
		    	}
		    	else {
		    		System.out.println("Heuristic is greater");
		    	}*/
		    	    	
		    	//since the value for amplitude is much greater than the value of HEURISTIC(current_city, help_city)
		    	//we must bring the value of amplitude on the same scale with value of HEURISTIC(current_city, help_city)
		    	if (amplitude > 0) {
		    		double raport = amplitude / HEURISTIC(current_city, help_city);
			    	int nrDigits = (int) Math.log10(raport) + 1;
					double additionalValue =  Math.pow(10, nrDigits);
					amplitude = amplitude / additionalValue;
					//raport = amplitude / HEURISTIC(current_city, help_city);
			    	//System.out.println("After Raportul intre cele 2 euristici este: " + raport);
		    	}
		    		    	
		    	//remove the node that we've previously added temporary
		    	int lastPos1 = a.tours[indexSalesman].size() - 1;
		    	a.tours[indexSalesman].remove(lastPos1);		    	
		    	
		    	heuristicValue = ((double)TSP_ACO.k) / (amplitude + HEURISTIC(current_city, help_city));
		    	
		    	help = Math.pow(weightedSum, alpha) * Math.pow(heuristicValue, beta);
				//help = total[current_city][help_city];
				if (help > value_best) {
				    value_best = help;
				    next_city = help_city - 1;
				}
		    }
		}
		
		if (next_city == MTsp.n)
		    // all cities in nearest neighbor list were already visited 
		    choose_best_next(a, indexSalesman);
		else {
			a.tours[indexSalesman].add(next_city);
		    a.visited[next_city] = true;
		    a.toVisit--;
		}
    }
    
    static void choose_closest_next(Ant a, int indexSalesman)
    {
		int current_city, next_city, city;
	    double  min_distance;
		
		next_city = MTsp.n;
		int lastPos = a.tours[indexSalesman].size() - 1;
		current_city = a.tours[indexSalesman].get(lastPos);
		current_city++;
		
		min_distance = Integer.MAX_VALUE;  //Search shortest edge 
		for (city = 0; city < MTsp.n; city++) {
		    if (a.visited[city])
			;  //city already visited 
		    else {
				if (MTsp.instance.distance[current_city][city + 1] < min_distance) {
				    next_city = city;
				    min_distance = MTsp.instance.distance[current_city][city + 1];
				}
		    }
		}
		a.tours[indexSalesman].add(next_city);
		a.visited[next_city] = true; 
		a.toVisit--;
    }

    //Choose for an ant probabilistically a next city among all unvisited cities in the current city's candidate list
    static void neighbour_choose_and_move_to_next(Ant a, int indexSalesman) {
		int i, help, city;
		int current_city;
		double rnd, partial_sum = 0., sum_prob = 0.0;
		double prob_ptr[];
	
		if ((q_0 > 0.0) && (Utilities.random01() < q_0)) {
		    /*
		     * with a probability q_0 make the best possible choice
		     * according to pheromone trails and heuristic information, this corresponds to exploitation
		     */
		    /*
		     * we first check whether q_0 > 0.0, to avoid the very common case
		     * of q_0 = 0.0 to have to compute a random number, which is
		     * expensive computationally
		     */
		    neighbour_choose_best_next(a, indexSalesman);
		    return;
		}
	
		prob_ptr = prob_of_selection;
	
	    /* current_city city of ant k */
		int lastPos = a.tours[indexSalesman].size() - 1;
		current_city = a.tours[indexSalesman].get(lastPos);
    	current_city++;
    	double weightedSum = 0, heuristicValue = 0, amplitude;
    
		for (i = 0; i < nn_ants; i++) {
			weightedSum = 0;
			heuristicValue = 0;
			amplitude = 0;
			city = MTsp.instance.nn_list[current_city][i];
		    if (a.visited[city - 1])
		    	prob_ptr[i] = 0.0; /* city already visited */
		    else {
		    	for (int indexObj = 0; indexObj < TSP_ACO.k; indexObj++) {
		    		weightedSum += a.weights[indexObj] * pheromone[indexObj][current_city][city];		    		
		    	}
		    	//add temporary the node under verification to its corresponding subtour to take it 
		    	//into account when computing the amplitude of the partial constructed subtours
		    	a.tours[indexSalesman].add(city - 1);
	
		    	//compute lengths for all partial subtours constructed so far so as to compute 
		    	//their amplitude
		    	for (int i1 = 0; i1 < MTsp.m; i1++) {
		    		a.tour_lengths[i1] = Tsp.compute_tour_length_(a.tours[i1]);
		    	}
		    	amplitude = computeToursAmplitude(a);
		    	
		    	//since the value for amplitude is much greater than the value of HEURISTIC(current_city, help_city)
		    	//we must bring the value of amplitude on the same scale with value of HEURISTIC(current_city, help_city)
		    	if (amplitude > 0) {
		    		double raport = amplitude / HEURISTIC(current_city, city);
			    	int nrDigits = (int) Math.log10(raport) + 1;
					double additionalValue =  Math.pow(10, nrDigits);
					amplitude = amplitude / additionalValue;
		    	}
		    	
		    	//remove the node that we've previously added temporary
		    	int lastPos1 = a.tours[indexSalesman].size() - 1;
		    	a.tours[indexSalesman].remove(lastPos1);
		    	
		    	heuristicValue = ((double)TSP_ACO.k) / (amplitude + HEURISTIC(current_city, city));
		    	
				prob_ptr[i] = Math.pow(weightedSum, alpha) * Math.pow(heuristicValue, beta);
				sum_prob += prob_ptr[i];
		    }
		}
	
		if (sum_prob <= 0.0) {
		    /* All cities from the candidate set are tabu (are already visited) */
		    choose_best_next(a, indexSalesman);
		} else {
		    /*
		     * at least one neighbor is eligible, choose one according to the
		     * selection probabilities
		     */
		    rnd = Utilities.random01();
		    rnd *= sum_prob;
		    i = 0;
		    partial_sum = prob_ptr[i];
		    /* This loop always stops because prob_ptr[nn_ants] == HUGE_VAL */
		    while (partial_sum <= rnd) {
				i++;
				partial_sum += prob_ptr[i];
		    }
		    /*
		     * This may very rarely happen because of rounding if rnd is close to 1.
		     */
		    if (i == nn_ants) {
				neighbour_choose_best_next(a, indexSalesman);
				return;
		    }
		    help = MTsp.instance.nn_list[current_city][i];
		    a.tours[indexSalesman].add(help - 1);
		    a.visited[help - 1] = true;
		    a.toVisit--;
		}
    }

    //reinforces the edges used in ant's solution as in ACS
    //global pheromone update performed by the best iteration ant
    static Map<Integer, Integer>[] global_acs_pheromone_update1(Ant a, int indexObj) {
		int i, j, h, k, size;
		double d_tau;
		Map<Integer, Integer>[] edgesMapArray = (HashMap<Integer, Integer>[])new HashMap[MTsp.m];
		
		for (i = 0; i < MTsp.m; i++) {
    		edgesMapArray[i] = new HashMap<Integer, Integer>();
    	}
		
		d_tau = 1.0 / (double) a.costObjectives[indexObj];
		if (indexObj == Ants.indexMin) {
			d_tau = d_tau * Ants.addedWeight;
		}
	
		for (i = 0; i < MTsp.m; i++) {
			size = a.tours[i].size();
			for (k = 0; k < size - 1; k++)  {
			    j = a.tours[i].get(k);
			    h = a.tours[i].get(k + 1);	
			    
			    edgesMapArray[i].put(j, h);
			    
			    j++;
			    h++;
		
			    pheromone[indexObj][j][h] = (1. - rho) * pheromone[indexObj][j][h] + rho * d_tau;
			    pheromone[indexObj][h][j] = pheromone[indexObj][j][h];
		
/*			    total[h][j] = Math.pow(pheromone[h][j], alpha) * Math.pow(HEURISTIC(h, j), beta);
			    total[j][h] = total[h][j];*/
			}
		}
		
		return edgesMapArray;
    }
    
    //check if edge (i, j) is used in the tours of the best iteration ant
    static boolean checkBestIterationTour(int startNode, int endNode, Map<Integer, Integer>[] edgesMapArray) {
    	boolean result = false;
    	int node;
    	
    	for (int i = 0; i < edgesMapArray.length; i++) {
    		if (edgesMapArray[i].containsKey(startNode)) {
    			node = edgesMapArray[i].get(startNode);
    			if (node == endNode) {
    				return true;
    			}
    		}
    	}
    	
    	return result;  	
    }
    
    //reinforces the edges used in ant's solution as in ACS
    //global pheromone update performed by the second best iteration ant
    static void global_acs_pheromone_update2(Ant a, int indexObj, Map<Integer, Integer>[] edgesMapArray) {
		int i, j, h, k, size;
		double d_tau;
		boolean belongs = false;
	
		d_tau = 1.0 / (double) a.costObjectives[indexObj];
		if (indexObj == Ants.indexMin) {
			d_tau = d_tau * Ants.addedWeight;
		}
	
		for (i = 0; i < MTsp.m; i++) {
			size = a.tours[i].size();
			for (k = 0; k < size - 1; k++)  {
			    j = a.tours[i].get(k);
			    h = a.tours[i].get(k + 1);	
			    
			    //check if edge (j, h) contained in a tour of the second best iteration ant it is 
			    //also used in the tours of the best iteration ant
			    belongs = checkBestIterationTour(j, h, edgesMapArray);
			 
			    j++;
			    h++;
	 
			    if (belongs) {
			    	pheromone[indexObj][j][h] += rho * d_tau;
				    pheromone[indexObj][h][j] += rho * d_tau;
			    }
			    else {
			    	pheromone[indexObj][j][h] = (1. - rho) * pheromone[indexObj][j][h] + rho * d_tau;
			    	pheromone[indexObj][h][j] = pheromone[indexObj][j][h];
			    }
		
/*			    total[h][j] = Math.pow(pheromone[h][j], alpha) * Math.pow(HEURISTIC(h, j), beta);
			    total[j][h] = total[h][j];*/
			}
		}
		
    }

    //removes some pheromone on edge just passed by the ant
    static void local_acs_pheromone_update(Ant a, int indexSalesman) {
		int h, j;
		
		int lastPos = a.tours[indexSalesman].size() - 1;
		j = a.tours[indexSalesman].get(lastPos);
		h = a.tours[indexSalesman].get(lastPos - 1);	
		
		j++;
		h++;

		//perform local pheromone update on each pheromone matrix corresponding to an objective
		for (int indexObj = 0; indexObj < TSP_ACO.k; indexObj++) {
			/* still additional parameter has to be introduced */
			pheromone[indexObj][h][j] = (1. - 0.1) * pheromone[indexObj][h][j] + 0.1 * trail_0[indexObj];
			pheromone[indexObj][j][h] = pheromone[indexObj][h][j];
			/*total[h][j] = Math.pow(pheromone[h][j], alpha) * Math.pow(HEURISTIC(h, j), beta);
			total[j][h] = total[h][j];	 */   		
    	}
		
    }
    
    //create a copy of the ant a and return the created copy at the output
    static Ant copyAnt(Ant a) {
    	//first create an empty ant
    	Ant copy = new Ant();
    	copy.tours = (ArrayList<Integer>[])new ArrayList[MTsp.m];
    	copy.tour_lengths = new double[MTsp.m];
	    for (int j = 0; j < MTsp.m; j++) {
	    	copy.tours[j] = new ArrayList<Integer>();
	    	copy.tour_lengths[j] = 0;
	    }
	    copy.visited = new boolean[MTsp.n];
	    copy.toVisit = MTsp.n;
	    
	    copy.costObjectives = new double[TSP_ACO.k];
	    copy.weights = new double[TSP_ACO.k];
	    /*for (int indexObj = 0; indexObj < TSP_ACO.k; indexObj++) {
	    	copy.costObjectives[indexObj] = 0;
	    	copy.weights[indexObj] = 0;
    	}*/
    	
	    //then copy the information from the ant a
	    copy.total_tour_length = a.total_tour_length;
		for (int indexObj = 0; indexObj < TSP_ACO.k; indexObj++) {
			copy.costObjectives[indexObj] = a.costObjectives[indexObj];
			copy.weights[indexObj] = a.weights[indexObj];
    	}
		for (int i = 0; i < MTsp.m; i++) {
			copy.tour_lengths[i] = a.tour_lengths[i];
			int size = a.tours[i].size();
			copy.tours[i] = new ArrayList<Integer>(size);
			for (int j = 0; j < size; j++) {
				int elem = a.tours[i].get(j);
				copy.tours[i].add(elem);
			}		
		}
    	
    	return copy;
    }

    //copy solution from ant a1 into ant a2
    static void copy_from_to(Ant a1, Ant a2) {
		int i, j;
	
		Ants.ant_empty_memory(a2);
		
		a2.total_tour_length = a1.total_tour_length;
		for (int indexObj = 0; indexObj < TSP_ACO.k; indexObj++) {
		   a2.costObjectives[indexObj] = a1.costObjectives[indexObj];
		   a2.weights[indexObj] = a1.weights[indexObj];
    	}
		for (i = 0; i < MTsp.m; i++) {
			a2.tour_lengths[i] = a1.tour_lengths[i];
			int size = a1.tours[i].size();
			a2.tours[i] = new ArrayList<Integer>(size);
			for (j = 0; j < size; j++) {
				int elem = a1.tours[i].get(j);
				a2.tours[i].add(elem);
			}
			
		}
    }
    
    static double computeToursAmplitude(Ant a) {
    	double min, max;
		int i;

		min = a.tour_lengths[0];
		max = a.tour_lengths[0];
		for (i = 1; i < a.tours.length; i++) {
		    if (a.tour_lengths[i] < min) {
				min = a.tour_lengths[i];
		    }
		    if (a.tour_lengths[i] > max) {
				max = a.tour_lengths[i];
		    }	    
		}
		
		return (max - min);
    }

    //generate some nearest neighbor tour and compute tour length according to a specific objective
    static double nn_tour(int indexObjective) {
    	int step, salesman;
    	double sum = 0;
    	double amplitude = 0;
    	double result = 0;
    	
    	ant_empty_memory(ants[0]);
    	step = 0;
    	
    	for (int i = 0; i < MTsp.m; i++) {
    		//place the ant on the depot city, which is the start city of each tour
			// -1 is a special marker for the deport city, so that it's not be confused with the rest of the cities
			// all the rest of the cities are represented by integer values > 0
    		ants[0].tours[i].add(-1);  	
    	}
		
		while (Ants.ants[0].toVisit > 0) {   //there are still left cities to be visited	
    		//choose in a random manner with equal probability which salesman to consider next, that will visit a city
	    	salesman = (int)(Math.random() * MTsp.m);
	    	choose_closest_next(ants[0], salesman);	 
		}
		
		for (int i = 0; i < MTsp.m; i++) {
			step = ants[0].tours[i].size();
			ants[0].tours[i].add(step, -1);
			
			ants[0].tour_lengths[i] = Tsp.compute_tour_length_(ants[0].tours[i]);
			ants[0].total_tour_length += ants[0].tour_lengths[i];
		}
		sum = ants[0].total_tour_length;
		
		if (indexObjective == 0) {
			result = sum;
		}
		else if (indexObjective == 1) {
			amplitude = computeToursAmplitude(ants[0]);
			result = amplitude;
		}
		ant_empty_memory(ants[0]);
		
		return result;
    }


    //generate some nearest neighbor tour and compute tour length
    /*static int nn_tour_(Partition[] clusters) 
    {
		int phase, sum = 0;
		
		//for each of the m crisp partitions obtained with fuzzy cMeans clustering algorithm, that define a 
		//TSP problem, a nearest neighbor tour should be constructed and at the end, the total sum of these m
		//tour lengths should be computed and returned as the length of the nn_tour
	    for (int i = 0; i < clusters.length; i++) {	    	
	    	ant_empty_memory(ants[0]);
	    	phase = 0;  //counter of the construction steps 
			//place_ant(i, ants[0], phase);
	    	ants[0].tours[i].add(phase, -1);
		
			while (phase < MTsp.clusterCities[i].size()) {
			    phase++;
			    choose_closest_next_(i, ants[0], phase);
			}
			phase = MTsp.clusterCities[i].size() + 1;
			ants[0].tours[i].add(phase, -1);
			//we assume the tour contains at the first and last position the depot city, so that the total length of the tour 
			//should take into account the depot city
			InOut.n_tours += 1;
			ants[0].total_tour_length = Tsp.compute_tour_length(ants[0].tours[i]);
		
			sum += ants[0].total_tour_length;
			ant_empty_memory(ants[0]);
	    }
		
		return sum;
    }*/


}
