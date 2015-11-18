package multiobjectiveACO;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

import multiobjectiveACO.ParetoFront.Solution;

public class Ants {

    static class Ant {
    	//id of the ant
    	int id;
    	//ids/indexes of the subproblems solved by this ant 
    	ArrayList<Integer> idSolvedSubproblems;
    	//denotes which solution from the ones constructed by an ant should be referred next when the ants
    	//from the next subcolony will solve their corresponding subproblem
    	int nextSubproblemId;
    	//for each subproblem this ant solves, for each of the m salesmen an ant will construct a tour, so that a candidate solution constructed by
    	//an ant will be represented by a list of tours, one for each salesman
    	ArrayList<ArrayList<Integer>[]> tours;
    	ArrayList<boolean[]> visited;
    	ArrayList<double[]> tour_lengths;
    	ArrayList<Double> total_tour_lengths;
		//for each subproblem this ant solves, cities left to be visited by an ant (initially toVisit = n, which is the number of cities from the mTSP instance)
    	ArrayList<Integer> toVisit;
		//for each subproblem this ant solves, stores the cost of each solution according to the considered objectives (2 in this case)
    	ArrayList<double[]> costObjectives;
		
		
		public Ant() {}
		
		public Ant(int id_) {
			id = id_;
		}
		
    }

    public static final int MAX_ANTS = 1024;
    public static final int MAX_NEIGHBOURS = 512;

    static Ant ants[];
    static Ant best_so_far_ants[];
    
    //the overlapped subcolonies of ants used to solve subproblems
    static Ant[][] subcolonies;

    static double pheromone[][][];
    static double heuristic[][][];
    //static double total[][][];  //keeps heuristic information times pheromone for each arc

    static double[][] prob_of_selection;
    
    /** the reference point */
    static double[] z;

    static int n_ants; /* number of ants */
 
    static int nn_ants; /* length of nearest neighbor lists for the ants' solution construction */

    static double rho; /* parameter for evaporation */
    static double alpha; /* importance of trail */
    static double beta; /* importance of heuristic evaluate */
    static double q_0; /* probability of best choice in tour construction */

    static boolean as_flag; /* ant system */
    static boolean acs_flag; /* ant colony system (ACS) */

    //for each subproblem there will be a different value for the initial trail
    static double trail_0[] = new double[TSP_ACO.N]; /* initial pheromone level in ACS */
    

    /*static double HEURISTIC(int m, int n) {
    	return (1.0 / (double) MTsp.instance.distance[m][n]);
    }*/

    //allocate the memory for the ant colony, the best-so-far ant
    static void allocate_ants(int[][] closestWeightVectors, ArrayList<Integer>[] indexSubproblems) {
		int i, j, size, indexTour, indexObj, indexAnt;
	
		ants = new Ant[n_ants];
		best_so_far_ants = new Ant[TSP_ACO.N];
		prob_of_selection = new double[TSP_ACO.N][];
	
		for (i = 0; i < n_ants; i++) {
			//how many subproblems the ant i will have to solve
			size = indexSubproblems[i].size();
		    ants[i] = new Ant(i);
		    ants[i].idSolvedSubproblems = new ArrayList<Integer>(size);
		    ants[i].nextSubproblemId = 0;
		    for (j = 0; j < size; j++) {
		    	ants[i].idSolvedSubproblems.add(indexSubproblems[i].get(j));
		    }
		    ants[i].tours = (ArrayList<ArrayList<Integer>[]>)new ArrayList(size);
		    ants[i].tour_lengths = (ArrayList<double[]>)new ArrayList(size);
		    ants[i].total_tour_lengths = new ArrayList<Double>(size);
		    ants[i].visited = (ArrayList<boolean[]>)new ArrayList(size);
		    ants[i].toVisit = new ArrayList<Integer>(size);
		    ants[i].costObjectives = (ArrayList<double[]>)new ArrayList(size);
		    for (j = 0; j < size; j++) {
		    	ants[i].tours.add(j, (ArrayList<Integer>[])new ArrayList[MTsp.m]);
		    	ants[i].tour_lengths.add(j, new double[MTsp.m]);
		    	ants[i].total_tour_lengths.add(j, 0.0);
		    	ants[i].visited.add(j, new boolean[MTsp.n]);
		    	ants[i].toVisit.add(j, MTsp.n);
		    	ants[i].costObjectives.add(j, new double[TSP_ACO.k]);	
		    	
		    	for (indexTour = 0; indexTour < MTsp.m; indexTour++) {
			    	ants[i].tours.get(j)[indexTour] = new ArrayList<Integer>();
			    	ants[i].tour_lengths.get(j)[indexTour] = 0;
		    	}
		    	
		    	for (indexObj = 0; indexObj < TSP_ACO.k; indexObj++) {
			    	ants[i].costObjectives.get(j)[indexObj] = 0;
		    	}	    	
		    }	    
		}
		
		for (i = 0; i < TSP_ACO.N; i++) {
			/*System.out.print("Subproblem " + i + " will be solved by ants: ");*/
			for (j = 0; j < TSP_ACO.T; j++) {
				indexAnt = closestWeightVectors[i][j];
				subcolonies[i][j] = ants[indexAnt];
				/*System.out.print(indexAnt + " ");*/
			}
			/*System.out.println();*/
			
			prob_of_selection[i] = new double[nn_ants + 1];
			for (int i1 = 0; i1 < nn_ants + 1; i1++) {
			    prob_of_selection[i][i1] = Double.POSITIVE_INFINITY;
			}
			
			//initialize best ant for each subproblem/subcolony
			best_so_far_ants[i] = new Ant();		
			best_so_far_ants[i].tours = (ArrayList<ArrayList<Integer>[]>)new ArrayList(1);
			best_so_far_ants[i].tour_lengths = (ArrayList<double[]>)new ArrayList(1);
			best_so_far_ants[i].total_tour_lengths = new ArrayList<Double>(1);
			best_so_far_ants[i].visited = (ArrayList<boolean[]>)new ArrayList(1);
			best_so_far_ants[i].toVisit = new ArrayList<Integer>(1);
			best_so_far_ants[i].costObjectives = (ArrayList<double[]>)new ArrayList(1);
			best_so_far_ants[i].tours.add(0, (ArrayList<Integer>[])new ArrayList[MTsp.m]);
			best_so_far_ants[i].tour_lengths.add(0, new double[MTsp.m]);
			best_so_far_ants[i].total_tour_lengths.add(0, 0.0);
			best_so_far_ants[i].visited.add(0, new boolean[MTsp.n]);
			best_so_far_ants[i].toVisit.add(0, MTsp.n);
			best_so_far_ants[i].costObjectives.add(0, new double[TSP_ACO.k]);	
	    	
	    	for (indexTour = 0; indexTour < MTsp.m; indexTour++) {
	    		best_so_far_ants[i].tours.get(0)[indexTour] = new ArrayList<Integer>();
	    		best_so_far_ants[i].tour_lengths.get(0)[indexTour] = 0;
	    	}
	    	
	    	for (indexObj = 0; indexObj < TSP_ACO.k; indexObj++) {
	    		best_so_far_ants[i].costObjectives.get(0)[indexObj] = 0;
	    	}	
		}

    }
    
    static int find_best(int indexSubproblem) {
    	int indexAnt = 0;
    	double min, value, absValue;
    	double max[] = new double[TSP_ACO.T];
    	
    	min = Double.MAX_VALUE;
    	for (int i = 0; i < TSP_ACO.T; i++) {
    		max[i] = Double.MIN_VALUE;
    	}
    	
		for (int j = 0; j < TSP_ACO.T; j++) {
 			  Ant a = subcolonies[indexSubproblem][j];
 			  for (int nrObj = 0; nrObj < TSP_ACO.k; nrObj++) {
 				  absValue = Math.abs(a.costObjectives.get(a.nextSubproblemId)[nrObj] - z[nrObj]);
 				  value = TSP_ACO.weightVectors[indexSubproblem][nrObj] * absValue;
 				  if (value > max[j]) {
 					 max[j] = value;
 				  }	  
 			  }
 			  a.nextSubproblemId++;
    	}
		for (int j = 0; j < TSP_ACO.T; j++) {
			if (max[j] < min) {
				min = max[j];
				indexAnt = j;
			}
		}
    	
    	return indexAnt;
    }

    //initialize pheromone trails
    //matricea cu urmele de feromoni trebuie sa se faca relativ la toate cele n orase
    static void init_pheromone_trails(int indexSubproblem, double initial_trail)
    {
		int i, j;
		
		/* Initialize pheromone trails */
		for (i = 0; i < MTsp.n + 1; i++) {
		    for (j = 0; j <= i; j++) {
				pheromone[indexSubproblem][i][j] = initial_trail;
				pheromone[indexSubproblem][j][i] = initial_trail;
		    }
		}
		
    }
    
    static void init_heuristic_matrix(int indexSubproblem)
    {
		int i, j;
		double value, sum = 0.0;
		
		/* Initialize matrix with heuristic information corresponding to indexSubproblem */
		for (i = 0; i < MTsp.n + 1; i++) {
		    for (j = 0; j <= i; j++) {
		    	sum = 0.0;
		    	for (int nrObj = 0; nrObj < TSP_ACO.k; nrObj++) {
		    		sum += TSP_ACO.weightVectors[indexSubproblem][nrObj] * MTsp.instance.distance[i][j];
		    	}
		    	value = 1. / sum;
				heuristic[indexSubproblem][i][j] = value;
				heuristic[indexSubproblem][j][i] = value;
		    }
		}
		
    }

    //implements the pheromone trail evaporation
    static void evaporation()
    {
		//for every subproblem/subcolony
		for (int l = 0; l < TSP_ACO.N; l++) {
			 for (int i = 0; i < MTsp.n + 1; i++) {
				    for (int j = 0; j <= i; j++) {
						pheromone[l][i][j] = (1 - rho) * pheromone[l][i][j];
						pheromone[l][j][i] = pheromone[l][i][j];
				    }
			 }
		}
    }

    //reinforces edges used in ant k's solution
    static void global_update_pheromone(int indexSubproblem)
    {
		int i1, i2, j, h, size;
		Map<Integer, Integer>[][] edgesMapArray = (HashMap<Integer, Integer>[][])new HashMap[TSP_ACO.T][MTsp.m];
		boolean belongs = false, alreadyConsidered = false;
		
		for (int l = 0; l < TSP_ACO.T; l++) {
			for (int i = 0; i < MTsp.m; i++) {
				edgesMapArray[l][i] = new HashMap<Integer, Integer>(); 
			}
		}
		
		for (int l = 0; l < TSP_ACO.T; l++) {
			Ant a = Ants.subcolonies[indexSubproblem][l];
			for (i1 = 0; i1 < MTsp.m; i1++) {
				size = a.tours.get(a.nextSubproblemId)[i1].size();
				for (i2 = 0; i2 < size - 1; i2++) {
				    j = a.tours.get(a.nextSubproblemId)[i1].get(i2); 
				    h = a.tours.get(a.nextSubproblemId)[i1].get(i2 + 1);
				    
				    edgesMapArray[l][i1].put(j, h);
				}
			}
			a.nextSubproblemId++;
		}
		
		//reset value of nextSubproblemId field of all the ants
		for (int l = 0; l < Ants.n_ants; l++) {
			Ants.ants[l].nextSubproblemId = 0;
		}
		
		double sum = 0, value, d_tau;
		for (int l = 0; l < TSP_ACO.T; l++) {
			Ant a = Ants.subcolonies[indexSubproblem][l];
			for (i1 = 0; i1 < MTsp.m; i1++) {
				size = a.tours.get(a.nextSubproblemId)[i1].size();
				for (i2 = 0; i2 < size - 1; i2++) {
				    j = a.tours.get(a.nextSubproblemId)[i1].get(i2); 
				    h = a.tours.get(a.nextSubproblemId)[i1].get(i2 + 1);
				    
				    sum = 0; d_tau = 0;
				    alreadyConsidered = false;		    	
				    
				    for (int j2 = 0; j2 < TSP_ACO.T && j2 != l; j2++) {
				    	Ant otherAnt = Ants.subcolonies[indexSubproblem][j2];
				    	//check if edge (j, h) contained in a tour of the second best iteration ant it is 
					    //also used in the tours of the best iteration ant
					    belongs = checkTour(j, h, edgesMapArray[j2]);
					    if (belongs) {
					    	if (j2 < l) {
					    		alreadyConsidered = true;
					    		break;
					    	}
					    	else {
					    		for (int nrObj = 0; nrObj < TSP_ACO.k; nrObj++) {
					  				sum += TSP_ACO.weightVectors[indexSubproblem][nrObj] * otherAnt.costObjectives.get(otherAnt.nextSubproblemId)[nrObj];
					  			}
					    		value = 1. / sum;
					    		d_tau += value;
					    		sum = 0;
					    	}				    		
					    }		    	
				    }
				    if (! alreadyConsidered) {
				    	for (int nrObj = 0; nrObj < TSP_ACO.k; nrObj++) {
			  				sum += TSP_ACO.weightVectors[indexSubproblem][nrObj] * a.costObjectives.get(a.nextSubproblemId)[nrObj];
			  			}
				    	value = 1. / sum;
			    		d_tau += value;
			    		sum = 0;
				    }
				    value = 1. / sum;	
				    pheromone[indexSubproblem][j][h] += d_tau;
				    pheromone[indexSubproblem][h][j] += pheromone[indexSubproblem][j][h];
				}
			}
		}		
		
		/*for (int indexObj = 0; indexObj < TSP_ACO.k; indexObj++) {
		    d_tau = 1.0 / (double) a.costObjectives[indexObj];
			for (i = 0; i < MTsp.m; i++) {
				size = a.tours[i].size();
				for (l = 0; l < size - 1; l++) {
				    j = a.tours[i].get(l); 
				    h = a.tours[i].get(l + 1);
				    
				    j++;
		            h++;
		            
				    pheromone[indexObj][j][h] += d_tau;
				    pheromone[indexObj][h][j] = pheromone[indexObj][j][h];
				}
			} 		
    	}*/
		
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

        if (a.idSolvedSubproblems != null) {
        	int size1 = a.idSolvedSubproblems.size();
        	int size2 = a.total_tour_lengths.size();
        	int size = Math.min(size1, size2);
	        a.nextSubproblemId = 0;
	        for (int index = 0; index < size; index++) {
	        	a.total_tour_lengths.set(index, 0.0);
			    for (int indexObj = 0; indexObj < TSP_ACO.k; indexObj++) {
			    	a.costObjectives.get(index)[indexObj] = 0;
			    }
				
				//clear all the elements (cities) from the m tours of an ant
				for (i = 0; i < MTsp.m; i++) {
				    a.tours.get(index)[i].clear();
				    a.tour_lengths.get(index)[i] = 0;
				}
				
				for (j = 0; j < MTsp.n; j++) {
				    a.visited.get(index)[j] = false;
				}
				a.toVisit.set(index, MTsp.n);
	        }
        }

    }

    //choose for an ant as the next city the one with maximal value of heuristic information times pheromone
    static void choose_best_next(Ant a, int indexSalesman, int indexSubproblem)
    {
		int city, current_city, next_city;
		double value_best;
		double help, heuristicValue = 0, amplitude, sum;
	
		next_city = MTsp.n;
		int lastPos = a.tours.get(a.nextSubproblemId)[indexSalesman].size() - 1;
		current_city = a.tours.get(a.nextSubproblemId)[indexSalesman].get(lastPos);
    	current_city++;
    	
		value_best = -1.; /* values in total matrix are always >= 0.0 */
		for (city = 0; city < MTsp.n; city++) {
			heuristicValue = 0;
			amplitude = 0;
			sum = 0.0;
		    if (a.visited.get(a.nextSubproblemId)[city])
			; /* city already visited, do nothing */
		    else {
		    	//add temporary the node under verification to its corresponding subtour to take it 
		    	//into account when computing the amplitude of the partial constructed subtours
		    	a.tours.get(a.nextSubproblemId)[indexSalesman].add(city);
	
		    	//compute lengths for all partial subtours constructed so far so as to compute 
		    	//their amplitude
		    	for (int i1 = 0; i1 < MTsp.m; i1++) {
		    		a.tour_lengths.get(a.nextSubproblemId)[i1] = Tsp.compute_tour_length_(a.tours.get(a.nextSubproblemId)[i1]);
		    	}
		    	amplitude = computeToursAmplitude(a);
		    	
		    	//remove the node that we've previously added temporary
		    	int lastPos1 = a.tours.get(a.nextSubproblemId)[indexSalesman].size() - 1;
		    	a.tours.get(a.nextSubproblemId)[indexSalesman].remove(lastPos1);
		    	
		    	sum = TSP_ACO.weightVectors[indexSubproblem][0] * MTsp.instance.distance[current_city][city + 1] + 
		    		+ TSP_ACO.weightVectors[indexSubproblem][1] * amplitude;
		    	heuristicValue = 1.0 / sum;
		    	
		    	help = Math.pow(pheromone[indexSubproblem][current_city][city + 1], alpha) * Math.pow(heuristicValue, beta);
				if (help > value_best) {
				    next_city = city;
				    value_best = help;
				}
		    }
		}
		a.tours.get(a.nextSubproblemId)[indexSalesman].add(next_city);
		a.visited.get(a.nextSubproblemId)[next_city] = true;
		int nr = a.toVisit.get(a.nextSubproblemId);
	    a.toVisit.set(a.nextSubproblemId, nr - 1);
    }

    //chooses for an ant as the next city the one with maximal value of heuristic information times pheromone
    static void neighbour_choose_best_next(Ant a, int indexSalesman, int indexSubproblem) {
		int i, current_city, next_city, help_city;
		double value_best, help, heuristicValue = 0, amplitude, sum;
		
		next_city = MTsp.n;   //next_city = Integer.MAX_VALUE;
		int lastPos = a.tours.get(a.nextSubproblemId)[indexSalesman].size() - 1;
		current_city = a.tours.get(a.nextSubproblemId)[indexSalesman].get(lastPos);
        current_city++;
 
		value_best = -1.;  //values in total matrix are always >= 0.0 
		for (i = 0; i < nn_ants; i++) {
			heuristicValue = 0;
			amplitude = 0;
			sum = 0.0;
		    help_city = MTsp.instance.nn_list[current_city][i];
		    if (a.visited.get(a.nextSubproblemId)[help_city - 1])
		    	; // city already visited, do nothing 
		    else {
		    	//add temporary the node under verification to its corresponding subtour to take it 
		    	//into account when computing the amplitude of the partial constructed subtours
		    	a.tours.get(a.nextSubproblemId)[indexSalesman].add(help_city - 1);
	
		    	//compute lengths for all partial subtours constructed so far so as to compute 
		    	//their amplitude
		    	for (int i1 = 0; i1 < MTsp.m; i1++) {
		    		a.tour_lengths.get(a.nextSubproblemId)[i1] = Tsp.compute_tour_length_(a.tours.get(a.nextSubproblemId)[i1]);
		    	}
		    	amplitude = computeToursAmplitude(a);
		    	
		    	//remove the node that we've previously added temporary
		    	int lastPos1 = a.tours.get(a.nextSubproblemId)[indexSalesman].size() - 1;
		    	a.tours.get(a.nextSubproblemId)[indexSalesman].remove(lastPos1);
		    	
		    	sum = TSP_ACO.weightVectors[indexSubproblem][0] * MTsp.instance.distance[current_city][help_city] + 
		    		+ TSP_ACO.weightVectors[indexSubproblem][1] * amplitude;
		    	heuristicValue = 1.0 / sum;
		    	
		    	help = Math.pow(pheromone[indexSubproblem][current_city][help_city], alpha) * Math.pow(heuristicValue, beta);
				//help = total[current_city][help_city];
				if (help > value_best) {
				    value_best = help;
				    next_city = help_city - 1;
				}
		    }
		}
		
		if (next_city == MTsp.n)
		    // all cities in nearest neighbor list were already visited 
		    choose_best_next(a, indexSalesman, indexSubproblem);
		else {
			a.tours.get(a.nextSubproblemId)[indexSalesman].add(next_city);
		    a.visited.get(a.nextSubproblemId)[next_city] = true;
		    int nr = a.toVisit.get(a.nextSubproblemId);
		    a.toVisit.set(a.nextSubproblemId, nr - 1);
		}
    }
    
    static void choose_closest_next(Ant a, int indexSalesman)
    {
		int current_city, next_city, city;
	    double  min_distance;
		
		next_city = MTsp.n;
		int lastPos = a.tours.get(0)[indexSalesman].size() - 1;
		current_city = a.tours.get(0)[indexSalesman].get(lastPos);
		current_city++;
		
		min_distance = Integer.MAX_VALUE;  //Search shortest edge 
		for (city = 0; city < MTsp.n; city++) {
		    if (a.visited.get(0)[city])
			;  //city already visited 
		    else {
				if (MTsp.instance.distance[current_city][city + 1] < min_distance) {
				    next_city = city;
				    min_distance = MTsp.instance.distance[current_city][city + 1];
				}
		    }
		}
		a.tours.get(0)[indexSalesman].add(next_city);
		a.visited.get(0)[next_city] = true; 
		int nr = a.toVisit.get(0);
		nr--;
		a.toVisit.set(0, nr);
    }

    //Choose for an ant probabilistically a next city among all unvisited cities in the current city's candidate list
    static void neighbour_choose_and_move_to_next(Ant a, int indexSalesman, int indexSubproblem) {
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
		    neighbour_choose_best_next(a, indexSalesman, indexSubproblem);
		    return;
		}
	
		prob_ptr = prob_of_selection[indexSubproblem];
	
	    /* current_city city of ant k */
		int lastPos = a.tours.get(a.nextSubproblemId)[indexSalesman].size() - 1;
		current_city = a.tours.get(a.nextSubproblemId)[indexSalesman].get(lastPos);
    	current_city++;
    	double heuristicValue = 0, amplitude, sum;
    
		for (i = 0; i < nn_ants; i++) {
			heuristicValue = 0;
			amplitude = 0;
			sum = 0.0;
			city = MTsp.instance.nn_list[current_city][i];
		    if (a.visited.get(a.nextSubproblemId)[city - 1])
		    	prob_ptr[i] = 0.0; /* city already visited */
		    else {
		    	//add temporary the node under verification to its corresponding subtour to take it 
		    	//into account when computing the amplitude of the partial constructed subtours
		    	a.tours.get(a.nextSubproblemId)[indexSalesman].add(city - 1);
	
		    	//compute lengths for all partial subtours constructed so far so as to compute 
		    	//their amplitude
		    	for (int i1 = 0; i1 < MTsp.m; i1++) {
		    		a.tour_lengths.get(a.nextSubproblemId)[i1] = Tsp.compute_tour_length_(a.tours.get(a.nextSubproblemId)[i1]);
		    	}
		    	amplitude = computeToursAmplitude(a);
		    	
		    	//remove the node that we've previously added temporary
		    	int lastPos1 = a.tours.get(a.nextSubproblemId)[indexSalesman].size() - 1;
		    	a.tours.get(a.nextSubproblemId)[indexSalesman].remove(lastPos1);
		    	
		    	sum = TSP_ACO.weightVectors[indexSubproblem][0] * MTsp.instance.distance[current_city][city] + 
		    		+ TSP_ACO.weightVectors[indexSubproblem][1] * amplitude;
		    	heuristicValue = 1.0 / sum;
		    	
				prob_ptr[i] = Math.pow(pheromone[indexSubproblem][current_city][city], alpha) * Math.pow(heuristicValue, beta);
				sum_prob += prob_ptr[i];
		    }
		}
	
		if (sum_prob <= 0.0) {
		    /* All cities from the candidate set are tabu (are already visited) */
		    choose_best_next(a, indexSalesman, indexSubproblem);
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
				if (i < prob_ptr.length) {
					partial_sum += prob_ptr[i];
				}
				else {
					System.out.println("Not reaching the sum; i=" + i);
					break;
				}
		    }
		    /*
		     * This may very rarely happen because of rounding if rnd is close to 1.
		     */
		    if (i == nn_ants || i == (nn_ants + 1)) {
		    	//System.out.println("Rare situation: got here..");
		    	neighbour_choose_best_next(a, indexSalesman, indexSubproblem);
				return;
		    }
		    help = MTsp.instance.nn_list[current_city][i];
		    a.tours.get(a.nextSubproblemId)[indexSalesman].add(help - 1);
		    a.visited.get(a.nextSubproblemId)[help - 1] = true;
			int nr = a.toVisit.get(a.nextSubproblemId);
		    a.toVisit.set(a.nextSubproblemId, nr - 1);
		}
    }

    //reinforces the edges used in ant's solution as in ACS
    //global pheromone update performed by the best iteration ant
    static void global_acs_pheromone_update(int indexSubproblem) {
		int i, j, h, l, size;
		double d_tau, sum = 0;
		
		Ant bestAnt = Ants.best_so_far_ants[indexSubproblem];
		for (int nrObj = 0; nrObj < TSP_ACO.k; nrObj++) {
			sum += TSP_ACO.weightVectors[indexSubproblem][nrObj] * bestAnt.costObjectives.get(0)[nrObj];
		}
		d_tau = 1.0 / sum;
	
		for (i = 0; i < MTsp.m; i++) {
			size = bestAnt.tours.get(0)[i].size();
			for (l = 0; l < size - 1; l++)  {
			    j = bestAnt.tours.get(0)[i].get(l);
			    h = bestAnt.tours.get(0)[i].get(l + 1);	
			    
			    j++;
			    h++;
		
			    pheromone[indexSubproblem][j][h] = (1. - rho) * pheromone[indexSubproblem][j][h] + rho * d_tau;
			    pheromone[indexSubproblem][h][j] = pheromone[indexSubproblem][j][h];
			}
		}
	
    }
    
    static void global_acs_share_pheromone_trails(int indexSubproblem) {
    	int i, j, h, l, size, idSuproblem;
		double d_tau, sum = 0;
		
		Ant bestAnt = Ants.best_so_far_ants[indexSubproblem];
		ArrayList<Integer> idSubproblems = bestAnt.idSolvedSubproblems;
		
		for (int i1 = 0; i1 < idSubproblems.size(); i1++) {
			idSuproblem = idSubproblems.get(i1);
			
			if (idSuproblem != indexSubproblem) {
				for (int nrObj = 0; nrObj < TSP_ACO.k; nrObj++) {
					sum += TSP_ACO.weightVectors[idSuproblem][nrObj] * bestAnt.costObjectives.get(0)[nrObj];
				}
				d_tau = 1. / sum;
				
				for (i = 0; i < MTsp.m; i++) {
					size = bestAnt.tours.get(0)[i].size();
					for (l = 0; l < size - 1; l++)  {
					    j = bestAnt.tours.get(0)[i].get(l);
					    h = bestAnt.tours.get(0)[i].get(l + 1);	
					    
					    j++;
					    h++;
				
					    pheromone[idSuproblem][j][h] = pheromone[idSuproblem][j][h] + d_tau;
					    pheromone[idSuproblem][h][j] = pheromone[idSuproblem][j][h];
					}
				}
				sum = 0;
				d_tau = 0;
			}
		}
		
    }
    
    //check if edge (i, j) is used in the tours of the best iteration ant
    static boolean checkTour(int startNode, int endNode, Map<Integer, Integer>[] edgesMapArray) {
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

    //removes some pheromone on edge just passed by the ant
    static void local_acs_pheromone_update(Ant a, int indexSalesman, int indexSubproblem) {
		int h, j;
		
		int lastPos = a.tours.get(a.nextSubproblemId)[indexSalesman].size() - 1;
		j = a.tours.get(a.nextSubproblemId)[indexSalesman].get(lastPos);
		h = a.tours.get(a.nextSubproblemId)[indexSalesman].get(lastPos - 1);	
		
		j++;
		h++;

		/* still additional parameter has to be introduced */
		pheromone[indexSubproblem][h][j] = (1. - 0.1) * pheromone[indexSubproblem][h][j] + 0.1 * trail_0[indexSubproblem];
		pheromone[indexSubproblem][j][h] = pheromone[indexSubproblem][h][j];	
    }
    
    //create a copy of the ant a and return the created copy at the output
    static Solution copyAnt(Ant a) {
    	//first create an empty solution
    	Solution copy = new Solution(a.id);
    	copy.tours = (ArrayList<Integer>[])new ArrayList[MTsp.m];
    	copy.tour_lengths = new double[MTsp.m];
	    for (int j = 0; j < MTsp.m; j++) {
	    	copy.tours[j] = new ArrayList<Integer>();
	    	copy.tour_lengths[j] = 0;
	    }
	   
	    copy.costObjectives = new double[TSP_ACO.k];
 	
	    //then copy the information from the ant a
		for (int indexObj = 0; indexObj < TSP_ACO.k; indexObj++) {
			copy.costObjectives[indexObj] = a.costObjectives.get(a.nextSubproblemId)[indexObj];
    	}
		for (int i = 0; i < MTsp.m; i++) {
			copy.tour_lengths[i] = a.tour_lengths.get(a.nextSubproblemId)[i];
			int size = a.tours.get(a.nextSubproblemId)[i].size();
			copy.tours[i] = new ArrayList<Integer>(size);
			for (int j = 0; j < size; j++) {
				int elem = a.tours.get(a.nextSubproblemId)[i].get(j);
				copy.tours[i].add(elem);
			}		
		}
    	
    	return copy;
    }

    //copy solution from ant a1 into ant a2
    static void copy_from_to(Ant a1, Ant a2) {
		int i, j;
	
		Ants.ant_empty_memory(a2);
		a2.id = a1.id;
		a2.nextSubproblemId = a1.nextSubproblemId;
		int size = a1.idSolvedSubproblems.size();
		a2.idSolvedSubproblems = new ArrayList<Integer>(size);
		for (int index = 0; index < size; index++) {
			a2.idSolvedSubproblems.add(a1.idSolvedSubproblems.get(index));
		}
		
		int index = a1.nextSubproblemId - 1;
		a2.total_tour_lengths.set(0, a1.total_tour_lengths.get(index));
		for (int indexObj = 0; indexObj < TSP_ACO.k; indexObj++) {
		    a2.costObjectives.get(0)[indexObj] = a1.costObjectives.get(index)[indexObj];
		}
		for (i = 0; i < MTsp.m; i++) {
			a2.tour_lengths.get(0)[i] = a1.tour_lengths.get(index)[i];
			int size1 = a1.tours.get(index)[i].size();
			a2.tours.get(0)[i] = new ArrayList<Integer>(size1);
			for (j = 0; j < size1; j++) {
				int elem = a1.tours.get(index)[i].get(j);
				a2.tours.get(0)[i].add(elem);
			}	
		}
		
    }
    
    static double computeToursAmplitude(Ant a) {
    	double min, max;
		int i;

		min = a.tour_lengths.get(a.nextSubproblemId)[0];
		max = a.tour_lengths.get(a.nextSubproblemId)[0];
		for (i = 1; i < a.tours.get(a.nextSubproblemId).length; i++) {
		    if (a.tour_lengths.get(a.nextSubproblemId)[i] < min) {
				min = a.tour_lengths.get(a.nextSubproblemId)[i];
		    }
		    if (a.tour_lengths.get(a.nextSubproblemId)[i] > max) {
				max = a.tour_lengths.get(a.nextSubproblemId)[i];
		    }	    
		}
		
		return (max - min);
    }

    //generate some nearest neighbor tour and return the cost of this tour according to each objective
    static double[] nn_tour() {
    	int step, salesman;
    	double sum = 0;
    	double amplitude = 0;
    	double[] costValues = new double[2];
    	
    	ant_empty_memory(ants[0]);
    	step = 0;
    	
    	for (int i = 0; i < MTsp.m; i++) {
    		//place the ant on the depot city, which is the start city of each tour
			// -1 is a special marker for the deport city, so that it's not be confused with the rest of the cities
			// all the rest of the cities are represented by integer values > 0
    		ants[0].tours.get(0)[i].add(-1);  	
    	}
		
		while (ants[0].toVisit.get(0) > 0) {   //there are still left cities to be visited	
    		//choose in a random manner with equal probability which salesman to consider next, that will visit a city
	    	salesman = (int)(Math.random() * MTsp.m);
	    	choose_closest_next(ants[0], salesman);	 
		}
		
		for (int i = 0; i < MTsp.m; i++) {
			step = ants[0].tours.get(0)[i].size();
			ants[0].tours.get(0)[i].add(step, -1);
			
			ants[0].tour_lengths.get(0)[i] = Tsp.compute_tour_length_(ants[0].tours.get(0)[i]);
			double nr = ants[0].total_tour_lengths.get(0);
			ants[0].total_tour_lengths.set (0, nr + ants[0].tour_lengths.get(0)[i]);
		}
		sum = ants[0].total_tour_lengths.get(0);
		amplitude = computeToursAmplitude(ants[0]);
		costValues[0] = sum;
		costValues[1] = amplitude;
		
		ant_empty_memory(ants[0]);
		
		return costValues;
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
