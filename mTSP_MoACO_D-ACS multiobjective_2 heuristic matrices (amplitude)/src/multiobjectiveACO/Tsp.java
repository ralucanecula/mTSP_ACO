package multiobjectiveACO;

import java.util.ArrayList;

import multiobjectiveACO.InOut.Distance_type;


public class Tsp {

    static class Point {
		double x;
		double y;
    }

    public static class Problem {
        int n; /* number of cities */
		String name; /* instance name */
		int n_near; /* number of nearest neighbors */
		Point[] nodes; /* array of classes containing coordinates of nodes */
		double[][] distance; /* distance matrix: distance[i][j] gives distance */
		int[][] nn_list; /* nearest neighbor list; contains for each node i a sorted list of n_near nearest neighbors */
    }

    static Problem instance;

    static double dtrunc(double x) {
		int k;
	
		k = (int) x;
		x = (double) k;
		return x;
    }

    /*
     * FUNCTION: the following four functions implement different ways of
     * computing distances for TSPLIB instances
     * INPUT: two node indices
     * OUTPUT: distance between the two nodes
     */
    
    //compute Euclidean distances between two nodes rounded to next integer for TSPLIB instances
    static double euclidianDistance(int i, int j) {
		double xd = 0, yd = 0; 
		
		//we need to compute the distance to the depot city
		if (i == 0 && j != 0) {
			xd = MTsp.depotCity.x - MTsp.instance.nodes[j - 1].x;
			yd = MTsp.depotCity.y - MTsp.instance.nodes[j - 1].y;
		}
		else if (j == 0 && i != 0) {
			xd = MTsp.instance.nodes[i - 1].x - MTsp.depotCity.x;
			yd = MTsp.instance.nodes[i - 1].y - MTsp.depotCity.y;
		}
		else if (i > 0 && j > 0) {
			xd = MTsp.instance.nodes[i - 1].x - MTsp.instance.nodes[j - 1].x;
			yd = MTsp.instance.nodes[i - 1].y - MTsp.instance.nodes[j - 1].y;
		}
		double r = Math.sqrt(xd * xd + yd * yd);
	
		return r;
    }
    
    //compute Euclidean distances between two weight vectors
    static double euclidianDistance(double[] v1, double[] v2) {
    	if (v1.length != v2.length) {
            throw new IllegalArgumentException(
                    "Vectors should be of equal length: x1 length = "
                            + v1.length + " x2 length = " + v2.length);
        }

        double val = 0;
        double temp;
        for (int i = 0; i < v1.length; i++) {
            temp = v1[i] - v2[i];
            temp *= temp;
            val += temp;
        }

        return Math.sqrt(val);
    }

    //compute ceiling distance between two nodes rounded to next integer for TSPLIB instances
    static int ceil_distance(int i, int j) {
		double xd = instance.nodes[i].x - instance.nodes[j].x;
		double yd = instance.nodes[i].y - instance.nodes[j].y;
		double r = Math.sqrt(xd * xd + yd * yd);
	
		return (int) Math.ceil(r);
    }

    //compute geometric distance between two nodes rounded to next integer for TSPLIB instances
    static int geo_distance(int i, int j)
    {
		double deg, min;
		double lati, latj, longi, longj;
		double q1, q2, q3;
		int dd;
		double x1 = instance.nodes[i].x, x2 = instance.nodes[j].x, y1 = instance.nodes[i].y, y2 = instance.nodes[j].y;
	
		deg = dtrunc(x1);
		min = x1 - deg;
		lati = Math.PI * (deg + 5.0 * min / 3.0) / 180.0;
		deg = dtrunc(x2);
		min = x2 - deg;
		latj = Math.PI * (deg + 5.0 * min / 3.0) / 180.0;
	
		deg = dtrunc(y1);
		min = y1 - deg;
		longi = Math.PI * (deg + 5.0 * min / 3.0) / 180.0;
		deg = dtrunc(y2);
		min = y2 - deg;
		longj = Math.PI * (deg + 5.0 * min / 3.0) / 180.0;
	
		q1 = Math.cos(longi - longj);
		q2 = Math.cos(lati - latj);
		q3 = Math.cos(lati + latj);
		dd = (int) (6378.388 * Math.acos(0.5 * ((1.0 + q1) * q2 - (1.0 - q1) * q3)) + 1.0);
		return dd;
    }

    //compute ATT distance between two nodes rounded to next integer for TSPLIB instances
    static int att_distance(int i, int j)
    {
		double xd = instance.nodes[i].x - instance.nodes[j].x;
		double yd = instance.nodes[i].y - instance.nodes[j].y;
		double rij = Math.sqrt((xd * xd + yd * yd) / 10.0);
		double tij = dtrunc(rij);
		int dij;
	
		if (tij < rij)
		    dij = (int) tij + 1;
		else
		    dij = (int) tij;
		return dij;
    }
   
    //computes the matrix of all intercity distances
    static double[][] compute_distances()
    {
		int i, j;
		int size = MTsp.n;
		//include also the depot city in the distance matrix: it will correspond to index 0 for row and column
		double matrix[][] = new double[size + 1][size + 1];
		
		for (i = 0; i < size + 1; i++) {			
		    for (j = 0; j < size + 1; j++) {
				if (InOut.distance_type == Distance_type.ATT) {
				    matrix[i][j] = att_distance(i, j);
				} else if (InOut.distance_type == Distance_type.CEIL_2D) {
				    matrix[i][j] = ceil_distance(i, j);
				} else if (InOut.distance_type == Distance_type.EUC_2D) {
				    matrix[i][j] = euclidianDistance(i, j);
				} else if (InOut.distance_type == Distance_type.GEO) {
				    matrix[i][j] = geo_distance(i, j);
				}
		    }
		}
		return matrix;
    }
    
     //computes the matrix of all intercity distances
    static double[][] computeWeightVectorsDistances(double[][] weightVectors)
    {
		int i, j;
		double matrix[][] = new double[TSP_ACO.N][TSP_ACO.N];
		
		for (i = 0; i < TSP_ACO.N; i++) {			
		    for (j = 0; j < TSP_ACO.N; j++) {	
		    	matrix[i][j] = euclidianDistance(weightVectors[i], weightVectors[j]);
		    }
		}
		return matrix;
    }

    //computes nearest neighbor lists of depth nn for each city
    static int[][] compute_nn_lists() {
		int i, node, nn, count;
	
		int size = MTsp.n;
		double[] distance_vector = new double[size + 1];
		int[] help_vector = new int[size + 1];
	
		nn = Ants.nn_ants;
		if (nn >= size + 1)
		    nn = size - 2;
		Ants.nn_ants = nn;
	
		int[][] m_nnear = new int[size + 1][nn];  //include also the depot city
	
		for (node = 0; node < size + 1; node++) { /* compute candidate-sets for all nodes */	
		    for (i = 0; i < size + 1; i++) { /* Copy distances from nodes to the others */
				distance_vector[i] = MTsp.instance.distance[node][i];
				help_vector[i] = i;
		    }
		    distance_vector[node] = Integer.MAX_VALUE; /* city itself is not nearest neighbor */
		    Utilities.sort2(distance_vector, help_vector, 0, size);
		    //node 0 (which is the depot city) shouldn't be considered in the nearest neighbor list (m_nnear) of a node
		    count = 0; i = -1;
		    while (count < nn) {
		    	i++;
		    	if (help_vector[i] != 0) {
		    		m_nnear[node][count] = help_vector[i];
		    		count++;
		    	}
		    	else {
		    		continue;
		    	}
		    }
		}
	
		return m_nnear;
    }
    
    //computes T nearest weight vectors for each weight vector 
    static int[][] computeNearestWeightVectors(double[][] distance) {
		int i, indexWeightVector, count;
	
		double[] distance_vector = new double[TSP_ACO.N];
		int[] help_vector = new int[TSP_ACO.N];
	
		int[][] m_nnear = new int[TSP_ACO.N][TSP_ACO.T];  
	
		for (indexWeightVector = 0; indexWeightVector < TSP_ACO.N; indexWeightVector++) { /* compute T closest for all weight vectors */	
		    for (i = 0; i < TSP_ACO.N; i++) { /* Copy distances from weight vectors to the others */
				distance_vector[i] = distance[indexWeightVector][i];
				help_vector[i] = i;
		    }
		    distance_vector[indexWeightVector] = 0; 
		    Utilities.sort2(distance_vector, help_vector, 0, TSP_ACO.N - 1);
		    count = 0; i = 0;
		    while (count < TSP_ACO.T) {
	    		m_nnear[indexWeightVector][count] = help_vector[i];
	    		count++;	
	    		i++;
		    }
		}
	
		return m_nnear;
    }
    
    

    //compute the tour length of tour t taking also into account the depot city
   /* static int compute_tour_length(ArrayList<Integer> t) {
		int i;
		int tour_length = 0;
	
		tour_length += MTsp.instance.distance[0][t.get(1) - 1];
		for (i = 1; i < t.size() - 2; i++) {
		    tour_length += MTsp.instance.distance[t.get(i) - 1][t.get(i + 1) - 1];
		}
		tour_length += MTsp.instance.distance[t.get(t.size() - 2) - 1][0];
		
		return tour_length;
    }*/
    
    static double compute_tour_length_(ArrayList<Integer> t) {
        int i;
		double sum = 0;
	
		if (t.size() > 1) {
			sum += MTsp.instance.distance[0][t.get(1) + 1];
			for (i = 1; i < t.size() - 2; i++) {
				sum += MTsp.instance.distance[t.get(i) + 1][t.get(i + 1) + 1];
			}
			sum += MTsp.instance.distance[t.get(t.size() - 2) + 1][0];
		}		
		
		return sum;
    }
      
  //compute the tour length of tour t taking also into account the depot city
    /*static int compute_tour_lengths(ArrayList<Integer>[] t) {
		int i, sum = 0;
		int tour_length = 0;
	
		for (int j = 0; j < t.length; j++) {
			tour_length = 0;
			if (t[j].size() > 1) {
				tour_length += MTsp.instance.distance[0][t[j].get(1) + 1];
				for (i = 1; i < t[j].size() - 2; i++) {
				    tour_length += MTsp.instance.distance[t[j].get(i) + 1][t[j].get(i + 1) + 1];
				}
				tour_length += MTsp.instance.distance[t[j].get(t[j].size() - 2) + 1][0];
				sum += tour_length;
			}		
		}
		return sum;
    }*/

    
}
