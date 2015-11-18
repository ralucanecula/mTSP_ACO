package multiobjectiveACO;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Random;

import multiobjectiveACO.Ants.Ant;

import org.apache.poi.ss.usermodel.Cell;
import org.apache.poi.ss.usermodel.CellStyle;
import org.apache.poi.ss.usermodel.Font;
import org.apache.poi.ss.usermodel.IndexedColors;
import org.apache.poi.ss.usermodel.Row;
import org.apache.poi.xssf.usermodel.XSSFSheet;
import org.apache.poi.xssf.usermodel.XSSFWorkbook;

public class Utilities {

    private static Random random;

    static int seed;
    static int rowIndex = 1;
    
    //private static String filePath = "output/Experiments Fuzzy cMeans.xlsx";
    
    private static String filePath = "../../../Conferinte/CEC 2015/Rulari ACO + k-means_" + TSP_ACO.instanceName + ".xlsx";
    //private static String filePath1 = "output/ParetoFront_" + TSP_ACO.instanceName + " (m=" + MTsp.m + ")" + ".xlsx";
    private static String filePath1 = "../../../Conferinte/GECCO_2015/multiobjective ACO/MACS (2 matrici cu informatie euristica)/ParetoFront_" + TSP_ACO.instanceName + " (m=" + MTsp.m + ")_amplitude_2200 iter (MACS)_2 heuristic" + ".xlsx";
    private static String filePath2 = "../../../Conferinte/GECCO_2015/fisiere intermediare/Rulari MACS_2 heuristic_" + TSP_ACO.instanceName + " (m=" + MTsp.m + ")_total, subtour costs and amplitude.txt";
    
    //for setting the content to be inserted in the Excel file
    private static String tours[];
    private static int nrCities[];
    private static double subtoursCost[];
    private static double totalCost;
    
    //for verbose output: save at each 5 iteration the best (minimum) total cost of all m subtours 
    private static ArrayList<Double> iterTotalCost;

    //auxiliary routine for sorting an integer array
    static void swap2(double v[], int v2[], int i, int j) {
    	double tmp1;
		int tmp2;
	
		tmp1 = v[i];
		v[i] = v[j];
		v[j] = tmp1;
		
		tmp2 = v2[i];
		v2[i] = v2[j];
		v2[j] = tmp2;
    }

    //recursive routine (quicksort) for sorting one array; second array does the same sequence of swaps
    static void sort2(double v[], int v2[], int left, int right) {
		    int k, last;
		
			if (left >= right)
			    return;
			swap2(v, v2, left, (left + right) / 2);
			last = left;
			for (k = left + 1; k <= right; k++)
			    if (v[k] < v[left])
			    	swap2(v, v2, ++last, k);
			swap2(v, v2, left, last);
			sort2(v, v2, left, last);
			sort2(v, v2, last + 1, right);
    }
    
    //auxiliary routine for sorting an integer array and a double array
    static void swap2_(Double v[], Integer v2[], int i, int j) {
    	double tmp1;
		int tmp2;
	
		tmp1 = v[i];
		v[i] = v[j];
		v[j] = tmp1;
		
		tmp2 = v2[i];
		v2[i] = v2[j];
		v2[j] = tmp2;
    }

    //recursive routine (quicksort) for sorting one array; second array does the same sequence of swaps
    static void sort2_(Double v[], Integer v2[], int left, int right) {
		    int k, last;
		
			if (left >= right)
			    return;
			swap2_(v, v2, left, (left + right) / 2);
			last = left;
			for (k = left + 1; k <= right; k++)
			    if (v[k] < v[left])
			    	swap2_(v, v2, ++last, k);
			swap2_(v, v2, left, last);
			sort2_(v, v2, left, last);
			sort2_(v, v2, last + 1, right);
    }

    //generate a random number that is uniformly distributed in [0,1]
    static double random01() {
		if (random == null) {
		    random = new Random();
		}
	
		return random.nextDouble();
    }
    
    static void writeExcel(int n, int m, int result) {
    	//the file already exists; we should add a new row as the last one in the Excel file
	    if (new File(filePath).canRead()) {
	    	//System.out.println("File already exists..");
	    	try {
		    	FileInputStream file = new FileInputStream(new File(filePath));
		    	
		    	//Create Workbook instance holding reference to .xlsx file
	            XSSFWorkbook workbook1 = new XSSFWorkbook(file);
	 
	            //Get first/desired sheet from the workbook
	            XSSFSheet sheet1 = workbook1.getSheetAt(2);
		    	int countRows = sheet1.getLastRowNum() + 1;
		    	Row newRow = sheet1.createRow(countRows++);
		    	
		    	int cellnum = 0;
		    	Cell cell = newRow.createCell(cellnum++);
		    	cell.setCellValue(n);
		    	cell = newRow.createCell(cellnum++);
		    	cell.setCellValue(m);
		    	cell = newRow.createCell(cellnum++);
		    	cell.setCellValue(result);
			    
				//Write the workbook in file system
			    FileOutputStream out = new FileOutputStream(new File(filePath));
			    workbook1.write(out);
			    out.close();
			    
			    //System.out.println("Written successfully on disk.");
	    	}
			catch (Exception e) {
			    e.printStackTrace();
			}

	    }
	    else {
	    	//Blank workbook
			XSSFWorkbook workbook2 = new XSSFWorkbook(); 
			
			//Create a blank sheet
			XSSFSheet sheet2 = workbook2.createSheet("Results - 51 cities");
		 
			//Iterate over data and write to sheet
			int rownum = 0, cellnum = 0;
			Row row = sheet2.createRow(rownum++);
			Cell cell = row.createCell(cellnum++);
			cell.setCellValue(n);
	    	cell = row.createCell(cellnum++);
	    	cell.setCellValue(m);
	    	cell = row.createCell(cellnum++);
	    	cell.setCellValue(result);
			
			try {
				//Write the workbook in file system
			    FileOutputStream out = new FileOutputStream(new File(filePath));
			    workbook2.write(out);
			    out.close();
			    
			    //System.out.println("Written successfully on disk.");
			} 
			catch (Exception e) {
			    e.printStackTrace();
			}
			    
	    }
    }
    
    static void writeResultsExcel(int trialNumber, boolean saveIterCosts) {
    	Row r, r1;
    	Cell c;
    	int index1 = 0, index2 = 0, index3 = 0, index4 = 0, index5 = 0;
    	
    	//the file already exists; we should add a new row as the last one in the Excel file
	    if (new File(filePath).canRead()) {
	    	//System.out.println("File already exists..");
	    	try {
		    	FileInputStream file = new FileInputStream(new File(filePath));
		    	
		    	//Create Workbook instance holding reference to .xlsx file
	            XSSFWorkbook workbook1 = new XSSFWorkbook(file);
	            
	            int startIndex = 0, rowIndex = 0;
	            switch (MTsp.m) {
	            	case 2: 
	            		startIndex = 0;
	            		rowIndex = 4;
	            		break;
	            	case 3: 
	            		startIndex = 5;
	            		rowIndex = 5;
	            		break;
	            	case 5: 
	            		startIndex = 10;
	            		rowIndex = 7;
	            		break;
	            	case 7: 
	            		startIndex = 15;
	            		rowIndex = 9;
	            		break;
	            	default:
	            		System.out.println("Unknown value for m");
	            		break;          		
	            }
	            
	            //Get desired sheet from the workbook
	            XSSFSheet sheet1 = workbook1.getSheetAt(startIndex);  //for tours
	            XSSFSheet sheet2 = workbook1.getSheetAt(startIndex + 1);  //for number of assigned cities
	            XSSFSheet sheet3 = workbook1.getSheetAt(startIndex + 2);  //for cost of individual subtours
	            XSSFSheet sheet4 = workbook1.getSheetAt(startIndex + 3);  //for total cost of subtours
	            XSSFSheet sheet5 = workbook1.getSheetAt(startIndex + 4);  //for verbose output of individual costs at each 5 iteration
	            
	            //define a cell style for bold font
	            CellStyle style = workbook1.createCellStyle();
	            Font font = workbook1.createFont();
	            font.setBoldweight(Font.BOLDWEIGHT_BOLD);
	            style.setFont(font);
	            
	            //define style with bold font and blue color for font
	            CellStyle styleBoldBlue = workbook1.createCellStyle();
	            font = workbook1.createFont();
	            font.setBoldweight(Font.BOLDWEIGHT_BOLD);
	            font.setColor(IndexedColors.BLUE.index);
	            styleBoldBlue.setFont(font);
	            
	            index1 = 467;
	            if (!saveIterCosts) {
	            	//write only once the name of the algorithm that was run
	            	if (trialNumber == 0) {
			            r = sheet1.getRow(index1); 
			            if (r == null) {
			               // First cell in the row, create
			               //System.out.println("Empty row, create new one");
			               r = sheet1.createRow(index1);
			            }
		
			            c = r.getCell(0); 
			            if (c == null) {
			                // New cell
			            	//System.out.println("Empty cell, create new one");
			                c = r.createCell(0);
			            }
			            c.setCellValue("Obtained tours after running vers. II (ACO global, fara clustering)");
			            c.setCellStyle(styleBoldBlue);
	            	}
	            	
	            	//write number of run
		            index1 = index1 + 2 + rowIndex * trialNumber;
		            r = sheet1.getRow(index1); 
		            if (r == null) {
		               // First cell in the row, create
		               //System.out.println("Empty row, create new one");
		               r = sheet1.createRow(index1);
		            }
	
		            c = r.getCell(0); 
		            if (c == null) {
		                // New cell
		            	//System.out.println("Empty cell, create new one");
		                c = r.createCell(0);
		            }
		            c.setCellValue("Run #" + (trialNumber + 1));
		            c.setCellStyle(style);
		            
		            //write obtained subtours
		            index1++;
		            String toursText[] = getTours();
		            for (int i = 0; i < toursText.length; i++) {
		                r = sheet1.getRow(index1); 
			            if (r == null) {
			               // First cell in the row, create
			               //System.out.println("Empty row, create new one");
			               r = sheet1.createRow(index1);
			            }
		
			            c = r.getCell(0); 
			            if (c == null) {
			                // New cell
			            	//System.out.println("Empty cell, create new one");
			                c = r.createCell(0);
			            }
			            c.setCellValue(toursText[i]);
			            index1++;
		            }
		            
		        //write only once the name of the algorithm that was run
		        index2 = 57;
            	if (trialNumber == 0) {
		            r = sheet2.getRow(index2); 
		            if (r == null) {
		               // First cell in the row, create
		               //System.out.println("Empty row, create new one");
		               r = sheet2.createRow(index2);
		            }
	
		            c = r.getCell(0); 
		            if (c == null) {
		                // New cell
		            	//System.out.println("Empty cell, create new one");
		                c = r.createCell(0);
		            }
		            c.setCellValue("Number of cities assigned in each tour after running vers. II (ACO global, fara clustering)");
		            c.setCellStyle(styleBoldBlue);
		            
		            //write only once the table header
		            index2 = index2 + 3;
		            r = sheet2.getRow(index2); 
		            if (r == null) {
		               // First cell in the row, create
		               //System.out.println("Empty row, create new one");
		               r = sheet2.createRow(index2);
		            }
	
		            c = r.getCell(0); 
		            if (c == null) {
		                // New cell
		            	//System.out.println("Empty cell, create new one");
		                c = r.createCell(0);
		            }
		            c.setCellValue("Run #");
		            c.setCellStyle(style);
		            
		            c = r.getCell(1); 
		            if (c == null) {
		                // New cell
		            	//System.out.println("Empty cell, create new one");
		                c = r.createCell(1);
		            }
		            c.setCellValue("List with number of assigned cities");
		            c.setCellStyle(style);
		            
            	}
      
		            //write number of run
		            index2 = 61 + trialNumber;
		            r = sheet2.getRow(index2); 
		            if (r == null) {
		               // First cell in the row, create
		               //System.out.println("Empty row, create new one");
		               r = sheet2.createRow(index2);
		            }
	
		            c = r.getCell(0); 
		            if (c == null) {
		                // New cell
		            	//System.out.println("Empty cell, create new one");
		                c = r.createCell(0);
		            }
		            c.setCellValue(trialNumber + 1);
		            
		            //write number of assigned cities in each subtour
		            int nrCities[] = getNrCities();
		            for (int i = 0; i < nrCities.length; i++) {
			            c = r.getCell(i + 1); 
			            if (c == null) {
			                // New cell
			            	//System.out.println("Empty cell, create new one");
			                c = r.createCell(i + 1);
			            }
			            c.setCellValue(nrCities[i]);
		            }
		            
		            index3 = 57;
		            //write only once the name of the algorithm that was run
	            	if (trialNumber == 0) {		            
			            r = sheet3.getRow(index3); 
			            if (r == null) {
			               // First cell in the row, create
			               //System.out.println("Empty row, create new one");
			               r = sheet3.createRow(index3);
			            }
		
			            c = r.getCell(0); 
			            if (c == null) {
			                // New cell
			            	//System.out.println("Empty cell, create new one");
			                c = r.createCell(0);
			            }
			            c.setCellValue("Cost of each subtour after running vers. II (ACO global, fara clustering)");
			            c.setCellStyle(styleBoldBlue);
			            
			            //write only once the table header
			            index3 = index3 + 3;
			            r = sheet3.getRow(index3); 
			            if (r == null) {
			               // First cell in the row, create
			               //System.out.println("Empty row, create new one");
			               r = sheet3.createRow(index3);
			            }
		
			            c = r.getCell(0); 
			            if (c == null) {
			                // New cell
			            	//System.out.println("Empty cell, create new one");
			                c = r.createCell(0);
			            }
			            c.setCellValue("Run #");
			            c.setCellStyle(style);
			            
			            c = r.getCell(1); 
			            if (c == null) {
			                // New cell
			            	//System.out.println("Empty cell, create new one");
			                c = r.createCell(1);
			            }
			            c.setCellValue("List with cost of subtours");
			            c.setCellStyle(style);  
	            	}
		            
		            //write number of run
		            index3 = 61 + trialNumber;
		            r = sheet3.getRow(index3); 
		            if (r == null) {
		               // First cell in the row, create
		               //System.out.println("Empty row, create new one");
		               r = sheet3.createRow(index3);
		            }
	
		            c = r.getCell(0); 
		            if (c == null) {
		                // New cell
		            	//System.out.println("Empty cell, create new one");
		                c = r.createCell(0);
		            }
		            c.setCellValue(trialNumber + 1);
		            
		            //write cost of each individual subtour
		            double subtoursCost[] = getSubtoursCost();
		            for (int i = 0; i < subtoursCost.length; i++) {
			            c = r.getCell(i + 1); 
			            if (c == null) {
			                // New cell
			            	//System.out.println("Empty cell, create new one");
			                c = r.createCell(i + 1);
			            }
			            c.setCellValue(subtoursCost[i]);
		            }
		            
		          index4 = 57;  
		          //write only once the name of the algorithm that was run
	            	if (trialNumber == 0) {			            
			            r = sheet4.getRow(index4); 
			            if (r == null) {
			               // First cell in the row, create
			               //System.out.println("Empty row, create new one");
			               r = sheet4.createRow(index4);
			            }
		
			            c = r.getCell(0); 
			            if (c == null) {
			                // New cell
			            	//System.out.println("Empty cell, create new one");
			                c = r.createCell(0);
			            }
			            c.setCellValue("Total cost of all subtours after running vers. II (ACO global, fara clustering)");
			            c.setCellStyle(styleBoldBlue);
			            
			            //write only once the table header
			            index4 = index4 + 3;
			            r = sheet4.getRow(index4); 
			            if (r == null) {
			               // First cell in the row, create
			               //System.out.println("Empty row, create new one");
			               r = sheet4.createRow(index4);
			            }
		
			            c = r.getCell(0); 
			            if (c == null) {
			                // New cell
			            	//System.out.println("Empty cell, create new one");
			                c = r.createCell(0);
			            }
			            c.setCellValue("Run #");
			            c.setCellStyle(style);
			            
			            c = r.getCell(1); 
			            if (c == null) {
			                // New cell
			            	//System.out.println("Empty cell, create new one");
			                c = r.createCell(1);
			            }
			            c.setCellValue("Total cost of subtours");
			            c.setCellStyle(style); 
	            	}
		            
		            //write number of run
		            index4 = 61 + trialNumber;
		            r = sheet4.getRow(index4); 
		            if (r == null) {
		               // First cell in the row, create
		               //System.out.println("Empty row, create new one");
		               r = sheet4.createRow(index4);
		            }
	
		            c = r.getCell(0); 
		            if (c == null) {
		                // New cell
		            	//System.out.println("Empty cell, create new one");
		                c = r.createCell(0);
		            }
		            c.setCellValue(trialNumber + 1);
		            
		            //write total cost of all subtours
		            double totalCost = getTotalCost();
		            c = r.getCell(1); 
		            if (c == null) {
		                // New cell
		            	//System.out.println("Empty cell, create new one");
		                c = r.createCell(1);
		            }
		            c.setCellValue(totalCost);
	            }
	            
	            index5 = 687;
	            if (saveIterCosts) {
	            	//write only once the name of the algorithm that was run
	            	if (trialNumber == 0) {			            
			            r = sheet5.getRow(index5); 
			            if (r == null) {
			               // First cell in the row, create
			               //System.out.println("Empty row, create new one");
			               r = sheet5.createRow(index5);
			            }
		
			            c = r.getCell(0); 
			            if (c == null) {
			                // New cell
			            	//System.out.println("Empty cell, create new one");
			                c = r.createCell(0);
			            }
			            c.setCellValue("Total cost of subtours at each 5 iteration after running vers. II (ACO global, fara clustering)");
			            c.setCellStyle(styleBoldBlue);
	            	}
	            	
	            	index5 = index5 + 3;
		            r = sheet5.getRow(index5); 
		            if (r == null) {
		               // First cell in the row, create
		               //System.out.println("Empty row, create new one");
		               r = sheet5.createRow(index5);
		            }
		            
		            //for each trial run save at each 5 iteration the best total cost so far
		            ArrayList<Double> iterTotalCost = getIterTotalCost();
		            int index;
		           
	            	//for each run write the table header cell
		            c = r.getCell(trialNumber); 
		            if (c == null) {
		                // New cell
		            	//System.out.println("Empty cell, create new one");
		                c = r.createCell(trialNumber);
		            }
		            c.setCellValue("Run " + (trialNumber + 1));
		            c.setCellStyle(style); 
	            	
	            	for (int j = 0; j < iterTotalCost.size(); j++) {
	            		index = index5 + 1 + j;
			            r1 = sheet5.getRow(index); 
			            if (r1 == null) {
			               // First cell in the row, create
			               //System.out.println("Empty row, create new one");
			               r1 = sheet5.createRow(index);
			            }
		
			            c = r1.getCell(trialNumber); 
			            if (c == null) {
			                // New cell
			            	//System.out.println("Empty cell, create new one");
			                c = r1.createCell(trialNumber);
			            }
			            c.setCellValue(iterTotalCost.get(j));     
	            	}	
		            
	            }
	           
				//Write the workbook in file system
			    FileOutputStream out = new FileOutputStream(new File(filePath));
			    workbook1.write(out);
			    out.close();
			    
			    int nrOfRun = trialNumber + 1;
			    System.out.println("\nRun #" + nrOfRun + " written successfully on disk.\n");
	    	}
			catch (Exception e) {
			    e.printStackTrace();
			}

	    }
	    else {
	    	//Blank workbook
	    	System.out.println("File " + filePath + " doesn't exists.."); 
			    
	    }
    }
    
    static void writeParetoSet(ArrayList<Ant> bestSoFarPareto, int trial) {
    	Row r;
    	Cell c;
    	int lineNumber = 0;
    	
    	//the file already exists; we should add a new row as the last one in the Excel file
	    if (new File(filePath1).canRead()) {
	    	//System.out.println("File already exists..");
	    	try {
		    	FileInputStream file = new FileInputStream(new File(filePath1));
		    	
		    	//Create Workbook instance holding reference to .xlsx file
	            XSSFWorkbook workbook1 = new XSSFWorkbook(file);
	 
	            //Get first/desired sheet from the workbook
	            XSSFSheet sheet1 = workbook1.getSheetAt(trial);
	            
	            //write table header cells
	            r = sheet1.getRow(lineNumber); 
	            if (r == null) {
	               // First cell in the row, create
	               r = sheet1.createRow(lineNumber);
	            }
	            c = r.getCell(0); 
	            if (c == null) {
	                // New cell
	                c = r.createCell(0);
	            }
	            c.setCellValue("Point #");   
	            c = r.getCell(1); 
	            if (c == null) {
	                // New cell
	                c = r.createCell(1);
	            }
	            c.setCellValue("Total tours length");            
	            c = r.getCell(2); 
	            if (c == null) {
	                // New cell
	                c = r.createCell(2);
	            }
	            c.setCellValue("Amplitude of tours");	            
	            c = r.getCell(3); 
	            if (c == null) {
	                // New cell
	                c = r.createCell(3);
	            }
	            c.setCellValue("List with cost of subtours");	
	            
	            lineNumber++;
	            for (int i = 0; i < bestSoFarPareto.size(); i++) {
	            	r = sheet1.getRow(i + lineNumber); 
		            if (r == null) {
		               // First cell in the row, create
		               r = sheet1.createRow(i + lineNumber);
		            }
		            //write point id
		            c = r.getCell(0); 
		            if (c == null) {
		                // New cell
		                c = r.createCell(0, Cell.CELL_TYPE_NUMERIC);
		            }
		            c.setCellValue(i + 1);
		            //write total cost and amplitude
	            	for (int indexObj = 0; indexObj < TSP_ACO.k; indexObj++) {
	            		c = r.getCell(indexObj + 1); 
			            if (c == null) {
			                // New cell
			                c = r.createCell(indexObj + 1, Cell.CELL_TYPE_NUMERIC);
			            }
			            c.setCellValue(bestSoFarPareto.get(i).costObjectives[indexObj]);
	            	}
	            	//write cost of each individual subtour
		            for (int j = 0; j < bestSoFarPareto.get(i).tour_lengths.length; j++) {
			            c = r.getCell(j + 3); 
			            if (c == null) {
			                // New cell
			                c = r.createCell(j + 3);
			            }
			            c.setCellValue(bestSoFarPareto.get(i).tour_lengths[j]);
		            }
	            }
    
				//Write the workbook in file system
			    FileOutputStream out = new FileOutputStream(new File(filePath1));
			    workbook1.write(out);
			    out.close();
			    
			    int nrOfRun = trial + 1;
			    System.out.println("\nRun #" + nrOfRun + " written Pareto front points successfully on disk.\n");
	    	}
			catch (Exception e) {
			    e.printStackTrace();
			}

	    }
	    else {
	    	System.out.println(" File " + filePath1 + " doesn't exists" );
	    }

    }
    
    //save in a .txt output file the best solution resulted after a run to be later used when
    //computing the Pareto front
    static void writeParetoSolutions(ArrayList<Ant> bestSoFarPareto) {
    	File f = new File(filePath2);
    	double[] objValues = new double[TSP_ACO.k];
    	
    	try {
    		BufferedWriter bw = new BufferedWriter(new FileWriter(f, true));
    			
    		for (int i = 0; i < bestSoFarPareto.size(); i++) {
    			if (rowIndex > 0) {
    				bw.newLine();
    			}
    			bw.write(rowIndex + "\t");
    			//get values total cost and amplitude
    			for (int indexObj = 0; indexObj < TSP_ACO.k; indexObj++) {
    				objValues[indexObj] = bestSoFarPareto.get(i).costObjectives[indexObj];
    			}
    			bw.write(objValues[0] + "\t");
    			//write cost of each individual subtour
    			for (int j = 0; j < bestSoFarPareto.get(i).tour_lengths.length; j++) {
    				bw.write(bestSoFarPareto.get(i).tour_lengths[j] + "\t");
    			}
    			bw.write(objValues[1] + "\t");
    			
    			rowIndex++;
    		}
    		bw.newLine();
    		bw.close();
        }
    	catch (IOException e) {
    		System.out.println("error writing file");
        }
    }

	public static String[] getTours() {
		return tours;
	}

	public static void setTours(String[] tours) {
		Utilities.tours = tours;
	}

	public static int[] getNrCities() {
		return nrCities;
	}

	public static void setNrCities(int[] nrCities) {
		Utilities.nrCities = nrCities;
	}

	public static double[] getSubtoursCost() {
		return subtoursCost;
	}

	public static void setSubtoursCost(double[] subtoursCost) {
		Utilities.subtoursCost = subtoursCost;
	}

	public static double getTotalCost() {
		return totalCost;
	}

	public static void setTotalCost(double totalCost) {
		Utilities.totalCost = totalCost;
	}

	public static ArrayList<Double> getIterTotalCost() {
		return iterTotalCost;
	}

	public static void setIterTotalCost(ArrayList<Double> iterTotalCost) {
		Utilities.iterTotalCost = iterTotalCost;
	}
    
    
    
}
