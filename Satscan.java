import java.io.*;
import java.util.*;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import javax.swing.plaf.basic.BasicInternalFrameTitlePane.MaximizeAction;

import java.math.*;





public class Satscan {
	static int mthreads = 1 ;
	static int cthreads =1 ; 
	//Distance function calculates the straight distance, same as we draw straight line between two points 
	//This is totally different from google map distance, since it uses road dataset for calculating distance.
	public static double distance(double lat1, double lon1, double lat2, double lon2) {
		  double p = 0.017453292519943295;    // Math.PI / 180
		  //char c = Math.cos;
		  double a = 0.5 - Math.cos((lat2 - lat1) * p)/2 + 
				  Math.cos(lat1 * p) * Math.cos(lat2 * p) * 
		          (1 - Math.cos((lon2 - lon1) * p))/2;

		  return 12742 * Math.asin(Math.sqrt(a)); //In KM
		}
	
	public static double loglikelihood(int point_count, int P, double circle_area, int study_area){
		
		double B = (P*circle_area)/study_area;// Expected Number of points the circle.
		
		if(point_count > B){
			return Math.log10(Math.pow(point_count/B, point_count) * Math.pow((P - point_count)/(P - B), P - point_count));//math.log is base e
		}
		else{
			//System.out.println("log0");
			return Math.log10(0);
		}
	}
	
	
	public static List<Coordinates> random_set(double min_lat, double max_lat, double min_long, double max_long, int number){
		
		List<Coordinates> Address = new ArrayList<>();
		
		for(int i=0;i<number;i++){
			double latitude = min_lat + (Math.random()*(max_lat - min_lat));
			double longitude = min_long + (Math.random()*(max_long - min_long));
			
			Address.add(new Coordinates(latitude, longitude));
		}
		
		return Address;
		
	}
	
	public static List<Double> MCSimulations(int simulations,final int P,final  double min_lat,final double max_lat,final double min_long,final  double max_long,final int study_area) throws InterruptedException, ExecutionException{
		//Integer threads = Runtime.getRuntime().availableProcessors();
	    ExecutorService service = Executors.newFixedThreadPool(mthreads);
	   // System.out.println(threads);
	    List<Future<Double>> futures = new ArrayList<Future<Double>>();
	    
		
		for(int a=0;a<simulations;a++){
			        Callable<Double> callable = new Callable<Double>() {
			            public Double call() throws Exception {
			            	 
			                Integer output = 1;
			                output *= 100 ; 
			                List<Coordinates> Address = random_set(min_lat, max_lat, min_long, max_long, P);
			    			
			    			double[][] radius_grid = new double[Address.size()][Address.size()];
			    		       
			    		       for(int i=0;i<Address.size();i++){
			    		    	   
			    		    	   for(int j=0;j<i;j++){
			    		    		  radius_grid[i][j] = distance(Address.get(i).latitude, Address.get(i).longitude, Address.get(j).latitude, Address.get(j).longitude);
			    		    		  radius_grid[j][i] = radius_grid[i][j];
			    		    	   
			    		    	   }
			    		    	   radius_grid[i][i] = 0.0;
			    		       }
			    			
			    		       double max_llh = Double.MIN_VALUE;
			    		       
			    		       
			    		       for(int i=0;i<Address.size();i++){//each point as centre of candidate circle
			    		    	   
			    		    	   for(int j=0;j<Address.size();j++){//All other points as radius 
			    		    		   if(i!=j){
			    		    			   
			    		    			   int point_count = 0;//counting number of points inside candidate circle
			    		    			   
			    		    			   for(int k=0;k<Address.size();k++){
			    		    				   
			    		    				   if(radius_grid[i][k] <= radius_grid[i][j]){
			    		    					   point_count++;
			    		    				   }
			    		    			   }
			    		    			   
			    		    			   //calculate log likelihood of each candidate circle
			    		    			   double current_llh = loglikelihood(point_count, Address.size(), 3.14*radius_grid[i][j]*radius_grid[i][j], study_area);
			    		    			   
			    		    			   if(!Double.isFinite(current_llh)&& current_llh > max_llh){
			    		    				   max_llh = current_llh;
			    		    			   }
			    		    		   }
			    		    	   }
			    		    	   
			    		       }//calculating max llh among all candidate circles
			    		       
			    		       return max_llh;
			    	
			               
			            }
			        };
			        futures.add(service.submit(callable));
				}//m silulations

	    service.shutdown();
		  List<Double> outputs = new ArrayList<Double>();
		    for (Future<Double> future : futures) {
		        outputs.add(future.get());
		    }
		    Collections.sort(outputs);
			Collections.reverse(outputs);
		    return outputs;
	}
	
	/**
	 * @param args
	 * @throws IOException
	 * @throws InterruptedException
	 * @throws ExecutionException
	 */
	public static void main(String[] args) throws IOException, InterruptedException, ExecutionException{
		
		/*
		 * Code for getting the input from the file
		 * */
         	 double lat_spacing=0.3;
		double long_spacing = 0.005 ; // This is the variable used while enumerating all possible points in that rectangle 
		long startTime = System.currentTimeMillis();
		BufferedReader br = new BufferedReader(new FileReader("break.csv"));
        String line = br.readLine();
        String splitBy = ",";
        double min_lat = 90, max_lat = -90;
        double min_long = 180, max_long = -180;
        List<Coordinates> Address = new ArrayList<>();
        double rmin = 0.1 ; 
        double rmax = 4 ; 
        
        while((line = br.readLine()) != null){
            String[] b = line.split(splitBy);
            //System.out.println(b[0]+": "+b[1]+" "+b[2]);
            double latitude = Double.parseDouble(b[13]);
            double longitude = Double.parseDouble(b[14]);
          //  System.out.println("lat = " + latitude + "long = " + longitude);
            Address.add(new Coordinates(latitude, longitude));
            //System.out.println(distance(-111.9606778, 33.38364645, lat2, long2));
            
            if(latitude < min_lat){
            	min_lat = latitude;
            }
            else if(latitude > max_lat){
            	max_lat = latitude;
            }
            
            if(longitude < min_long){
            	min_long =longitude; 
            }
            else if(longitude > max_long){
            	 max_long = longitude;
            }
            
       }
       br.close();
      // System.out.println("min = "+min_long + "max = " + max_long);
	int study_area = (int)(distance(min_lat , min_long , min_lat , max_long)* distance(min_lat,min_long , max_lat , min_long ));
       int numberofsimulations = 800 ;
       List<Double> simulations_llh = MCSimulations(numberofsimulations,Address.size(), min_lat, max_lat, min_long, max_long,study_area);
       double alpha = 0.1;
       double llh_limit = simulations_llh.get((int) Math.round(alpha*(numberofsimulations+1)));
       	
       // -------------------- Total Latitudes Possible ------------------- Storing it into an arrayList 
       ArrayList<Double> allPossibleLatitutes = new ArrayList<Double>() ; 
       // Enumerated all possible latitudes in the arraylist 
       enumerateAllpossible(allPossibleLatitutes, min_lat, max_lat, lat_spacing); 
       // -------------------- Total Longitutes Possible ------------------- Storing it into an arrayList 
       ArrayList<Double> allPossibleLongitutes = new ArrayList<Double>();

	   System.out.println("Size of all pos lats = " +allPossibleLatitutes.size());
       enumerateAllpossible(allPossibleLongitutes, min_long, max_long , long_spacing);

       System.out.println("Size of all pos lats = "+allPossibleLongitutes.size());
       //1) ENUMERATE ALL POSSIBLE CENTERS OF THE GRID 
       ArrayList <Coordinates> centres = new ArrayList<Coordinates>();
       
       //2) Add all possible centers of the grid into the created array-list using the below function 
       addPossibleCentres(centres,allPossibleLatitutes,allPossibleLongitutes);
       
       //3) Now find the radius of all possible circles from each other and store their radius in the 2Darray
       double [][] finalGrid = new double [centres.size()][centres.size()];
       
       //3.1) populate the grid with the radius values ... 
       populateFinalGrid(finalGrid, centres);
       
       //4.0) Create candidate circles List and now populate the candidate circles 
       List<Circles> Candidate_circles = new ArrayList<>();
       
       //4.1) Populate the list with the Candidate circles 
       populateCandidateCirclesCopy(Candidate_circles,finalGrid,centres,Address,study_area,llh_limit,rmin,rmax);
       
       System.out.println("Size Of Candidate Circles = " + Candidate_circles.size());
       
       //Upto here we can get candidate circle list and MCS.

      
      
      
      //Upto here we can get candidate circle list and MCS.
      int[] overlapping_state = new int[Candidate_circles.size()];
      Collections.sort( Candidate_circles);
      for(int i=0;i<Candidate_circles.size();i++){
   	   Circles c1 = Candidate_circles.get(i);
   	  
	    	   if(overlapping_state[i] == 0){
	    		   overlapping_state[i] = 1;
	    		   Circles Hotspot = c1;

	    		   for(int j=0;j<Candidate_circles.size();j++){
	    			   
	    			   if(overlapping_state[j] == 0){
	    				   Circles c2 = Candidate_circles.get(j);
	    				   
	    				   if(c1.radius + c2.radius > finalGrid[c1.centre][c2.centre]){
	    					   overlapping_state[j] = 1;
	    				   }
	    			   }
	    		   }
	    		   System.out.println("Centre: "+ centres.get(Hotspot.centre).latitude + ", "+ centres.get(Hotspot.centre).longitude+" Radius: "+Hotspot.radius
						+"llr :" + Hotspot.llh);
	    	   }
      }//all candidate circles
       long endTime   = System.currentTimeMillis();
       long totalTime = endTime - startTime;
       System.out.println("Run Time for Satscan: "+totalTime);
       
	}//main

	private static void populateCandidateCircles(List<Circles> candidate_circles, double[][] finalGrid, ArrayList<Coordinates> centres, List<Coordinates> address,int study_area,double llh_limit) {
	//	System.out.println("entered"); // Not entered 
		for (int i = 0; i < centres.size() ; i++) {	// For Each Center in the Grid 
			double max_llh = Double.MIN_VALUE;
			int endpoint = 0 ;
			for (int j = 0; j < centres.size() ; j++) {	// For All Other centers in the list 
				if(i != j) {	// If they are different centers
					 int numOfPoints = 0 ; // Finding number of points in each circle......
					 for (int k = 0; k < address.size() ; k++) { // For each point in the dataset find if the point lies in the circle or not ? 
						if(isInsideaCircle(centres.get(i),finalGrid[i][j] , address.get(k))) {
							numOfPoints++;
						} // end inner if 
					} // End the K thing 
					//Now Calculate log likelihood of each circle. 
					// System.out.println(numOfPoints);
					  double current_llh = loglikelihood(numOfPoints, address.size(), 3.14*finalGrid[i][j]*finalGrid[i][j], study_area);	    			   
	    			   if(!Double.isInfinite(current_llh)&&current_llh > max_llh){
	    				   max_llh = current_llh;
	    				   endpoint = j;
	    			   }
	    			   
				} // end if
			}// end j for
			
			candidate_circles.add(new Circles(max_llh,i, endpoint, finalGrid[i][endpoint],true));
			
		}// end i for 
	}
	
	private static void populateCandidateCirclesCopy(List<Circles> candidate_circles,double[][] finalGrid, ArrayList<Coordinates> centres, List<Coordinates> address,int study_area, double llh_limit,
			double rmin,double rmax) {
		//	System.out.println("entered"); // Not entered 
		 
	    ExecutorService service = Executors.newFixedThreadPool(cthreads);
	    List <Future<Circles>> list = new ArrayList<Future<Circles>>(); 
			for (int i = 0; i < centres.size() ; i++) {	// For Each Center in the Grid 
				Future <Circles> future = service.submit(new CallableAdder(i, finalGrid, centres, address,study_area,rmin,rmax,llh_limit));
				list.add(future);
			}// end i for
			service.shutdown();
			//System.out.println(CallableAdder.numOfTimes);
			// = new ArrayList<Circles>();
			for(Future<Circles> fut : list) {
				try {
					if(fut.get().verified==true)
						{candidate_circles.add(fut.get());}
				} catch (InterruptedException | ExecutionException e) {
					System.out.println("Exception arised due to the future object ");
					e.printStackTrace();
				}
			}
		
			
		}
	/**
	 * This Function is created for checking whether the point 
	 * @param coordinates : Center 
	 * @param d	: Radius of the circle 
	 * @param coordinates2	: Point to be checked 
	 * @return	: True if Point inside , else false !!
	 */
	private static boolean isInsideaCircle(Coordinates coordinates, double d, Coordinates coordinates2) {
		if(distance(coordinates.latitude, coordinates.longitude, coordinates2.latitude, coordinates2.longitude) <= d ) {// If distance between the centers <= 
			// rad then inside the circle 
			return true ; 
		}
		else {
			return false;
		}
		
	}

	private static void populateFinalGrid(double[][] finalGrid, ArrayList<Coordinates> centres) {
		for (int i = 0; i < centres.size() ; i++) {
			for (int j = 0; j < i ; j++) {
				// calculating the radius 
				finalGrid[i][j] = distance(centres.get(i).latitude, centres.get(i).longitude, centres.get(j).latitude, centres.get(j).longitude);
				//System.out.println(finalGrid[i][j]);
				finalGrid[j][i] = finalGrid[i][j];
			}
		 finalGrid[i][i] = 0 ;
		}
	}

	/**
	 * This Function is Just a utility function which creates all possible centers in the grid  
	 * @param centres
	 * @param allPossibleLatitutes
	 * @param allPossibleLongitutes
	 */
	private static void addPossibleCentres(ArrayList<Coordinates> centres, ArrayList<Double> allPossibleLatitutes,
			ArrayList<Double> allPossibleLongitutes) {
			for(int i = 0 ; i < allPossibleLatitutes.size() ; i++ ) {
				for (int j = 0; j < allPossibleLongitutes.size(); j++) {
					centres.add(new Coordinates(allPossibleLatitutes.get(i), allPossibleLongitutes.get(j)));
				}
			}	
	}

	
	/**
	 * Function fills the Arraylist with allPossibleLatitutes | allPossibleLongitutes  Based on the minimum and maximum parameter 
	 * @param allPossibleLatitutes : Array List to be Filled 
	 * @param min_lat 	: minimum possible latitude or longitude
	 * @param max_lat	: max possible latitude or longitude 
	 */
	private static void enumerateAllpossible(ArrayList<Double> allPossibleLatitutes, double min_lat, double max_lat,double spacing) {
		for(double i = min_lat ; i<=max_lat ; i+= spacing) {
			allPossibleLatitutes.add(i);
		}
	}
}