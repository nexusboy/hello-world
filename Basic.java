import java.io.*;
import java.util.*;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.math.*;

class Coordinates {
    final double latitude;
    final double longitude;

    Coordinates(double lat, double lon) {
        latitude = lat;
        longitude = lon;
    }
}

class Circles implements Comparable<Circles>{
	double llh;
	int centre;
	int endpoint;
	double radius;
	
	Circles(double llh, int centre,int endpoint, double radius) {
		this.llh = llh;
		this.centre = centre;
		this.endpoint = endpoint;
		this.radius = radius;
	}
	
	public int compareTo(Circles compareCircle){
		
		double curr_llh = ((Circles) compareCircle).llh;
		
		int j =  (curr_llh - this.llh) > 0 ? 1 : -1 ;
		return j ; 
	}
}

public class Basic {
	
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
	
	public static List<Double> MCSimulations(int simulations,final int P, final double min_lat, final double max_lat, final double min_long, final double max_long, final int study_area) throws InterruptedException, ExecutionException{
		Integer threads = Runtime.getRuntime().availableProcessors();
	    ExecutorService service = Executors.newFixedThreadPool(threads);
	   // System.out.println(threads);
	    List<Future<Double>> futures = new ArrayList<Future<Double>>();
	    
		List<Double> simulations_llh = new ArrayList<>();
		
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
			    		    			   
			    		    			   if(!Double.isInfinite(current_llh) && current_llh > max_llh){
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
		long startTime = System.currentTimeMillis();
		BufferedReader br = new BufferedReader(new FileReader("beep.csv"));
        String line = br.readLine();
        String splitBy = ",";
        double min_lat = 90, max_lat = -90;
        double min_long = 180, max_long = -180;
        List<Coordinates> Address = new ArrayList<>();
        
        while((line = br.readLine()) != null){
            String[] b = line.split(splitBy);
            //System.out.println(b[0]+": "+b[1]+" "+b[2]);
            double latitude = Double.parseDouble(b[1]);//13 for nanaimo, 1 for tempe
            double longitude = Double.parseDouble(b[0]);//14 for nanaimo, 0 for tempe
            
            Address.add(new Coordinates(latitude, longitude));
            //System.out.println(distance(-111.9606778, 33.38364645, lat2, long2));
            
            if(latitude < min_lat){
            	min_lat = latitude;
            }
            else if(latitude > max_lat){
            	max_lat = latitude;
            }
            
            if(longitude < min_long){
            	 min_long = longitude;
            }
            else if(longitude > max_long){
            	 max_long = longitude;
            }
            
       }
       br.close();
       //System.out.println("Border points: "+min_lat+" "+min_long+" "+max_lat+" "+max_long);
       int study_area = (int) (distance(min_lat, min_long, min_lat, max_long)*distance(min_lat, min_long, max_lat, min_long));
       //study_area = 2500;
      // System.out.println("Study area: "+study_area);
       
       //System.out.println("Distance testing: "+distance(33.4049868949, -111.9281948145, 33.4229380496, -111.9508124747));
       
       int numberofsimulations = 800;
       List<Double> simulations_llh = MCSimulations(numberofsimulations,Address.size(), min_lat, max_lat, min_long, max_long,study_area);
//       for(int i=0;i<simulations_llh.size();i++){
//    	   System.out.println(simulations_llh.get(i));
//       }
       double alpha = 0.001;
       double llh_limit = simulations_llh.get((int) Math.round(alpha*(numberofsimulations+1)));
       
       double[][] radius_grid = new double[Address.size()][Address.size()];
       
       for(int i=0;i<Address.size();i++){
    	   
    	   for(int j=0;j<i;j++){
    		  radius_grid[i][j] = distance(Address.get(i).latitude, Address.get(i).longitude, Address.get(j).latitude, Address.get(j).longitude);
    		  radius_grid[j][i] = radius_grid[i][j];
    	   }
    	   
    	   radius_grid[i][i] = 0.0;
       }
       
       List<Circles> Candidate_circles = new ArrayList<>();
       double rmin = 0.1;
       int point_count = 0;
       System.out.println("Number of points: "+Address.size()+" Study area: "+study_area);
       //System.out.println("LLH Testing(infinity) "+loglikelihood(632, 1720, 42.20795626, 137));
       
       
      
       
       
       for(int i=0;i<Address.size();i++){//each point as centre of candidate circle
    	   
    	   double max_llh = Double.MIN_VALUE;
    	   int endpoint = 0;
    	   int flag = 0;
    	   
    	   for(int j=0;j<Address.size();j++){//All other points as radius 
    		   if(i!=j && radius_grid[i][j] > rmin && radius_grid[i][j] < 4){
    			   
    			   point_count = 0;//counting number of points inside candidate circle
    			   flag = 1;
    			   
    			   for(int k=0;k<Address.size();k++){
    				   
    				   if(radius_grid[i][k] <= radius_grid[i][j]){
    					   point_count++;
    				   }
    			   }
    			   
    			   //calculate log likelihood of each candidate circle
    			   double current_llh = loglikelihood(point_count, Address.size(), 3.14*(radius_grid[i][j])*(radius_grid[i][j]), study_area);
    			   //System.out.println(current_llh);
    			   
    			   if(!Double.isInfinite(current_llh) && current_llh > max_llh){
    				   max_llh = current_llh;
    				   endpoint = j;
    			   }
    		   }
    	   }//All other points as radius 
    	   
    	   //for each point as centre and all other points as radius, we have got candidate circle of max log likelihood among them
    	   //Since all circles are overlapping for a single point as centre for all, we have taken only one circle of max llh
    	   
    	   if(flag == 1 && max_llh > llh_limit){
    		   Candidate_circles.add(new Circles(max_llh,i, endpoint, radius_grid[i][endpoint]));
    		   //System.out.println(max_llh+" "+ radius_grid[i][endpoint]+" "+(Address.size()/study_area)*3.14*(radius_grid[i][endpoint])*(radius_grid[i][endpoint])+" "+point_count);
    		
    	   }
    	   
       }
       
    
       
       //Upto here we can get candidate circle list and MCS.
       int[] overlapping_state = new int[Candidate_circles.size()];
       Collections.sort(Candidate_circles);
     
       
       for(int i=0;i<Candidate_circles.size();i++){
    	   Circles c1 = Candidate_circles.get(i);
	    	   if(overlapping_state[i] == 0){
	    		   overlapping_state[i] = 1;
	    		   Circles Hotspot = c1;
	    		   //double max_llh_overlapping = c1.llh;
	    		   int flag_print = 0;
	    		   flag_print = 1;
	    		   
	    		   for(int j=0;j<Candidate_circles.size();j++){
	    			   
	    			   if(overlapping_state[j] == 0){
	    				   Circles c2 = Candidate_circles.get(j);
	    				   
	    				   if(c1.radius + c2.radius > radius_grid[c1.centre][c2.centre]){
	    					   overlapping_state[j] = 1;
	    				   }
	    			   }
	    		   }//iterating through all overlapping circles
	    		   
	    		
	    		   
	    		   
	    	   }//array == 0
       }//all candidate circles
       
       
       long endTime   = System.currentTimeMillis();
       long totalTime = endTime - startTime;
       System.out.println("Run Time for Satscan(milli seconds): "+totalTime);
       
	}//main
}
