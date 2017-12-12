import java.io.*;
import java.util.*;
import java.math.*;



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
			return Math.log(Math.pow(point_count/B, point_count) * Math.pow((P - point_count)/(P - B), P - point_count));
		}
		else{
			return Math.log(0);
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
	
	public static List<Double> MCSimulations(int simulations,int P, double min_lat, double max_lat, double min_long, double max_long){
		
		List<Double> simulations_llh = new ArrayList<>();
		
		for(int a=0;a<simulations;a++){
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
		    			   double current_llh = loglikelihood(point_count, Address.size(), 3.14*radius_grid[i][j]*radius_grid[i][j], 2500);
		    			   
		    			   if(current_llh > max_llh){
		    				   max_llh = current_llh;
		    			   }
		    		   }
		    	   }
		    	   
		       }//calculating max llh among all candidate circles
		       
		       simulations_llh.add(max_llh);
		}//m silulations
		
		Collections.sort(simulations_llh);
		Collections.reverse(simulations_llh);
		
		return simulations_llh;
	}
	
	public static void main(String[] args) throws IOException{
		long startTime = System.currentTimeMillis();
		//System.out.println("Started");
		BufferedReader br = new BufferedReader(new FileReader("test.csv"));
        String line = br.readLine();
        String splitBy = ",";
        double min_lat = 90, max_lat = 0;
        double min_long = 180, max_long = 0;
        List<Coordinates> Address = new ArrayList<>();
        
        while((line = br.readLine()) != null){
            String[] b = line.split(splitBy);
            //System.out.println(b[0]+": "+b[1]+" "+b[2]);
            double latitude = Double.parseDouble(b[0]);
            double longitude = Double.parseDouble(b[1]);
            
            Address.add(new Coordinates(latitude, longitude));
            //System.out.println(distance(-111.9606778, 33.38364645, lat2, long2));
            
            if(latitude < min_lat){
            	min_lat = latitude;
            }
            else if(latitude > max_lat){
            	max_lat = latitude;
            }
            
            if(longitude < min_long){
            	longitude = min_long;
            }
            else if(longitude > max_long){
            	longitude = max_long;
            }
            
       }
       br.close();
        
       //System.out.println("Nirmal to Bhainsa(testing):"+distance(19.096412, 78.342975, 19.103066, 77.965299));
       //Actual distance is 41.5 km, but it is showing 39.69 km, because there are two curves.
       //System.out.println("Testing loglikelihood with paper values "+loglikelihood(199, 800, 186.62, 2500));
       
       int numberofsimulations = 500;
       List<Double> simulations_llh = MCSimulations(numberofsimulations,Address.size(), min_lat, max_lat, min_long, max_long);
       double alpha = 0.5;
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
       
       for(int i=0;i<Address.size();i++){//each point as centre of candidate circle
    	   
    	   double max_llh = Double.MIN_VALUE;
    	   int endpoint = 0;
    	   
    	   for(int j=0;j<Address.size();j++){//All other points as radius 
    		   if(i!=j){
    			   
    			   int point_count = 0;//counting number of points inside candidate circle
    			   
    			   for(int k=0;k<Address.size();k++){
    				   
    				   if(radius_grid[i][k] <= radius_grid[i][j]){
    					   point_count++;
    				   }
    			   }
    			   
    			   //calculate log likelihood of each candidate circle
    			   double current_llh = loglikelihood(point_count, Address.size(), 3.14*radius_grid[i][j]*radius_grid[i][j], 2500);
    			   
    			   if(current_llh > max_llh){
    				   max_llh = current_llh;
    				   endpoint = j;
    			   }
    		   }
    	   }
    	   
    	   //for each point as centre and all other points as radius, we have got candidate circle of max log likelihood among them
    	   //Since all circles are overlapping for a single point as centre for all, we have taken only one circle of max llh
    	   
    	   Candidate_circles.add(new Circles(max_llh,i, endpoint, radius_grid[i][endpoint]));
       }
       
       //Upto here we can get candidate circle list and MCS.
       int[] overlapping_state = new int[Candidate_circles.size()];
       
       for(int i=0;i<Candidate_circles.size();i++){
    	   Circles c1 = Candidate_circles.get(i);
    	   if(c1.llh > llh_limit){
	    	   if(overlapping_state[i] == 0){
	    		   overlapping_state[i] = 1;
	    		   Circles Hotspot = c1;
	    		   double max_llh_overlapping = c1.llh;
	    		   
	    		   for(int j=0;j<Candidate_circles.size();j++){
	    			   
	    			   if(overlapping_state[j] == 0){
	    				   Circles c2 = Candidate_circles.get(j);
	    				   
	    				   if(c1.radius + c2.radius > radius_grid[c1.centre][c2.centre]){
	    					   overlapping_state[j] = 1;
	    					   if(c2.llh > max_llh_overlapping){
	    						   max_llh_overlapping = c2.llh;
	    						   Hotspot = c2;
	    					   }
	    				   }
	    			   }
	    		   }
	    		   
	    		   //Print hotspot
	    		   System.out.println("Centre: "+Address.get(Hotspot.centre).latitude + ", "+ Address.get(Hotspot.centre).longitude+" Radius: "+Hotspot.radius);
	    	   }
    	   }
    	   else{
    		   overlapping_state[i] = 1;
    	   }
       }//all candidate circles
       
       long endTime   = System.currentTimeMillis();
       long totalTime = endTime - startTime;
       System.out.println("Run Time for Satscan: "+totalTime);
       
	}//main
}