import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;
public class CallableAdder implements Callable<Circles> {
       public static int numOfTimes = 0 ; 
       double [][] finalGrid ; 
       ArrayList<Coordinates> centres;
       List<Coordinates> address;
       int num_i ; 
       int end_point ; 
       double max_llh ; 
       int flag ; 
       double rmin ; 
       double rmax ; 
	int study_area;
	double llh_limit ; 
       CallableAdder(int i  ,double[][] finalGrid , ArrayList<Coordinates> centres, List<Coordinates> address,int stud,double rmin , double rmax, double llh_limit)
       {	 this.num_i = i ; 
             this.finalGrid=finalGrid;
             this.centres = centres;             
             this.address = address ; 
             this.end_point = 0 ;
             this.max_llh = Double.MIN_VALUE;
      	    this.study_area = stud;
      	    this.flag = 0 ; 
      	    this.rmin = rmin ; 
      	    this.rmax = rmax ; 
	 }          
       public Circles call() throws Exception {
    	   for (int j = 0; j < centres.size() ; j++) {	// For All Other centers in the list 
				if(num_i != j && finalGrid[num_i][j] > rmin && finalGrid[num_i][j] < rmax ) {	// If they are different centers
					 int numOfPoints = 0 ; // Finding number of points in each circle......
					 flag = 1; 
					 for (int k = 0; k < address.size() ; k++) { // For each point in the dataset find if the point lies in the circle or not ? 
						if(isInsideaCircle(centres.get(num_i),finalGrid[num_i][j] , address.get(k))) {
							numOfPoints++;
						} // end inner if 
					} // End the K thing 
					//Now Calculate log likelihood of each circle. 
					// System.out.println(numOfPoints);
					  double current_llh = loglikelihood(numOfPoints, address.size(), 3.14*finalGrid[num_i][j]*finalGrid[num_i][j], study_area);	    			   
	    			   if(!Double.isInfinite(current_llh)&&current_llh > max_llh){
	    				   max_llh = current_llh;
	    				   end_point = j;
	    			   }
	    			   
				} // end if
			}// end j for
    	   numOfTimes++;
    	   if(max_llh > llh_limit ) {
    	     return new Circles(max_llh, num_i, end_point, finalGrid[num_i][end_point],true);
    	   }
    	   return new Circles(max_llh, num_i, end_point, finalGrid[num_i][end_point],false);
       }
       private static boolean isInsideaCircle(Coordinates coordinates, double d, Coordinates coordinates2) {
   		if(distance(coordinates.latitude, coordinates.longitude, coordinates2.latitude, coordinates2.longitude) <= d ) {// If distance between the centers <= 
   			// rad then inside the circle 
   			return true ; 
   		}
   		else {
   			return false;
   		}
   		
   	}
    private static double loglikelihood(int point_count, int P, double circle_area, int study_area){
   		
   		double B = (P*circle_area)/study_area;// Expected Number of points the circle.
   		
   		if(point_count > B){
   			return Math.log(Math.pow(point_count/B, point_count) * Math.pow((P - point_count)/(P - B), P - point_count));
   		}
   		else{
   			return Math.log(0);
   		}
   	}
    private static double distance(double lat1, double lon1, double lat2, double lon2) {
		  double p = 0.017453292519943295;    // Math.PI / 180
		  //char c = Math.cos;
		  double a = 0.5 - Math.cos((lat2 - lat1) * p)/2 + 
				  Math.cos(lat1 * p) * Math.cos(lat2 * p) * 
		          (1 - Math.cos((lon2 - lon1) * p))/2;

		  return 12742 * Math.asin(Math.sqrt(a)); //In KM
		}
}