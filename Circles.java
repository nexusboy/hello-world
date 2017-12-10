class Circles implements Comparable<Circles>{
	double llh;
	int centre;
	int endpoint;
	double radius;
	boolean verified ; 
	Circles(double llh, int centre,int endpoint, double radius, boolean b) {
		this.llh = llh;
		this.centre = centre;
		this.endpoint = endpoint;
		this.radius = radius;
		this.verified = b ; 
	}
	
	public int compareTo(Circles compareCircle){
		
		double curr_llh = ((Circles) compareCircle).llh;
		
		int j =  (curr_llh - this.llh) > 0 ? 1 : -1 ;
		return j ; 
	}
}