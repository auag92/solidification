void rb_solve(double *y,double *f,int start,double *error)
{
	int i;
	double y_star,new_val;
	

	for(i=start;i<nodes_x;i+=2)//I am not solving at the left end; //see the 'start' values passed    
	{
		if(i==nodes_x-1)
			y_star=y[i-1]-dx*dx*f[i];	
			
		else
			y_star=0.5*(y[i-1]+y[i+1]-dx*dx*f[i]);
		

		new_val=y[i]+OMEGA*(y_star-y[i]);	

		*error+=(new_val-y[i])*(new_val-y[i]);
		

		y[i]=new_val;
	}

	
}				
