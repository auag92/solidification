void rb_solve(double *y,double *f,double *a,int start,double *error)
{
	int i;
	double y_star,new_val;

	double a_for_back;
	

	for(i=start;i<=nodes_x-2;i+=2)//I am not solving at the left end; //see the 'start' values passed    
	{
		//creating 'a_for_back'
		a_for_back=a[i-1]+a[i];

		y_star=(y[i-1]*a[i-1]+y[i+1]*a[i]-dx*dx*f[i])/a_for_back;
		

		new_val=y[i]+OMEGA*(y_star-y[i]);	

		*error+=(new_val-y[i])*(new_val-y[i]);
		

		y[i]=new_val;
	}

	
}				
