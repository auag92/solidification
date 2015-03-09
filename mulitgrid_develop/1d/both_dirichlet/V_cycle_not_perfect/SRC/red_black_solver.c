void rb_solve(double *y,double *f,int start,double *error,int count)
{
	int i;
	double y_star,new_val;
	
	for(i=start+left[count];i<right[count];i+=2)//leaving out the boundary points for the Dirichlet boundary condition 
	{
		y_star=0.5*(y[i-1]+y[i+1]-dx[count]*dx[count]*f[i]);
		new_val=y[i]+omega*(y_star-y[i]);	

		*error+=(new_val-y[i])*(new_val-y[i]);

		y[i]=new_val;
	}

	
}				
