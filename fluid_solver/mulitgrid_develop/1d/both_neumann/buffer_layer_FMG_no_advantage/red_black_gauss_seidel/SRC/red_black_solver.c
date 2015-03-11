void rb_solve(double *y,double *f,int start,double *error)
{
	int i;
	double y_star,new_val;
	

	for(i=start;i<nodes_x-1;i+=2)//leaving out the left and right boundary points corresponding to array indices '0' and 'nodes_x-1' as they will be fixed by the boundary condition   
	{
		y_star=0.5*(y[i-1]+y[i+1]-dx*dx*f[i]);
		new_val=y[i]+OMEGA*(y_star-y[i]);	

		*error+=(new_val-y[i])*(new_val-y[i]);

		y[i]=new_val;
	}

	
}				
