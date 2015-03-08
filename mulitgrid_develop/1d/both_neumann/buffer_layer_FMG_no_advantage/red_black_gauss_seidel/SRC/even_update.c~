void jacobi(double *y_old,double *y,double *f)
{
	int i;
	double y_star;
	for(i=1;i<nodes_x-1;i++)//leaving out the boundary points for the Dirichlet boundary condition 
	{
		y_star=0.5*(y_old[i-1]+y_old[i+1]-dx*dx*f[i]);
		y[i]=y_old[i]+omega*(y_star-y_old[i]);	
	}

	//imposing the boundary condition on the recently computed array
	y[0]=LEFT;
	y[nodes_x-1]=RIGHT;
}				
