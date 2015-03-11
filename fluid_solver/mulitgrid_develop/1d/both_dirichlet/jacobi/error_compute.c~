double comp_error(double *y_old,double *y)
{
	double error=0.0;
	int i;
	for(i=1;i<nodes_x-1;i++)//leaving out the boundary points for the Dirichlet boundary condition 	
	{
		error+=(y[i]-y_old[i])*(y[i]-y_old[i]);
	}

	error=sqrt(error*dx);
	
	return error;
}
