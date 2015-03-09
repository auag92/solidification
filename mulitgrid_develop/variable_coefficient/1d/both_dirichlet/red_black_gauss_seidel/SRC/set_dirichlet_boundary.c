void dirichlet_boundary(double *y)
{
		
	//at the left end (the Dirichlet BC)	
	y[0]=LEFT;	 	
	y[nodes_x-1]=RIGHT;

}
