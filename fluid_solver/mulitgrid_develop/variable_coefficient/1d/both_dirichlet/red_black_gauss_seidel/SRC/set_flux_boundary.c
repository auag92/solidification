void flux_boundary(double *f)
{
	
	//at the right end (the Neumann BC) 
	f[nodes_x-1]-=(RIGHT/dx);	
}
