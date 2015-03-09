void boundary(double *f)
{
	//the specified fluxes appear at the rhs nodes corresponding to the boundary points
	
	//at the left end	
	f[0]+=(LEFT/dx);

	//at the right end
	f[nodes_x-1]-=(RIGHT/dx);		 	


}
