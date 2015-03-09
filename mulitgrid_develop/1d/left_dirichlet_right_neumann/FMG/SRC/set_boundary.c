void boundary(double *y,double *f)
{
	//the specified fluxes alter the rhs
	//the specified fluxes appear at the rhs nodes corresponding to the boundary points in the finest grid only
	
	//at the left end (the Dirichlet BC)	
	y[start[0]]=LEFT;

	//at the right end
	f[end[0]]-=(RIGHT/dx[0]);		
	
	
}
