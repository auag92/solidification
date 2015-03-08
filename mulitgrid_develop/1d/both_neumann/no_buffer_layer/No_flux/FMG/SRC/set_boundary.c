void boundary(double *y)
{
	//setting the boundary conditon for the finest grid only
		
	//the boundary conditions for the other meshes are zero by default 
	
	y[start[0]]=LEFT;
	y[end[0]]=RIGHT;
	
	
}
