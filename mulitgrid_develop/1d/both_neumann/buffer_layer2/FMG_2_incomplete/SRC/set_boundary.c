void boundary(double *y, int mesh_tracker)
{
	int i;
	//========================================================================================//
	//the finest grid=========================================================================//
	//I am going to impose a NEUMANN boundary condition
	
	y[start[0]]=y[start[0]+1]-(LEFT*(dx[0]));//I am trying to impose a forward difference at y[1]; the physical left boundary of the system in the finest grid 
	
	y[end[0]]=y[end[0]-1]+(RIGHT*(dx[0]));//I am trying to impose a forward difference at y[nodes_x-2]; the physical right boundary of the system in the finest grid
	
	//========================================================================================//
	

	
	
}
