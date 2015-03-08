void boundary(double *y, int mesh_tracker)
{
	int i;

	if(mesh_tracker==0)
	{
		//========================================================================================//
		//the finest grid=========================================================================//
		//I am going to impose a NEUMANN boundary condition
	
		y[start[0]]=y[start[0]+2]-(LEFT*(2*dx[0]));//I am trying to impose a central difference at y[1]; the physical left boundary of the system in the finest grid 
	
		y[end[0]]=y[end[0]-2]+(RIGHT*(2*dx[0]));//I am trying to impose a central difference at y[nodes_x-2]; the physical right boundary of the system in the finest grid
	
		//========================================================================================//
	}

	else
	{
		//========================================================================================//
		//for coarser grids======================================================================// 
		//the boundary conditions for the other meshes are dirichlet's
	
	
		y[start[mesh_tracker]]=0.0;

		y[end[mesh_tracker]]=0.0;
	}
	
}
