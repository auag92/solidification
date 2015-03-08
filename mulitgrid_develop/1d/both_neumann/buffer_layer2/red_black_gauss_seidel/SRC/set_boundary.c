void boundary(double *y)
{
	//I am going to impose a NEUMANN boundary condition
	y[0]=y[1]-(LEFT*dx);//I am trying to impose a forward difference at y[1]; the physical left boundary of the system 
	y[nodes_x-1]=y[nodes_x-2]+(RIGHT*dx);//I am trying to impose a forward difference at y[nodes_x-2]; the physical right boundary of the system
}
