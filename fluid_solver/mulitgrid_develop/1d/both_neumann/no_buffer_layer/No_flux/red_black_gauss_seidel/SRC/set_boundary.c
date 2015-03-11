void boundary(double *y)
{
	//I am going to impose a NEUMANN boundary condition
	y[0]=y[2]-(LEFT*(2*dx));//I am trying to impose a central difference at y[1]; the physical left boundary of the system 
	y[nodes_x-1]=y[nodes_x-3]+(RIGHT*2*dx);//I am trying to impose a central difference at y[nodes_x-2]; the physical right boundary of the system
}
