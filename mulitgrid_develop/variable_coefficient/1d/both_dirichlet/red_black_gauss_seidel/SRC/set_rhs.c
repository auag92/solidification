void rhs(double *f)
{	
	int i;

	for(i=1;i<nodes_x;i++)//I am following the convention of "Briggs" where the equations at the boundary nodes at which Neumann bundary conditons are imposed are scaled by 2   
	{
		#ifdef ZERO_RHS
		f[i]=0.0;
		#endif

		#ifdef X_RHS
		f[i]=SYS_LEFT_END+i*dx;
		#endif
	}

	
}						
