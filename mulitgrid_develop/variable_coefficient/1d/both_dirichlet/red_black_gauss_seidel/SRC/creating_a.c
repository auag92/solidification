void cr_a(double *a)
{
	int i;
	
	for(i=0;i<NO_OF_INTERVALS;i++)//'a' is known a-priori
	//I am going to set it everywhere regardless of the BCs   
	{
		#ifdef A_ONE//mimics the simpler problem we have started with
 		a[i]=1.0;
		#endif

		#ifdef A_SIX_TIMES_X_SQ
		a[i]=6.0*(SYS_LEFT_END+((i+0.5)*dx))*(SYS_LEFT_END+((i+0.5)*dx));//as the '0' node coincides with the left end of the system
		#endif	
	}
}		
