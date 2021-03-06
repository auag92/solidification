void cr_a(double *a)
{
	int i,j;

	for(j=0;j<NO_OF_LEVELS;j++)
	{

		printf("at level=%d\n",j);
	
		for(i=0;i<mesh_intervals[j];i++)//'a' is known a-priori
		//I am going to set it everywhere regardless of the BCs   
		{

			if(j==0) //for the finest grid
			{
				#ifdef A_ONE//mimics the simpler problem we have started with
 				a[i]=1.0;
				#endif

				#ifdef A_SIX_TIMES_X_SQ
				a[i]=6.0*(SYS_LEFT_END+((i+0.5)*dx[j]))*(SYS_LEFT_END+((i+0.5)*dx[j]));//as the '0' node coincides with the left end of the system
				#endif

				#ifdef A_X_BY_100
				a[i]=(SYS_LEFT_END+((i+0.5)*dx[j]))/100.0;//as the '0' node coincides with the left end of the system
				#endif
			}

			else
			{
				a[i+start_mesh_inter[j]]=0.5*(a[2*i+start_mesh_inter[j-1]]+a[2*i+1+start_mesh_inter[j-1]]);//this is the averaging scheme gven in the book
			}	

			//if(j==1)
			//	printf("a1=%lf\ta2=%lf\n",a[2*i+start_mesh_inter[j-1]],a[2*i+1+start_mesh_inter[j-1]]);				

			printf("i=%d\ta=%lf\n",i,a[i+start_mesh_inter[j]]);
		}
	}

	getchar();
}		
