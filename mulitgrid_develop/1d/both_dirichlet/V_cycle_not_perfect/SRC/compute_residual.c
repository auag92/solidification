void comp_res(double *y,double *f,int count)
{
	double res_2i,res_2i_min_1,res_2i_plus_1;
	int array_index_2i,array_index_2i_min_1,array_index_2i_plus_1,array_index_2i_min_2,array_index_2i_plus_2;
	int i;	

	int array_index_i;
	
	for(i=1;i<nodes_x[count+1]-1;i++)//the range of iterations is selected by the algorithm 
	//I am iterating over the points in the coarser grid	
	{
		//===========================================================================//
		//computing the residuals in the finer grid
		array_index_2i=left[count]+2*i;
		array_index_2i_min_1=left[count]+2*i-1;
		array_index_2i_plus_1=left[count]+2*i+1;
		
		array_index_2i_min_2=left[count]+2*i-2;
		array_index_2i_plus_2=left[count]+2*i+2;	

		res_2i=f[array_index_2i]-((y[array_index_2i_min_1]-2.0*y[array_index_2i]+y[array_index_2i_plus_1])/(dx[count]*dx[count]));

		res_2i_min_1=f[array_index_2i_min_1]-((y[array_index_2i_min_2]-2.0*y[array_index_2i_min_1]+y[array_index_2i])/(dx[count]*dx[count]));	

		res_2i_plus_1=f[array_index_2i_plus_1]-((y[array_index_2i]-2.0*y[array_index_2i_plus_1]+y[array_index_2i_plus_2])/(dx[count]*dx[count]));
		//===========================================================================//	


		//===========================================================================//		
		//moving the errors to the coarser grid using the full weighing operator  	
		array_index_i=i+left[count+1];
		f[array_index_i]=0.25*(res_2i_min_1+2.0*res_2i+res_2i_plus_1);
		//===========================================================================//
	
	}
}				
		
