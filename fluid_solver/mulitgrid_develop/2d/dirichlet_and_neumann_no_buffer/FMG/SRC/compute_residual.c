double comp_res(double *y,double *f,int mesh_tracker,long temp_x,long temp_y,long i,long j)
{
	
	
	long array_index,array_index_up,array_index_below,array_index_left,array_index_right;

	double residual;

	//int tag;

	
	//this are all calculations in the fine meshes
	array_index=(start[mesh_tracker]+temp_x)+nodes[mesh_tracker].x*temp_y;

	array_index_up=(start[mesh_tracker]+temp_x)+nodes[mesh_tracker].x*(temp_y+1);
	
	array_index_below=(start[mesh_tracker]+temp_x)+nodes[mesh_tracker].x*(temp_y-1);

	array_index_left=(start[mesh_tracker]+(temp_x-1))+nodes[mesh_tracker].x*temp_y;
	
	array_index_right=(start[mesh_tracker]+(temp_x+1))+nodes[mesh_tracker].x*temp_y;


	if(temp_x==0 && temp_y==0)
	{
		//tag=1;	
		residual=f[array_index]-(((y[array_index_right]-y[array_index])/sq_dx[mesh_tracker]) + ((y[array_index_up]-y[array_index])/sq_dy[mesh_tracker]));
	}

	else if(temp_x==0 && temp_y!=0)
	{
		//tag=2;
		residual=f[array_index]-(((-2.0*y[array_index]+2.0*y[array_index_right])/sq_dx[mesh_tracker]) + ((y[array_index_up]-2.0*y[array_index]+y[array_index_below])/sq_dy[mesh_tracker]));	
	}

	else if(temp_x!=0 && temp_y==0)
	{
		//tag=3;
		residual=f[array_index]-(((y[array_index_left]-2.0*y[array_index]+y[array_index_right])/sq_dx[mesh_tracker]) + ((2.0*y[array_index_up]-2.0*y[array_index])/sq_dy[mesh_tracker]));
	}
		
	else
	{
		//tag=4;			 
		residual=f[array_index]-(((y[array_index_left]-2.0*y[array_index]+y[array_index_right])/sq_dx[mesh_tracker]) + ((y[array_index_up]-2.0*y[array_index]+y[array_index_below])/sq_dy[mesh_tracker]));
	}


	//printf("res=%lf\n",residual);

	//getchar();

	//if(fabs(residual)>1e-4)
//	{
//		printf("res=%lf at i=%ld and j=%ld tag=%d\n",residual,i,j,tag);
//		printf("temp_x=%ld temp_y=%ld\n",temp_x,temp_y);
//		getchar();
//	}	

	return residual;
}
