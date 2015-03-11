double comp_res(double *y,double *f,int mesh_tracker,long temp_x,long temp_y)
{
	long i,j;
	
	long array_index,array_index_up,array_index_below,array_index_left,array_index_right;

	double residual;

	
	//this are all calculations in the fine meshes
	array_index=(start[mesh_tracker]+temp_x)+nodes[mesh_tracker].x*temp_y;

	array_index_up=(start[mesh_tracker]+temp_x)+nodes[mesh_tracker].x*(temp_y+1);
	
	array_index_below=(start[mesh_tracker]+temp_x)+nodes[mesh_tracker].x*(temp_y-1);

	array_index_left=(start[mesh_tracker]+(temp_x-1))+nodes[mesh_tracker].x*temp_y;
	
	array_index_right=(start[mesh_tracker]+(temp_x+1))+nodes[mesh_tracker].x*temp_y;
		
		 
	residual=f[array_index]-(((y[array_index_left]-2.0*y[array_index]+y[array_index_right])/sq_dx[mesh_tracker]) + ((y[array_index_up]-2.0*y[array_index]+y[array_index_below])/sq_dy[mesh_tracker]));


	//printf("res=%lf\n",residual);

	return residual;
}
