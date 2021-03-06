void ana(double *ana_y,double *g,double *y)
{
	int p;

	int i,j;

	long array_index;

	double lambda_p;

	//==============================================================================//
	char fname[300];
	FILE *fp;
	//==============================================================================//


	sprintf(fname,"../DATA/ana_prof.dat");
	if((fp=fopen(fname,"w"))==NULL)
	{
		printf("analytical profile can't be written\n");
		exit(1);
	}



	for(j=0;j<=nodes.y-1;j++)
	{
		for(i=0;i<=nodes.x-1;i++)	
		{
			array_index=i+nodes.x*j;
			
			ana_y[array_index]=(1.0-(i*d_sp.x*i*d_sp.x));		

			for(p=0;p<NOFTER;p++)
			{
				lambda_p=((2*p+1)*M_PI)/2.0;

				ana_y[array_index]-=(4.0*pow(-1.0,p)*cos(lambda_p*i*d_sp.x)*cosh(lambda_p*j*d_sp.y))/(lambda_p*lambda_p*lambda_p*cosh(lambda_p));
			}

			ana_y[array_index]*=(-g[array_index]/2.0);	

			fprintf(fp,"%lf\t%lf\t%lf\t%lf\n",i*d_sp.x,j*d_sp.y,ana_y[array_index],fabs(y[array_index]-ana_y[array_index]));
		}
		fprintf(fp,"\n");	
	}

	fclose(fp);

}					
