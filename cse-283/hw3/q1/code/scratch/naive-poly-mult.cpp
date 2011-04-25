#include<stdio.h>


void generate_data(double *data, int length)
{

	int i;
	for(i=0;i<length;i++)
	{
		data[i] = i+1;
	}

}

#define LENGTH 2
int main(int argc, char** argv)
{
	double *A = new double[2*LENGTH];
	double *B = new double[2*LENGTH];
	double *C = new double[2*LENGTH];

//	generate_data(A,LENGTH);
//	generate_data(B,LENGTH);
	A[0] = 0.98893;
	A[1] = 0.01107;

	B[0] = 0.98893;
	B[1] = 0.01107;
	
	int i,j;
	for(i=0;i<LENGTH;i++)
	{
		for(j=0;j<LENGTH;j++)
		{
			C[i+j] += A[i]*B[j];
		}

	}

	for(i=0;i<2*LENGTH;i++)
		printf("%lf\n",C[i]);
		
	return 0;
}
