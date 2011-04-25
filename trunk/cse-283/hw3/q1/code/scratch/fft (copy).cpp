#include<complex>
#include<iostream>
#include<cstdio>

using namespace std;

#define PI 3.14

typedef complex<double> dcomplex;

int floorLog2(unsigned int n) {
// Function courtesy: Wikipedia
  if (n == 0)
    return -1;
 
  int pos = 0;
  if (n >= 1<<16) { n >>= 16; pos += 16; }
  if (n >= 1<< 8) { n >>=  8; pos +=  8; }
  if (n >= 1<< 4) { n >>=  4; pos +=  4; }
  if (n >= 1<< 2) { n >>=  2; pos +=  2; }
  if (n >= 1<< 1) {           pos +=  1; }
  return pos;
}

void ifft(dcomplex* A,unsigned  int n)
{

	unsigned int s,k,j;
	unsigned int m;
	dcomplex omegaM, omega, t, u;
	for(s=1; s<= floorLog2(n) ; s++)
	{
		m = 1<<s;
		omegaM = exp(  dcomplex(0, -2*PI/m) );

		for(k=0; k < (n-1)/m; k++)
		{
			omega = 1;
			for(j=0; j< (m/2); j++)
			{
				t = omega * A[k+j+m/2];
				u = A[k+j];

				A[k+j]= u+t;
				A[k+j+m/2] = u-t;

				omega = omega * omegaM;
			}

		}
	}

}

void fft(dcomplex* A,unsigned  int n, bool flag)
{
	unsigned int s,k,j;
	unsigned int m;
	dcomplex omegaM, omega, t, u;
	for(s=1; s<= floorLog2(n) ; s++)
	{
		m = 1<<s;
		if(flag == false)
			omegaM = exp(  dcomplex(0, 2*PI/m) );
		else
			omegaM = exp(  dcomplex(0, -2*PI/m) );	

		for(k=0; k <(n-1); k+=m)
		{
			omega = 1;
			for(j=0; j< (m/2); j++)
			{
				t = omega * A[k+j+m/2];
				u = A[k+j];

				A[k+j]= u+t;
				A[k+j+m/2] = u-t;

				omega = omega * omegaM;
			}

		}
	}
}

void bit_reverse_copy(dcomplex* input, dcomplex* output, unsigned int n)
{
    
    unsigned int rev,act;
    unsigned int i;
    int nbits = floorLog2(n);
 //   printf("%d\n",nbits);
    int j = nbits;
	unsigned int N = 1<<nbits;

	bool *reversed = new bool[n];
	
    for(i=0;i<n;i++)
    {
    	if(reversed[i]==true)
    		continue;
    		
	    nbits = j-1;
		act = rev = i;
    	while(nbits>0)
    	{
			act = act>>1;
			rev = rev<<1;
			rev = rev | (act&1);
			nbits--;
		}
		rev = rev & (N-1);

		printf("%d %d\n",i,rev);

//		output[rev] = input[i];

		dcomplex temp = input[rev];
		input[rev] = input[i];
		input[i] = temp;
		reversed[i]=true;
		reversed[rev]=true;
    }
}

void generate_input(dcomplex *input, int length)
{

	int i;
	for(i=0;i<length;i++)
	{
		input[i] = i;
	}

}

int next_powerof2(int length)
{
	int nbits = floorLog2(length);
	printf("%d\n",nbits);
	if( length == (1<<nbits))
		return length;
	else
		return (1 << (nbits+1));

}


void adjust_for_powerof2(dcomplex *input, int old_length, int new_length)
{
	// Doing Mirroring as suggested in Cormen
	int i;
	for(i=0; i<(new_length-old_length); i++)
	{
		input[old_length+i] =  input[old_length-i-2];

	}
}

void double_polynomial_data(dcomplex *input, int old_length)
{

	int i;
	for(i=old_length; i<(2*old_length); i++)
	{
		// Make the higher coefficients zero.
		input[i] = 0;
	}

}


#define LENGTH 8

int main(int argc, char** argv)
{

// orig_length = length of actual polynomial.
// data_length = double the size of input polynomial
//               coz after multiplying output will
//               have double length as input.
// power2_length = Adjust data_length to be powers of two
//                 coz FFT needs data to be powers of two.

    unsigned int orig_length=LENGTH;
    unsigned int data_length = orig_length*2;
    // We allocate double size polynomials because
    // after multiplication, we get double size
    // output
    unsigned int power2_length = next_powerof2(data_length);
    printf("%d\n",power2_length);
    
    dcomplex *input = new dcomplex[power2_length];
    dcomplex *interm = new dcomplex[power2_length];
    dcomplex *output = new dcomplex[power2_length];

	generate_input(input,orig_length);
	double_polynomial_data(input,orig_length);
	adjust_for_powerof2(input,data_length,power2_length);

    
    bit_reverse_copy(input, interm,power2_length);
    fft(input,power2_length,false);
    
    bit_reverse_copy(input, output,power2_length);
    fft(input,power2_length,true);    


	int i;
	for(i=0;i<power2_length;i++)
	{
		//printf("%lf, %lf, %lf\n",output[i].real()/power2_length,output[i].imag()/power2_length, abs(output[i])/power2_length);
		printf("%lf, %lf, %lf\n",input[i].real()/power2_length,input[i].imag()/power2_length, abs(input[i])/power2_length);
	}

	delete[] input;
	delete[] interm;
	delete[] output;
    //printf("%d\n",ceilLog2(7));
    return 0;
}
