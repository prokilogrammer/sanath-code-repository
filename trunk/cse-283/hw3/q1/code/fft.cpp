#include<complex>
#include<iostream>
#include<cstdio>

#include "polynomial.h"
#include "input.h"

using namespace std;

void multiply_polynomial(polynomial &A, polynomial &B, int a_length, int b_length, int fft_length)
{

	A.adjust_data(a_length, fft_length);
	A.do_fft();


	B.adjust_data(b_length, fft_length);
	B.do_fft();

	// Multiplying the FFT values	
	int i;
	for(i=0;i<A.power2_length;i++)
	{
		A.data[i] = A.data[i]*B.data[i];
	}
	
	A.do_ifft();	

	// Dividing ifft output by length to get actual value
	for(i=0;i<A.power2_length;i++)
	{
		A.data[i] /= A.power2_length; 
	}	
}

void print_output(polynomial &C)
{
	int i;
	printf("\n\n");
	printf("Isotope Profile\n");
	printf("-----------------------------\n");
	for(i=0;i<C.fft_length;i++)
	{
		printf("%d  %0.15lf\n",i, C.data[i].real()) ;//,C.data[i].imag(),  abs(C.data[i])  );
	}
	printf("\n\n");
}

void generate_multiply_polynomial(int length, moleculeType peptide_composition)
{
	polynomial A(length);
	polynomial B(length);

	int a_length=0;
	int b_length=0;
	int fft_length;
	bool flag=false;
	int max=100000; //some large number
	
	vector<struct isotope_st *>::iterator isoit;
	moleculeType::iterator it;
	for(it = peptide_composition.begin(); it != peptide_composition.end(); it++)
	{
		if(flag == false)
		{
			//executed only once.

			// populate A
			for(isoit = (isotope_prob[it->first]).begin(); isoit != (isotope_prob[it->first]).end(); isoit++)
			{
				A.data[(*isoit)->pos] = (*isoit)->prob;
				if(max==100000)
					max = (*isoit)->pos;
				if(max < (*isoit)->pos)
					max = (*isoit)->pos;
			}
			a_length = max+1;
			it->second--;
			flag=true;
		}
		
		for( ; it->second >0; it->second--)
		{
			max = 100000;
			// populate B
			for(isoit = isotope_prob[it->first].begin(); isoit != isotope_prob[it->first].end(); isoit++)
			{
				B.data[(*isoit)->pos] = (*isoit)->prob;
				if(max==100000)
					max = (*isoit)->pos;
				if(max < (*isoit)->pos)
					max = (*isoit)->pos;
			}
			b_length = max+1;
			
			multiply_polynomial(A, B,  a_length, b_length,  a_length + b_length - 1);

			// A holds the output. So it shouldn't be cleared
			B.clear_data();	
			a_length = a_length + b_length -1;			

		}
	}
print_output(A);

}

#define LENGTH 7
int main(int argc, char** argv)
{
	moleculeType peptide_composition;

	populate_isotope_probability(); //this must come first
	populate_amino_acid();
	int max_len = populate_peptide_composition("SLAMMER", &peptide_composition);
//	printf("max_len  %d\n",max_len );

	generate_multiply_polynomial(max_len, peptide_composition);
	
/*	polynomial A(2*LENGTH);
	polynomial B(2*LENGTH);

	A.generate_data(LENGTH);
	B.generate_data(LENGTH);
	
	multiply_polynomial(A, B, LENGTH, LENGTH, 2*LENGTH-1);		
*/		


    return 0;
}
