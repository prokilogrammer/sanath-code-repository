#ifndef __POLYNOMIAL__
#define __POLYNOMIAL__
using namespace std;

#define PI 3.14159265

typedef complex<double> dcomplex;

class polynomial
{
	public:
		dcomplex *data;
		unsigned int fft_length;
		unsigned int power2_length;
		unsigned int data_length;

	public:
		polynomial(int length)
		{
			data_length=0;
			fft_length=0;
			power2_length = 0; 

	        data = new dcomplex[next_powerof2(length)];
		}
		
		void adjust_data(int data_len, int fft_len)
		{
			data_length = data_len;
			fft_length = fft_len;
			power2_length = next_powerof2(fft_length);
			//printf("%d %d %d\n",data_length, fft_length, power2_length);
			clear_higher_coef();
			adjust_for_powerof2();
		}
		
		void do_fft()
		{
			bit_reverse_copy();
			fft(false);
		}
		void do_ifft()
		{	
			bit_reverse_copy();
			fft(true);    
		}
		
		
		void generate_data(int length)
		{
			int i;
			for(i=0;i<length;i++)
			{
				data[i] = i+1;
			}
		}

		void clear_data()
		{
			int i;
			for(i=0;i<power2_length;i++)
				data[i]=0;
		}
		
		void clear_higher_coef()
		{
			int i;
			for(i=data_length; i<fft_length;i++)
			{
				data[i]=0;
			}

		}
		void adjust_for_powerof2()
		{
			// Doing Mirroring as suggested in Cormen
			int i;
			for(i=0; i<(power2_length - fft_length); i++)
			{
				data[fft_length+i] = data[fft_length-i-2]; // this can be zero also. 

			}
		}


		~polynomial()
		{
			delete[] data;
		}
		
	private:
		int floorLog2(unsigned int n) 
		{
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


		int next_powerof2(int length)
		{
			int nbits = floorLog2(length);
			if( length == (1<<nbits))
				return length;
			else
				return (1 << (nbits+1));

		}

		
		void fft(bool flag)
		{
			//variables just for convenience
			dcomplex* A;
			unsigned  int n;
			A = data;
			n = power2_length;
			
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

		void bit_reverse_copy()
		{
			 unsigned int n = power2_length;

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

		//		printf("%d %d\n",i,rev);
				dcomplex temp = data[rev];
				data[rev] = data[i];
				data[i] = temp;
				reversed[i]=true;
				reversed[rev]=true;
				
			}
		}
};

#endif //__POLYNOMIAL__
