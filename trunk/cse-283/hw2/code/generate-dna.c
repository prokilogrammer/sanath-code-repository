#include<stdio.h>
#include<math.h>
#include<stdlib.h>


void generate(int length, char* filename)
{
    int i,j,k;
    
    FILE *fp = fopen(filename,"w");
    for(j=0;j<10;j++)
    {
        for(k=length; k>0;k--)
        {
            i=rand()%4;
            if( i == 0 )
                fputc( 'a', fp);
            else if( i == 1 )
                fputc('g',fp);
            else if( i == 2 )  
                fputc('t',fp);
            else
                fputc('c',fp);
            
        }
        fputc('\n',fp);
    }
    fclose(fp);

}

int main(int argc, char** argv)
{

    srand( time(NULL) );
    generate(pow(4,7), "dna-4-7.txt");
    generate(pow(4,8), "dna-4-8.txt");
    generate(pow(4,9), "dna-4-9.txt");
    generate(pow(4,10), "dna-4-10.txt");
    generate(pow(4,11), "dna-4-11.txt");
    generate(pow(4,12), "dna-4-12.txt");
    
}

