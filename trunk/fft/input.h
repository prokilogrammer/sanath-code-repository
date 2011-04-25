#ifndef __INPUT__
#define __INPUT__

#include<vector>
#include<utility>
#include<map>
#include<cstdlib>
#include "polynomial.h"
using namespace std;

enum ELEMENTS{
C, 
H,
N, 
O, 
P, 
S, 
D, 
Cl,
F, 
Na,
K, 
Fe,
Cu,
Zn,
Se,
Br,
Li,
I,
Hg,
Si
};


struct isotope_st
{
	int pos;
	double prob;
};
vector<struct isotope_st *> isotope_prob[10];
int isotope_prob_count;

void populate_isotope_probability()
{
	struct isotope_st *temp;
	// C  2 12.0 .98893 13.003354838 .01107
	temp =  new struct isotope_st();
	temp->pos=0;
	temp->prob = .98893;
	isotope_prob[C].push_back(temp);

	temp =  new struct isotope_st();
	temp->pos=1;
	temp->prob = .01107;
	isotope_prob[C].push_back(temp);

	// H 2 1.007825032 .99985 2.014101778 0.00015
	temp =  new struct isotope_st();
	temp->pos=0;
	temp->prob = .99985;
	isotope_prob[H].push_back(temp);

	temp =  new struct isotope_st();
	temp->pos=1;
	temp->prob = 0.00015;
	isotope_prob[H].push_back(temp);

	// N 2 14.003074005 .996337 15.000108898 0.003663
	temp =  new struct isotope_st();
	temp->pos=0;
	temp->prob = .996337;
	isotope_prob[N].push_back(temp);

	temp =  new struct isotope_st();
	temp->pos=1;
	temp->prob = 0.003663;
	isotope_prob[N].push_back(temp);	

	// O 2 15.994914622 .99759 17.999160419 0.002036
	temp =  new struct isotope_st();
	temp->pos=0;
	temp->prob = .99759;
	isotope_prob[O].push_back(temp);

	temp =  new struct isotope_st();
	temp->pos=2;
	temp->prob = 0.002036;
	isotope_prob[O].push_back(temp);	

	// S 3 31.972070690 .9502 33.967866831 .0421 32.971458497 0.0075
	temp =  new struct isotope_st();
	temp->pos=0;
	temp->prob = .9502;
	isotope_prob[S].push_back(temp);

	temp =  new struct isotope_st();
	temp->pos=1;
	temp->prob = .0421;
	isotope_prob[S].push_back(temp);

	temp =  new struct isotope_st();
	temp->pos=1;
	temp->prob = 0.0075;
	isotope_prob[S].push_back(temp);				
		
}


typedef map<ELEMENTS,int> moleculeType;

struct amino_acid_st
{
	char name;
	moleculeType *composition;

}amino_acid_database[20];
int amino_acid_dbcount;

void populate_amino_acid()
{
	struct element_count *temp_ele;
	moleculeType * temp_aa;
	amino_acid_dbcount=0;
		
	// A
	// C3 H5 N1 O1
	temp_aa = new moleculeType();
	temp_aa->insert( pair<ELEMENTS,int>(C,3) );
	temp_aa->insert( pair<ELEMENTS,int>(H,5) );
	temp_aa->insert( pair<ELEMENTS,int>(N,1) );
	temp_aa->insert( pair<ELEMENTS,int>(O,1) );
	
	amino_acid_database[amino_acid_dbcount].name='A';
	amino_acid_database[amino_acid_dbcount].composition=temp_aa;
	amino_acid_dbcount++;
		
	// C
	// C3 H5 N1 O1 S1
	temp_aa = new moleculeType();
	temp_aa->insert( pair<ELEMENTS,int>(C,3) );
	temp_aa->insert( pair<ELEMENTS,int>(H,5) );
	temp_aa->insert( pair<ELEMENTS,int>(N,1) );
	temp_aa->insert( pair<ELEMENTS,int>(O,1) );
	temp_aa->insert( pair<ELEMENTS,int>(S,1) );

	amino_acid_database[amino_acid_dbcount].name='C';
	amino_acid_database[amino_acid_dbcount].composition=temp_aa;
	amino_acid_dbcount++;

	// D
	// C4 H5 N1 O3	
	temp_aa = new moleculeType();
	temp_aa->insert( pair<ELEMENTS,int>(C,4) );
	temp_aa->insert( pair<ELEMENTS,int>(H,5) );
	temp_aa->insert( pair<ELEMENTS,int>(N,1) );
	temp_aa->insert( pair<ELEMENTS,int>(O,3) );

	amino_acid_database[amino_acid_dbcount].name='D';
	amino_acid_database[amino_acid_dbcount].composition=temp_aa;
	amino_acid_dbcount++;
	
	// E
	// C5 H7 N1 O3
	temp_aa = new moleculeType();
	temp_aa->insert( pair<ELEMENTS,int>(C,5) );
	temp_aa->insert( pair<ELEMENTS,int>(H,7) );
	temp_aa->insert( pair<ELEMENTS,int>(N,1) );
	temp_aa->insert( pair<ELEMENTS,int>(O,3) );

	amino_acid_database[amino_acid_dbcount].name='E';
	amino_acid_database[amino_acid_dbcount].composition=temp_aa;
	amino_acid_dbcount++;		
	
	// L
	// C6 H11 N1 O1
	temp_aa = new moleculeType();
	temp_aa->insert( pair<ELEMENTS,int>(C,6) );
	temp_aa->insert( pair<ELEMENTS,int>(H,11) );
	temp_aa->insert( pair<ELEMENTS,int>(N,1) );
	temp_aa->insert( pair<ELEMENTS,int>(O,1) );

	amino_acid_database[amino_acid_dbcount].name='L';
	amino_acid_database[amino_acid_dbcount].composition=temp_aa;
	amino_acid_dbcount++;
			
	// M
	// C5 H9 N1 O1 S1
	temp_aa = new moleculeType();
	temp_aa->insert( pair<ELEMENTS,int>(C,5) );
	temp_aa->insert( pair<ELEMENTS,int>(H,9) );
	temp_aa->insert( pair<ELEMENTS,int>(N,1) );
	temp_aa->insert( pair<ELEMENTS,int>(O,1) );
	temp_aa->insert( pair<ELEMENTS,int>(S,1) );

	amino_acid_database[amino_acid_dbcount].name='M';
	amino_acid_database[amino_acid_dbcount].composition=temp_aa;
	amino_acid_dbcount++;
	
	// R
	// C6 H12 N4 O1
	temp_aa = new moleculeType();
	temp_aa->insert( pair<ELEMENTS,int>(C,6) );
	temp_aa->insert( pair<ELEMENTS,int>(H,12) );
	temp_aa->insert( pair<ELEMENTS,int>(N,4) );
	temp_aa->insert( pair<ELEMENTS,int>(O,1) );

	amino_acid_database[amino_acid_dbcount].name='R';
	amino_acid_database[amino_acid_dbcount].composition=temp_aa;
	amino_acid_dbcount++;

	// S
	// C3 H5 N1 O2
	temp_aa = new moleculeType();
	temp_aa->insert( pair<ELEMENTS,int>(C,3) );
	temp_aa->insert( pair<ELEMENTS,int>(H,5) );
	temp_aa->insert( pair<ELEMENTS,int>(N,1) );
	temp_aa->insert( pair<ELEMENTS,int>(O,2) );

	amino_acid_database[amino_acid_dbcount].name='S';
	amino_acid_database[amino_acid_dbcount].composition=temp_aa;
	amino_acid_dbcount++;	
	
	/*
	F
	C9 H9 N1 O1
	G
	C2 H3 N1 O1
	H
	C6 H7 N3 O1
	I
	C6 H11 N1 O1
	K
	C6 H12 N2 O1
	N
	C4 H6 N2 O2
	P
	C5 H7 N1 O1
	Q
	C5 H8 N2 O2
	T
	C4 H7 N1 O2
	V
	C5 H9 N1 O1
	W
	C11 H10 N2 O1
	Y
	C9 H9 N1 O2
	h
	C4 H5 N O
	m
	C5 H9 N O2 S
	s
	C3 H6 N1 O5 P
	t
	C4 H8 N1 O5 P
	y
	C9 H10 N1 O5 P
	U
	C3 H5 N O Se
	*/
}

moleculeType *find_aminoacid(char c)
{
	int i;

	for(i=0;i<amino_acid_dbcount;i++)
	{
		if( amino_acid_database[i].name == c)
			return amino_acid_database[i].composition;
	}
	
	return NULL;
	
}

int populate_peptide_composition(string peptide, moleculeType *peptide_composition)
{
	moleculeType *amino;
	moleculeType::iterator it;
	int count=0;
	int i,j;
	for(i=0;i<peptide.length();i++)
	{
		amino = find_aminoacid(peptide[i]); //will return the molecular composition of peptide
		
		if( amino== NULL)
		{
			printf("ERROR: Wrong Amino Acid \n");exit(1);
		}
			
		for ( it=amino->begin() ; it != amino->end(); it++ )
		{
			(*peptide_composition)[  it->first ] += it->second;
			count += it->second * (isotope_prob[it->first].size()+1); //just for an upper bound on polynomial size
		}

	}

// for deuterium exchange
//	(*peptide_composition)[  H ] = (*peptide_composition)[H] - 0.60*(*peptide_composition)[H];

/*
	for(it=(*peptide_composition).begin(); it != (*peptide_composition).end(); it++)
	{
		printf("%d %d\n", it->first, it->second);
	}	
	*/
	return count;
}
#endif // __INPUT__
