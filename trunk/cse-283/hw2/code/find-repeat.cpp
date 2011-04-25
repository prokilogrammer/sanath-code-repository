#include<cstdio>
#include<vector>
#include<bitset>
#include<cmath>
#include<cstdlib>
#include<cstring>
#include<map>
#include<ctime>
#include<iostream>

using namespace std;
#define BIT char

struct repeat_list_st
{
    BIT sequence[17];
    vector<int> pos;
};

vector<struct repeat_list_st *> output;

struct trie_st
{
    // Its assumed that next[] should be
    // initialized to zero
    // a - g - t - c
    struct trie_st* next[4];
  //  struct repeat_list_st* repeat_list_ptr;
};
struct trie_st *trie_root;

vector<void*> allocated_tries;

char** read_data(string filename,int length, char** filedata)
{

    
    int i,j;
    char c;
    FILE *fp = fopen(filename.c_str(),"r");
    if(fp == NULL)
        perror("Unable to open file\n");
    for(i=0;i<10;i++)
    {
        filedata[i] = (char*)malloc(length+1*sizeof(char));
    }

    
    for(i=0;i<10;i++)
    {
        j=0;
        while(1)
        {   
            c = fgetc(fp);
            if(c == '\n')
                break;
            filedata[i][j++] = c;  

        }    
        filedata[i][j]='\0';
    }
    fclose(fp);
    return filedata;   
}

void freememory()
{
	vector<struct repeat_list_st *>::iterator it;
	for(it = output.begin(); it < output.end(); it++)
	{
		(*it)->pos.clear();
		free( (*it) );
	}
	output.clear();

	vector<void*>::iterator itt;
	for(itt = allocated_tries.begin(); itt < allocated_tries.end(); itt++)
	{
		free( (*itt) );
	
	}
	allocated_tries.clear();
	trie_root->next[0]=NULL;
	trie_root->next[1]=NULL;
	trie_root->next[2]=NULL;
	trie_root->next[3]=NULL;
	//trie_root->repeat_list_ptr=NULL;
}

void construct_and_find(char* data, int length)
{
    
    int start,i, index;
    struct trie_st* current_node;
    clock_t start_time = clock();

    for(start=0; start <=(length-16); start++)
    {
        BIT sequence[17]= {0};
        current_node = trie_root;
        for(i=0;i<16;i++)
        {
            switch(data[i+start])
            {
                case 'a':
                    index = 0;
                    sequence[i] = 'a';
                    break;
                case 'g':
                    index=1;
                    sequence[i] = 'g';
                    break;
                case 't':
                    index=2;
                    sequence[i] = 't';
                    break;
                case 'c':
                    index=3;
                    sequence[i] = 'c';
                    break;
                default:
                    printf("***** ERROR IN DNA *****\n");
            }
            if( current_node->next[index] == NULL )
            {
                current_node->next[index] = (struct trie_st*)calloc(1, sizeof(struct trie_st)); 
                allocated_tries.push_back(current_node->next[index]);
            }
            current_node = current_node->next[index];
        }
        if( i==16)
        {//reached end of a string
        	// next[0] will hold the repeat_list_st
        	repeat_list_st *tempptr;
            if(current_node->next[0] == NULL)
            {
                // first time i'm encountering this string; create new repeat_list
                // this is going to be a play with pointers. I'll just store my
                // pointer to repeat_list_st inside next[0] but later type cast it
                // to repeat_list_st when using.(Memory sainvg trick)
                current_node->next[0] = (struct trie_st*)calloc(1, sizeof(struct repeat_list_st));
                output.push_back((struct repeat_list_st*)(current_node->next[0]));
                
            }
            tempptr = (struct repeat_list_st*)current_node->next[0];
            memcpy( tempptr->sequence, sequence, sizeof(BIT)*17);
            tempptr->pos.push_back(start);
        }
        
    }

     clock_t ends = clock();
     cout << "Running Time : " << (double) (ends - start_time) / CLOCKS_PER_SEC << endl;
/*
	vector<struct repeat_list_st *>::iterator it;
	for(it = output.begin(); it < output.end(); it++)
	{
		vector<int>::iterator itt;
		if((*it)->pos.size() > 1)
		{
			printf("%s:\n",(*it)->sequence);
			for(itt = (*it)->pos.begin(); itt< (*it)->pos.end(); itt++)
			{
				printf("\t%d\n",*itt);		
				//(*it)->pos.erase(itt);
			}
		}
		//printf("%d\n",(int)(*it)->pos.size());
		
	}*/
	
	freememory();
}

int main(int argc, char** argv)
{
    char* filedata[10];
    int i, length;
    length=pow(4,7);
    read_data("dna-4-7.txt",length,filedata);
     
    trie_root = (struct trie_st*)calloc(1, sizeof(struct trie_st)); 
    // main loop 
//    for(i=0;i<10;i++)
 //   {
        if(argc==1) perror("dog");
       i=atoi(argv[1]);
    	//printf("Computing for test case #%d\n",i);
        construct_and_find(filedata[i],length);
		free(filedata[i]);        
//    }
    return 0;
}

