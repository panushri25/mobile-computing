//kernel
#pragma OPENCL EXTENSION cl_khr_fp64: enable

struct Node {
unsigned int numbranches ;
double *dicindx;
double *centroids;
double *dismat;
double *bounds;
struct Node *branches;
};

struct Reference{
int *input_vec;
unsigned int num;
struct Node noderef;
};

__kernel void product(__global double* h_inputvec,__global struct Reference* h_ref ,__global double* h_res,const int num_ipvec,const int N,const int numb,const int numb2)
  {	
	const int indicator =(get_global_id(0)) ;
	//printf(" size of h_ref = %d \n" , sizeof(h_ref));
	const int add = (int)(indicator*num_ipvec)/N;
	unsigned int ipv_size = h_ref[add].num ;
	
	//printf("indicator = %d sz = %d \n" , indicator, ipv_size);
	//printf("indicator = %d add = %d \n" , indicator, add);
						
	const int fract = (int)(ipv_size*N/num_ipvec); // fract = no. of compute units
	int fractn = (int) fract/numb;
	const int div = (int)ipv_size/fractn  ; // assume no.of vector greater than no. of fracts
						
	struct Node temp = h_ref[add].noderef ; // which node to refer and which input vec to work on 
	//printf("indicator = %d cent = %f \n" , indicator, temp.centroids[0]);
						
	double centroid[441];
	double value;
	int tempo ;
	
	for(int i=0; i<441 ; i++) // assuming that number of work items allotted are greater than centroids
	{
		int val =  (441*(indicator % numb));
		centroid[i] = temp.centroids[i + val] ;
	}	
	
	//printf( " centroid = %f  indicator = %d \n", centroid[0],indicator);
	
	for( int i = 0 ; i<= div ; i++)
	{			
		value = 0;
		//printf("i = %d  indicator = %d \n", i,indicator);
		int inc = (int) (indicator / numb);
		tempo = h_ref[add].input_vec[(i*fractn+inc)%ipv_size]; // can be improved
		
		//printf("centroid = %f input val = %f  tempo = %d element = %d indicator = %d \n" ,centroid[0],h_inputvec[441*tempo],tempo,(i*fractn+inc)%ipv_size,indicator);
		//printf("h_ref test 0 4 %d \n ", h_ref[0].input_vec[4]);
		//printf("h_ref test 1 2 %d \n ", h_ref[1].input_vec[2]);
		
		for(int j=0 ; j<441 ;j++)
		{						
			//printf("j = %d  indicator = %d \n", j,indicator);
			double diff = centroid[j]- h_inputvec[441*tempo + j];
			//printf("diff = %f , indicaror = %d \n",diff,indicator); 
			value += pow(diff,2);
		}
		double norm =  sqrt(value);
		h_res[numb2*tempo+(indicator%numb)] =  norm;	// check
		//printf("i = %d centroid = %f norm = %f , indicaror = %d tempo = %d location = %d element = %d \n",i,centroid[0],norm,indicator,tempo,N*tempo+(indicator%numb),(i*fractn+inc)%ipv_size); 
	}
																						
 }

__kernel void minimum(__global double* h_res,const int num_ipvec,__global int* h_minres,const int numb2,const int offset)
{	
	
	const int id =(get_global_id(0)) ;
	
	
			for(int i=0 ; i<num_ipvec ; i++)
			{
				float other = h_res[i*numb2+id+offset];
				float mine = h_res[i*numb2+id];
				if(mine <= other)
				{
					h_res[i*numb2+id] = mine ;
					
					//printf("min value at mine %f than other %f iteration %d indct %d  \n",h_res[i*numb2+id],other,i*numb2+id,id);
					//printf(" val at %d : %f numb2 %d  %d\n", i*numb2+id,h_res[i*numb2+id],numb2,i);
					
				}
				else
				{
					h_res[i*numb2+id] = other;
					int temp = (i*numb2+id);
					
					
						h_minres[temp] = h_minres[temp+ offset] ;
					//printf(" val at %d : %f \n", i*numb2+id,h_res[i*numb2+id]);
					
					//printf("min value at other %f than mine %f iteration %d indct %d  \n",h_res[i*numb2+id],mine,(i*numb2+id),id);
				}
				
			}
	
	
}
  
  
 

	
