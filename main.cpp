// includes

#include <vector>
#include <matrix.h>
#include <mat.h>
#include <mex.h>
#include <iostream>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <fstream>
#include <stdlib.h>
#include <mpi.h>
#include <sys/stat.h>
#include <time.h>
#include <limits>
#include <string.h>

#ifdef __APPLE__

#include <OpenCL/OpenCl.h>

#else

#include <CL/cl.h>

#endif

using namespace std;

struct Node {
public:
unsigned int numbranches ;
double *dicindx;
double *centroids;
double *dismat;
double *bounds;
Node *branches;
};

int cnt = 0 ;
Node parent;
const int num_level = 3;
const int num_ipvec =  500;
int num_dict = 1000;
double inputs[441*num_ipvec];
unsigned int num_clusters[num_level] = {10,10,10};
cl_uint N = 5024;
int Nt = 5024;

char *inputfile = "ainput500.mat" ;
char *dictfile = "10x10x10.mat" ;

int correct  ;
double error2;

int store[num_level] = {1000,100,10};
int lvlcnt = 0;

void setsize(unsigned int st,unsigned int size, Node &temp)
{
	unsigned int numbranchs = size;
	printf(" test st test = %u ",st);
	printf(" test numbranches = %u " ,size);
	printf(" test address = %p ///////////////////////////////////////////////////////////////////////////////////////////// \n",&temp.numbranches);
	
	temp.numbranches = numbranchs;
	temp.dicindx = new double[st];   
	temp.centroids = new double[441*numbranchs]; 
	temp.dismat = new double[numbranchs*numbranchs];
	temp.bounds = new double[numbranchs];	  
	temp.branches = new Node[numbranchs];
}
  				
mxClassID   analyze_class(const mxArray *array_ptr , Node &temp);

void analyze_structure(const mxArray *structure_array_ptr, Node &temp)
{
	mwSize total_num_of_elements;
	mwIndex index;
	int number_of_fields, field_index;
	const char  *field_name;
	const mxArray *field_array_ptr;
  
	cnt=0;  
	total_num_of_elements = mxGetNumberOfElements(structure_array_ptr); 
	number_of_fields = mxGetNumberOfFields(structure_array_ptr);
  
	/* Walk through each structure element. */
	for (index=0; index<total_num_of_elements; index++)  {
    
		/* For the given index, walk through each field. */ 
		for (field_index=0; field_index<number_of_fields; field_index++)  {
			mexPrintf("\n\t\t");
			// display_subscript(structure_array_ptr, index);
			field_name = mxGetFieldNameByNumber(structure_array_ptr, field_index);
			mexPrintf(".%s\n", field_name);
			field_array_ptr = mxGetFieldByNumber(structure_array_ptr,index, field_index);
			if (field_array_ptr == NULL)
			{
				mexPrintf("\tEmpty Field\n");
			}	 
			else 
			{
				/* Display a top banner. */
				mexPrintf("------------------------------------------------\n");
				//  get_characteristics(field_array_ptr);
				analyze_class(field_array_ptr ,temp);
				mexPrintf("\n");
			}
		}
    mexPrintf("\n\n");
	} 
}


void analyze_double(const mxArray *array_ptr , Node &temp)
{
	double *pr,*pi,*temp_elm; 
	mwSize total_num_of_elements, index; 
  
	pr = mxGetPr(array_ptr);
	pi = mxGetPi(array_ptr);
	total_num_of_elements = mxGetNumberOfElements(array_ptr);
  
	cnt = cnt + 1;
	mexPrintf ("cnt = %d \n",cnt);  
	mexPrintf ("total_num_of_elements = %d \n",total_num_of_elements);

	switch(cnt)
	{
		case 1 : {setsize(store[lvlcnt],num_clusters[lvlcnt],temp);printf("num = %u \n",temp.numbranches); break;}
		case 2 : {temp_elm = temp.dicindx; break;}
		case 3 : {temp_elm = temp.centroids; break;}
		case 4 : {temp_elm = temp.dismat; break;}
		case 5 : { temp_elm = temp.bounds; break;}
		case 10 : {temp_elm = inputs ; break ; }
		default : {cnt=0;break;}
	}
   
	for (index=0; index<total_num_of_elements; index++)  
	{
		//mexPrintf(" %d = %g\n",index, *pr);
		if(cnt != 1)
		{
			*temp_elm = *pr ;
			temp_elm = temp_elm+1;
			pr=pr+1;
		}
	} 
}

void analyze_full(const mxArray *numeric_array_ptr, Node &temp)
{
	mxClassID   category;
	category = mxGetClassID(numeric_array_ptr);
	switch (category)  
	{ 
		case mxDOUBLE_CLASS: analyze_double(numeric_array_ptr, temp); break;
		default: break;
	}
}

/* Pass analyze_cell a pointer to a cell mxArray.  Each element
   in a cell mxArray is called a "cell"; each cell holds zero
   or one mxArray.  analyze_cell accesses each cell and displays
   information about it. */  
 void analyze_cell(const mxArray *cell_array_ptr , Node &temp)
{
	mwSize total_num_of_cells;
	mwIndex index;
	const mxArray *cell_element_ptr;
  
	total_num_of_cells = mxGetNumberOfElements(cell_array_ptr); 
	mexPrintf("total num of cells = %d\n", total_num_of_cells);
	mexPrintf("******************************************************************************************************* \n");
	int temp_cluster_index ;
	int temp_store ;
	/* Each cell mxArray contains m-by-n cells; Each of these cells
     is an mxArray. */ 
	for (index=0; index<total_num_of_cells; index++) 
	{	
		
		if(index ==0)
		{
			lvlcnt += 1;
			printf( " inc cluster index %d  total_num_of_cells %d \n" , lvlcnt,total_num_of_cells);
			//cluster_index += 1 ;
			//store[lvlcnt] = int(store / temp.numbranches) ;
		}
		
	
		mexPrintf("\n\n\t\tCell Element: ");
		//display_subscript(cell_array_ptr, index);
		mexPrintf("\n");
		cell_element_ptr = mxGetCell(cell_array_ptr, index);
		if (cell_element_ptr == NULL) 
		{
			mexPrintf("\tEmpty Cell\n");
		}
		else {
		/* Display a top banner. */
		mexPrintf("------------------------------------------------\n");
		//get_characteristics(cell_element_ptr);
		printf(" test address = %p ",&temp.branches[index]);
		analyze_class(cell_element_ptr,temp.branches[index]);
		mexPrintf("\n");
		}
		if(index == total_num_of_cells - 1)
			{
				lvlcnt = lvlcnt - 1;
			}
	}

	mexPrintf("\n");
}

/* Determine the category (class) of the input array_ptr, and then
   branch to the appropriate analysis routine. */
mxClassID analyze_class(const mxArray *array_ptr, Node &temp)
{
    mxClassID  category;
    
    category = mxGetClassID(array_ptr);
    
    switch (category)
	{
		case mxSTRUCT_CLASS:  analyze_structure(array_ptr , temp); break;
        case mxCELL_CLASS:    analyze_cell(array_ptr , temp); break;
        case mxUNKNOWN_CLASS: mexPrintf("Unknown class."); break;
        default:              analyze_full(array_ptr , temp ); break;
    }
    return(category);
}


void pfn_notify(const char *errinfo, const void *private_info, size_t cb, void *user_data) // what is this for ??????
{
	fprintf(stderr, "OpenCL Error (via pfn_notify): %s\n", errinfo);
}

void oclInit(int plat, int dev,cl_context &context,cl_device_id &device,cl_command_queue &queue)
{
	/* set up CL */
	cl_int            err; 
	cl_platform_id    platforms[100];
	cl_uint           platforms_n;
	cl_device_id      devices[100];
	cl_uint           devices_n ;
	
	/* get list of platform IDs (platform == implementation of OpenCL) */

	err = clGetPlatformIDs(100, platforms, &platforms_n); //' Get the platform ID , platform_n number returned

	printf("platforms_n = %d \n ", platforms_n );
	printf("err = %d ", err );
	
	if( plat > platforms_n) 
	{
		printf("platforms_n %d \n ", platforms_n );
		printf("ERROR: platform %d unavailable \n", plat);
		exit(-1); 
	}

	// find all available device IDs on chosen platform (could restrict to CPU or GPU)

	cl_uint dtype = CL_DEVICE_TYPE_ALL;

	//clGetDeviceInfo(dev, CL_DEVICE_NAME, sizeof(buffer), buffer, NULL);
    // std::cout << "  DEVICE_NAME = " << buffer << std::endl;

	clGetDeviceIDs( platforms[plat], dtype, 100, devices, &devices_n);
	printf("devices_n = %d\n", devices_n);

	if(dev>=devices_n)
	{
		printf("invalid device number for this platform\n");
		exit(0);
	}
	
	// choose user specified device
	device = devices[dev];

	// make compute context on device
	context = clCreateContext((cl_context_properties *)NULL, 1, &device, &pfn_notify, (void*)NULL, &err); // pfn_notify , arguments ?????

	// create command queue
	queue   = clCreateCommandQueue(context, device, CL_QUEUE_PROFILING_ENABLE, &err);
}


void oclBuildKernel(const char *sourceFileName,
    const char *functionName,
    cl_context &context,
    cl_device_id &device,
    cl_kernel &kernel,
    const char *flags
    ){

  // read in text from source file
		cl_int    err;
		struct stat statbuf;  // stat statbuf ??????
  FILE *fh = fopen(sourceFileName, "r");
  if (fh == 0){
    	printf("Failed to open: %s\n", sourceFileName);
    	throw 1;
 }

  
  /* create program from source */

  /* get stats for source file */

  stat(sourceFileName, &statbuf); // ???????????????

  /* read text from source file and add terminator */

  char *source = (char *) malloc(statbuf.st_size + 1);
  fread(source, statbuf.st_size, 1, fh);
  source[statbuf.st_size] = '\0';

  /* create program from source */

  cl_program program = clCreateProgramWithSource(context, 1, (const char **) & source, (size_t*) NULL, &err);
  if (!program){
    printf("Error: Failed to create compute program!\n");
    throw 1;
  }

  /* compile and build program */
  err = clBuildProgram(program, 1, &device, flags, (void (*)(cl_program, void*))  NULL, NULL);
  /* check for compilation errors */
  char *build_log;
  size_t ret_val_size;
  err = clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, 0, NULL, &ret_val_size);
  build_log = (char*) malloc(ret_val_size+1);
  err = clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, ret_val_size, build_log, (size_t*) NULL);

  /* to be carefully, terminate with \0

     there's no information in the reference whether the string is 0 terminated or not */

  build_log[ret_val_size] = '\0';

  /* print out compilation log */
  fprintf(stderr, "%s", build_log );
  /* create runnable kernel */
  kernel = clCreateKernel(program, functionName, &err);
  if (! kernel || err != CL_SUCCESS){
    printf("Error: Failed to create compute kernel!\n");
    throw 1;
  }

}

struct Reference{
public:
cl_uint* input_vec;
cl_uint num;
Node noderef;
};

Reference *ref = new Reference[num_ipvec];
Reference temp_ref[num_ipvec];

int main(int argc, char** argv){

	float avgTimePerIteration = 0;
	float diff = 0;
	int iter = 0;
	int num_iter = 50;
	cl_uint node_taken[num_level+1][num_ipvec]; 
	MPI_Init(&argc, &argv);
  
for( iter = 0; iter < num_iter ; iter++)
{	
	// read the .mat file output of the k-means balanced clustered MATLAB code 
	
	
	MATFile *datafile;
	datafile = matOpen( dictfile, "r");
	
    mxArray  *structure_array_ptr ;
    structure_array_ptr = matGetVariable(datafile, "streeClus");
    analyze_structure(structure_array_ptr , parent ); 
    matClose(datafile);
	
	cnt = 9;
	cout << " done dict reading mat \n";

	// read the input file : here the input file is the dictionary elements itself
	//  trying to check the accuracy of algoritm
	
	MATFile *datafile1; 
	datafile1 = matOpen( inputfile, "r");
	mxArray  *input_array_ptr ;
    input_array_ptr = matGetVariable(datafile1, "input");
    analyze_double(input_array_ptr , parent );
    matClose(datafile1);	
	
	//cout <<"test 1 "<< parent.centroids[0]<<endl;
	//cout <<"test 2 "<< parent.branches[1].centroids[0]<<endl;
	//cout <<"test 3 "<< parent.branches[1].branches[3].centroids[0]<<endl;
	cout << " done input reading mat \n";
	
	Node *level_node[num_level]; // all nodes in each level ;
			
	int temp_store = 1;
	level_node[0] = new Node[temp_store]; 
	level_node[0][0] = parent ;
			
	//cout<<"test this? "<<level_node[0][0].dicindx[3]<<"\n";
	 
	// copy to level_node for easy access in subsequent code 
	// one time operation wont effect the runtime

	for(int i =1 ; i<num_level ; i++)
	{
		printf("i= %d \n",i);
		int save = temp_store ;
		temp_store = temp_store * level_node[i-1][0].numbranches ;
		level_node[i] = new Node[temp_store] ;
		for( int po = 0 ; po< save ; po++)
		{
			for (unsigned int t = 0 ; t < level_node[i-1][0].numbranches ; t++ )
			{
				printf("po = %d \n",po*(level_node[i-1][0].numbranches)+t);
				level_node[i][po*(level_node[i-1][0].numbranches)+t]=(level_node[i-1][po].branches[t]); 
			}
		}
	}
	
	
	cl_int  err;
	int plat = 0;
	int dev  = 0;
	cl_context context;
	cl_device_id device;
	cl_command_queue queue;
	cl_kernel kernel;
	cl_kernel kernel1;
	cl_event event_time; 

	// initialize platform
	oclInit(plat, dev, context, device, queue);

	const char *sourceFileName = "new_kernel1.cl";
	const char *functionName = "product";
	const char *functionName1 = "minimum";
	int BDIM = 16;
	char flags[BUFSIZ];
	sprintf(flags, "-DBDIM=%d", BDIM); // what is this ?????????

	// bulid kernel : check for syntax errors
	
	oclBuildKernel(sourceFileName,
	functionName,
	context,
	device,
	kernel,
	flags);
	
	oclBuildKernel(sourceFileName,
	functionName1,
	context,
	device,
	kernel1,
	flags);
			
  	// Determine maximum number of compute units available
  	//clGetDeviceInfo(device, CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(N), &N, NULL);
    //cout << "  DEVICE_MAX_COMPUTE_UNITS = " << (unsigned int)N << endl;

	/* set thread array */
 	int dim = 1;
	
 	
	// define local and global dimensions
  	size_t local[1] = {Nt};
  	size_t global[1] = {N};
  
	// store the cluster node taken by the input vector at each level of the hierarchy 
			
					
	// Initialize the Reference hash memory				
	for(int i=0 ; i<num_ipvec ; i++)
	{
		ref[i].noderef = level_node[0][0] ;
		ref[i].input_vec = (cl_uint*)malloc(num_ipvec*sizeof(int));
		for(int j=0 ; j<num_ipvec ; j++)
		{
			ref[i].input_vec[j]=j; // set size of ref
			node_taken[0][j] = 0;					
		}
		ref[i].num = num_ipvec ;							
	}
			
			
	printf("start loop \n");
	int id = 1;
	
	// traverse all the levels in the hierarchical tree structure
	double startTime = MPI_Wtime();
	for(int level=1 ; level <= (num_level); level++)
	{ 			
		printf( "level = %d  ****************************************************************************************** \n",level); 
															
		int numb = num_clusters[level-1];	
		id = id * numb ; // number of clusters in the level 								

		int pops = ceil(log2(numb));
		int numb2 = pow(2,pops);
		*temp_ref= *ref ; // can optimize this loop
		
				
		size_t sz_ipv = 441*num_ipvec*sizeof(double);
		size_t sz_ref = sizeof(temp_ref);
		size_t  sz_res = numb2*num_ipvec*sizeof(double); 
		size_t sz_minres = num_ipvec*sizeof(int)*(numb2);
		
		/* create host array */
		
	 	double *h_inputvec = (double*) calloc(441*num_ipvec,sizeof(double));
	 	Reference *h_ref = (Reference*) calloc(1,sizeof(temp_ref)); // change
		double *h_res = (double*) calloc(numb2*num_ipvec,sizeof(double)); // N,X data forward
		int *h_minres = (int*) calloc(num_ipvec*(numb2),sizeof(int));
		
		h_ref = ref ;
		
		
		printf("h_ref test 0 4  main %d \n ", sizeof(h_ref));
		printf("h_ref test 0 4  main %d \n ", sz_ref);
		//printf("h_ref test 1 2 main %d \n ", h_ref[1].input_vec[2]);
		//printf("h_ref test 0 4  main %d \n ", ref[0].input_vec[4]);
		//printf("h_ref test 1 2 main %d \n ", ref[1].input_vec[2]);
		
		
		unsigned int count[id]; // store number of inputs taking particular cluster in a level
 		vector<unsigned int> inputs_atnode[id];
		
		  for(int i=0 ; i< id ; i++) // initialize
		{
			count[i] =0;
		}	
		
		h_inputvec = inputs ;
		
		for(int in =0 ; in<numb2*num_ipvec ; in++) // initialze result to maximum value (padding) ; essential to find the minimum MMSE
		{
				h_res[in] = numeric_limits<int>::max();
				h_minres[in] = in % numb2 ;
		}
							
  		cl_mem c_inputvec = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sz_ipv, h_inputvec, &err);
		cl_mem c_ref = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,sz_ref, h_ref, &err);
		cl_mem c_res = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,sz_res, h_res, &err);
		cl_mem c_minres =  clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,sz_minres, h_minres, &err);
 								
		printf(" buffer created \n");
		clFinish(queue); 
		

		// set arguments 	
  							
		clSetKernelArg(kernel, 0, sizeof(cl_mem), &c_inputvec);
		clSetKernelArg(kernel, 1, sizeof(cl_mem), &c_ref);
		clSetKernelArg(kernel, 2, sizeof(cl_mem), &c_res);
		clSetKernelArg(kernel, 3, sizeof(const int), &num_ipvec);
		clSetKernelArg(kernel, 4, sizeof(const int), &N);
		clSetKernelArg(kernel, 5, sizeof(const int), &numb);
		clSetKernelArg(kernel, 6, sizeof(const int), &numb2);
		printf(" arguments to kernel \n");
       	
		local[0] = {Nt};
		global[0] = {N};
	
		//launch kernel			
       	clEnqueueNDRangeKernel(queue, kernel, dim, 0, global, local, 0, (cl_event*)NULL, &event_time);  
		
		clEnqueueReadBuffer(queue, c_res, CL_TRUE, 0,sz_res, h_res, 0, 0, 0); // how not necessary how will we refernece it?
		printf("reading ans \n");
      	clFinish(queue);
		
		//free(h_inputvec);
		
		for(int i=0 ; i<numb2*num_ipvec ; i++)
		{
			printf("%d val %f \n",i,h_res[i]);
		}
		
		printf("min cal \n");
		int offset=numb2/2;
			
		// set arguments for kernel2	
		clSetKernelArg(kernel1, 0,sizeof(cl_mem), &c_res);
		clSetKernelArg(kernel1, 1, sizeof(const int), &num_ipvec);
		clSetKernelArg(kernel1, 2,sizeof(cl_mem), &c_minres);
		clSetKernelArg(kernel1, 3, sizeof(const int), &numb2);
			
		// launch kernel 
		for(offset = numb2/2 ; offset>0 ; offset >>= 1)
		{
			printf("offset %d \n",offset);
			
			clSetKernelArg(kernel1, 4, sizeof(const int), &offset);
			
			global[0] = offset ;
			local[0] = offset ;
			
			clEnqueueNDRangeKernel(queue, kernel1, dim, 0 , global, local, 0, (cl_event*)NULL, &event_time);
			
			//clEnqueueReadBuffer(queue, c_minres, CL_TRUE, 0,sz_minres, h_minres, 0, 0, 0);
			clFinish(queue);
		}
								
		clEnqueueReadBuffer(queue, c_minres, CL_TRUE, 0,sz_minres, h_minres, 0, 0, 0);
		clEnqueueReadBuffer(queue, c_res, CL_TRUE, 0,sz_res, h_res, 0, 0, 0);
		
	
		
		//clFinish(queue);		
		// minimum calculation
		// count initialize to zero , even the othr two
		
		// based on the MMSE determin the next level node to be taken
		correct  =0 ;
		error2 = 0;
		for(int i=0; i<=(num_ipvec-1) ; i=i+1)
		{
			node_taken[level][i] = numb*node_taken[level-1][i] + h_minres[i*numb2];
			count[node_taken[level][i]] += 1;
			inputs_atnode[node_taken[level][i]].push_back(i);
			//printf( "node_taken by input %d = %d \n ",i, node_taken[level][i]);
			
			error2 = error2 + h_res[numb2*i];
			
		if (h_res[numb2*i] == double(0.000))
		{
			correct = correct + 1 ;
		}
		}
		
		error2 = error2 / num_ipvec ;
		
												
		int sum = 0;
		
		free(h_res);
		free(h_minres);
		
	
		
		if(level != num_level)
		{		
			for(int i=0 ; i<id; i++)
			{
				printf("i = %d \n",count[0]);
				for(int j=0 ; j<count[i]; j++)
				{
			
					free(ref[sum+j].input_vec);
					ref[sum+j].input_vec = (cl_uint*)malloc(count[i]*sizeof(unsigned int));
					//cl_uint *check = (unsigned int*) realloc(ref[sum+j].input_vec , count[i]);
					//if ( check != NULL )
					//{ 
						//ref[sum+j].input_vec = check ;
					//}
					//else
					//{ 
						//printf("error");
					//}
						
															
					ref[sum+j].noderef = level_node[level][i];
					memcpy (ref[sum+j].input_vec,&inputs_atnode[i][0],count[i]*sizeof(unsigned int));
					ref[sum+j].num = count[i];																				
				}
				
				
			sum += count[i];
			inputs_atnode[i].clear();
					//printf("h_ref test 0 4  main %d \n ", ref[0].input_vec[4]);
					//printf("h_ref test 1 2 main %d \n ", ref[1].input_vec[2]);
																					
			}
		}
		clFinish(queue);
		
		
	}
	double endTime = MPI_Wtime();
	diff = endTime - startTime ;
	avgTimePerIteration += diff;
}
	avgTimePerIteration = avgTimePerIteration / num_iter;
	
								
	for(int i=0 ; i <num_ipvec ; i++)
	{
	//	printf("node_taken = %d \n ",node_taken[num_level][i]);
	
	
	}
	
	double error1 = double((num_ipvec - correct))/double(num_ipvec) ;
	printf(" correct mapped : %d \n"  , correct );
	printf(" error : %f \n"  , error1 );
	printf("avg error : %f \n" , error2);

	
	printf("Average time per iteration : %3.5e s\n"  , avgTimePerIteration );
	printf("end loop \n");
	/*	
								int g = (int)(node_taken[num_level][0]/2) ;
								printf("closest search results \n");
								
								printf(" dict closest to input 0 \n");
								
								//int g = (int)(node_taken[num_level][i]/2);
								printf(" dict closest to input 0 %f \n",level_node[num_level-1][g].centroids[(node_taken[num_level][0]%2)]);
								
								 g = (int)(node_taken[num_level][1]/2) ;
								//printf(" dict closest to input 0 %f \n",level_node[2][g].centroids[0]);
								
								
								printf(" %d dict closest to input 1 %f \n",g,level_node[num_level-1][g].centroids[(node_taken[num_level][1]%2)]);
								printf(" %d dict closest to input 1 %f \n",g,level_node[num_level-1][0].centroids[441]);
								*/		


 }
 
 


