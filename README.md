# mobile-computing

Very basic version of the code is shared, for advanced version which uses alpha parameter, shared memory, Snapdragon version etc
contact p[dot]anushri25[at]gmail[dot]com

Steps in the code

1) Extract patches from images 
2) Build dictionaries from patches
3) Build a balanced hierarchical shallow tree using the dictionaries (matlab-files)
4) Use the Shallow tree for searching on CPU-GPU platform (snapdragon)

inorder to establish mex integration between matlab codes and OpenCL codes linking need to be performed

make using: (use corresponding folder)

 mpicxx -L/opt/intel/opencl-1.2-4.4.0.117/lib64 -L/usr/local/MATLAB/MATLAB_Compiler_Runtime/v83/bin/glnxa64 main.cpp -o main -lOpenCL -I/usr/include/CL -lmat -lmx -lmex -I/usr/local/MATLAB/MATLAB_Compiler_Runtime/v83/extern/include -Wl,-rpath,/opt/intel/opencl-1.2-4.4.0.117/lib64,-rpath,/usr/local/MATLAB/MATLAB_Compiler_Runtime/v83/bin/glnxa64 -fpermissive -g

execute using:

 LD_PRELOAD=/opt/intel/opencl-1.2-4.4.0.117/lib64/libtbb.so.2  ./main
