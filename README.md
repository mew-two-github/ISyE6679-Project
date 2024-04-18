### Short description of important files:
- driver.cpp: the main function which solves the optimization problem
- solverfile.cu: computes solution for a given sparse matrix linear system
- read_data.h: loads data into global variables accessible across functions
- utils_derivs.h: contains the main utilities, mathematical functions (including Lagrangian) and their first derivatives
- second_derivs.h: setup to compute the second derivative of the Lagrangian

### Requirements
- CoDiPack
- fast-csv-cpp-parser
- Boost Library for BLAS
#### Setting up the project
```Shell
mkdir project
cd project
mkdir CoDiPack
mkdir fast-csv-cpp-parser
mkdir C_code
cd C_code
git clone https://github.com/mew-two-github/ISyE6679-Project.git
cd ../fast_csv_parser
git clone https://github.com/ben-strasser/fast-cpp-csv-parser
cd ../CoDiPack
run git clone https://github.com/SciCompKL/CoDiPack.git
```
Download boost from here: https://www.boost.org/users/download/ and put it in the project folder
```Shell
tar -xvf "../boost_1_84_0"
```
And you are all set!
### Running the code
Commands to run the code locally (comment out line 185 in driver.cpp):
```Shell
g++ driver.cpp -I "./" -I "../CoDiPack/include" -I "../fast-csv-cpp-parser" -I "../boost_1_84_0"
```

Commands to run it on PACE:
```Shell
g++ -o main.o -c driver.cpp -I "./" -I "../CoDiPack/include" -I "../fast-csv-cpp-parser" -I "../boost_1_84_0"
nvcc -c solverfile.cu -o solverfile.o -lcusolver -lcusparse
g++ main.o solverfile.o -o executable -L/usr/local/cuda/lib64 -lcudart -lcusolver -lcusparse -pthread
```
Submit the executable to the queue.
### Other comments

Currently, PACE gives a seg fault in the matrix prep stage. We identified the function where the error occurs, but we are not able to zero in on which line is causing the error. More strangely, the code runs fine locally, 
but unfortunately, we can't use CuSolver locally.

All the modules have been extensively tested separately, including the CUDA code. You might be able to find codes that help you do that in either ./cusolver_trials/ or in ./old_code/
