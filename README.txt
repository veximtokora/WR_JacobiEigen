1. How to build:

> mkdir build
> cd build
> cmake ..
> make

2. the executable file of the code is supposed to be used in command line environment;

3. the command line has the form: ./WR_Jacobi_Eigen N eivec tr_param
where

N is the size of a symmetric square matrix;
eivec - a parameter equal to 0 or 1. 0 implies only eigen values computations, 1 includes both eigen vectors and eigen values;
tr-param - may have two values: 1) real; 2) test;
"test" uses one of the test matrices provided with the code; "real" implies a user's matrix would be used.

In the latter case a user has to create a matrix for which (s)he wants to compute eigen values and
eigen vectors manually. The matrix has to be written into the file matrix.dat which has to be put in the same folder
as the executable file. The matrix in the file matrix.dat has to be written column-wise in one long column whithout
empty lines. For example, if your matrix is 
     | 1 2 3 |
A =  | 2 1 2 |
     | 3 2 1 |

then this matrix has to be written into the file matrix.dat as follows

1
2
3
2
1
2
3
2
1

where first three numbers are the first column of A, second three numbers are the second column of A, etc. 

4. the output of the executable is written in one or two files depending on the parameter eivec. The names
of the files are
eigenValues.dat
eigenVectors.dat

5. Windows users have to make minor changes in the code before trying to compile it. It has to do with paths representation
(using back slashes, c:\ etc.)

