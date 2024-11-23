# C programs for solving a system of linear equations
 File "input.txt" contains the number of variables and the system represented in the form of an augmented matrix
### Algorithm:
1. Read data from the input.txt
2. Check if the given system is consistent and has a unique solution (using rank)
3. Calculate matrices for the equivalent system
4. Check the convergence of the method (using m-norm)
5. Define the initial approximation of X
6. Perform iterations to achieve a given accuracy (using iteration method or Seidel's method)
7. Print results
