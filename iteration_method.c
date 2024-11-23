#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>

#include "additional.h"

int main()
{
	double **ab = NULL; //the augmented matrix of the given system
	int e_out; //number of digits after the decimal point
	int n; //number of variables
	int i, j, number_iter = 0; //counters
	int r_a, r_ab; //ranks of the coefficient and the augmented matrices
	double **c = NULL, *d = NULL; //the coefficient matrix and the column vector of constants of the equivalent system
    double *x0 = NULL, *x1 = NULL; //the variable matrices
    double e = 1, delta; //accuracy and calculation error
    double norm_c; //m-norm of the matrix c
	//read the augmented matrix from the input file
	ab = input("input.txt", &n);
    printf("ITERATION METHOD for solving system of linear equations\n");
    for (i = 0; i < n; i++)
	{
		for (j = 0; j < n - 1; j++)
		{
			printf("%.2lf*x%d + ", ab[i][j], j);
		}
        printf("%.2lf*x%d = %.2lf\n", ab[i][j], j, ab[i][j + 1]);
	}
	//check if the given system is consistent and has a unique solution
	if (!check(ab, n))
	{
		exit(3);
	}
    do
    {
        printf("enter desired accuracy - number of digits after the decimal point -> ");
        scanf("%d", &e_out);
    } while (e_out <= 0);
	//allocate memory for the arrays
	c = (double**)malloc(n * sizeof(double*));
	d = (double*)malloc(n * sizeof(double));
	x0 = (double*)malloc(n * sizeof(double));
	x1 = (double*)malloc(n * sizeof(double));
	if (!c || !d || !x0 || !x1)
	{
		exit(2);
	}
	for (i = 0; i < n; i++)
	{
		c[i] = (double*)malloc(n * sizeof(double));
		if (!c[i])
		{
			exit(2);
		}
	}
	//calculate matrices for the equivalent system
	equivalent(ab, n, c, d);
	//calculate m-norm of the matrix c
	norm_c = norm(c, n, n);
	//check the convergence of the iteration method
	if (norm_c >= 1)
	{
		printf("iteration method can't be used\n");
		return 0;
	}
	//calculate initial accuracy
	for (i = 0; i < e_out; i++)
	{
		e /= 10;
	}
    //calculate accuracy for comparing with delta during the iteration process
	e = ((1 - norm_c) / norm_c) * e;
	//define the initial approximation of x
	for (i = 0; i < n; i++)
	{
		x1[i] = d[i];
	}
	//perform iterations to achieve a given accuracy
	do
	{
		for (i = 0; i < n; i++)
		{
			x0[i] = x1[i];
		}
		//x1=c*x0+d
		for (i = 0; i < n; i++)
		{
			x1[i] = 0;
			for (j = 0; j < n; j++)
			{
				x1[i] += c[i][j] * x0[j];
			}
			x1[i] += d[i];
		}
		delta = 0;
		//calculate delta
		for (i = 1; i < n; i++)
		{
			if (delta < fabs(x1[i] - x0[i]))
			{
				delta = fabs(x1[i] - x0[i]);
			}
		}
		number_iter++;
	} while (delta >= e);
	//print results
	printf("X = (");
	for (i = 0; i < n - 1; i++)
	{
		printf("%.*lf; ", e_out, x1[i]);
	}
	printf("%.*lf)\ndelta = %lg (e = %lg)\nnumber of iterations: %d\n", e_out, x1[i], delta, e, number_iter);
	//free the memory
	for (i = 0; i < n; i++)
	{
		free(c[i]);
		free(ab[i]);
	}
	free(ab);
	free(c);
	free(d);
	free(x0);
	free(x1);
	return 0;
}
