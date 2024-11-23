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
	double **c = NULL, **c1 = NULL, *d = NULL; //the coefficient matrix, the coefficient matrix with zeroes below the main diagonal and the column vector of constants of the equivalent system
    double *x0 = NULL, *x1 = NULL; //the variable matrices
    double e = 1, delta; //accuracy and calculation error
    double norm_c, norm_c1; //m-norm of the matrices c and c1
	//read the augmented matrix from the input file
	ab = input("input.txt", &n);
    printf("SEIDEL'S METHOD for solving system of linear equations\n");
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
    c1 = (double**)malloc(n * sizeof(double*));
	d = (double*)malloc(n * sizeof(double));
	x0 = (double*)malloc(n * sizeof(double));
	x1 = (double*)malloc(n * sizeof(double));
	if (!c || !c1 || !d || !x0 || !x1)
	{
		exit(2);
	}
	for (i = 0; i < n; i++)
	{
		c[i] = (double*)malloc(n * sizeof(double));
		c1[i] = (double*)malloc(n * sizeof(double));
		if (!c[i] || !c1[i])
		{
			exit(2);
		}
	}
	//calculate matrices for the equivalent system
	equivalent(ab, n, c, d);
	//calculate m-norm of the matrix c
	norm_c = norm(c, n, n);
	//check the convergence of the seidel's method
	if (norm_c >= 1)
	{
		printf("seidel's method can't be used\n");
		return 0;
	}
	//calculate matrix с1
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			if (j < i)
			{
				c1[i][j] = 0;
			}
			else
			{
				c1[i][j] = c[i][j];
			}
		}
	}
	//calculate m-norm of the matrix с1
	norm_c1 = norm(c1, n, n);
	//calculate initial accuracy
	for (i = 0; i < e_out; i++)
	{
		e /= 10;
	}
    //calculate accuracy for comparing with delta during the iteration process
	e = ((1 - norm_c) / norm_c1) * e;
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
		//calculate х1
		for (i = 0; i < n; i++)
		{
			x1[i] = 0;
			for (j = 0; j <= i - 1; j++)
			{
				x1[i] += c[i][j] * x1[j];
			}
			for (j = i + 1; j < n; j++)
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
		printf("%.*f; ", e_out, x1[i]);
	}
	printf("%.*lf)\ndelta = %lg (e = %lg)\nnumber of iterations: %d\n", e_out, x1[i], delta, e, number_iter);
	//free the memory
	for (i = 0; i < n; i++)
	{
		free(c[i]);
	}
	free(c);
	free(d);
	free(x0);
	free(x1);
	return 0;
}
