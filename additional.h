double** input(const char *filename, int *n); //function to read data from input file
int rank(double **matrix, int n_row, int n_col); //function to calculate rank of the matrix
double norm(double **matrix, int n_row, int n_col); //function to calculate m-norm of the matrix
void swap(double **matrix, int n_col, int k1, int k2); //function to swap two rows of the matrix
int check(double **ab, int n); //function to check if the given system is consistent and has a unique solution
void equivalent(double **ab, int n, double **c, double *d); //function to calculate matrices of the equivalent system

double** input(const char *filename, int *n)
{
	int i, j;
	double **ab = NULL;
	FILE *p = NULL;
	p = fopen(filename, "r"); //open file
	if (!p)
	{
		printf("file can't be opened\n");
		exit(1);
	}
	fscanf(p, "%d", n); //read the number of variables
	//allocate memory for the matrix ab
	ab = (double**)malloc((*n) * sizeof(double*));
	for (i = 0; i < *n; i++)
	{
		ab[i] = (double*)malloc(((*n) + 1) * sizeof(double));
		if (!ab[i])
		{
			exit(2);
		}
		//read a row from the file to ab
		for (j = 0; j < *n + 1; j++)
		{
			fscanf(p, "%lf", &ab[i][j]);
		}
	}
	return ab;
}

int rank(double **matrix0, int n_row, int n_col)
{
	int rank = 0, k, i, j;
	//create a copy of the given matrix
	double **matrix = NULL;
	matrix = (double**)malloc(n_row * sizeof(double*));
	if (!matrix)
	{
		exit(2);
	}
	for (i = 0; i < n_row; i++)
	{
		matrix[i] = (double*)malloc((n_col) * sizeof(double));
		if (!matrix[i])
		{
			exit(2);
		}
		for (j = 0; j < n_col; j++)
		{
			matrix[i][j] = matrix0[i][j];
		}
	}
	//transform matrix into row echelon form using Gaussian method
	for (k = 0; k < n_row - 1; k++)
	{
		if (matrix[k][k] == 0)
		{
			for (i = k + 1; i < n_row; i++)
			{
				if (matrix[i][k] != 0)
				{
					swap(matrix, n_col, k, i);
					break;
				}
			}
			if (matrix[k][k] == 0)
			{
				k++;
			}
		}
		for (i = k + 1; i < n_row; i++)
		{
			for (j = k; j < n_col; j++)
			{
				if (j == k)
				{
					matrix[i][j] = 0;
				}
				else
				{
					matrix[i][j] = matrix[k][k] * matrix[i][j] - matrix[i][k] * matrix[k][j];
				}
			}
		}
	}
	//calculate the number of non-zero rows
	for (i = 0; i < n_row; i++)
	{
		for (j = i; j < n_col; j++)
		{
			if (matrix[i][j])
			{
				rank++;
				break;
			}
		}
	}
    //free the memory
	for (i = 0; i < n_row; i++)
	{
		free(matrix[i]);
	}
	free(matrix);
	return rank;
}

double norm(double **matrix, int n_row, int n_col)
{
	int i, j;
	double sumrow, norm = 0;
	for (i = 0; i < n_row; i++)
	{
		sumrow = 0;
		//calculate the sum of the absolute values of the row elements
		for (j = 0; j < n_col; j++)
		{
			sumrow += fabs(matrix[i][j]);
		}
		//find the maximum sum
		if (norm < sumrow)
		{
			norm = sumrow;
		}
	}
	return norm;
}

void swap(double **matrix, int n_col, int k1, int k2)
{
	int j;
	double temp;
	for (j = 0; j < n_col; j++)
	{
		temp = matrix[k1][j];
		matrix[k1][j] = matrix[k2][j];
		matrix[k2][j] = temp;
	}
}

int check(double **ab, int n)
{
	int r_a, r_ab;
	//calculate ranks of the coefficient and the augmented matrices
	r_a = rank(ab, n, n);
	r_ab = rank(ab, n, n + 1);
	if (r_a != r_ab)
	{
		printf("system has no solution\n");
		return 0;
	}
	if (r_a < n)
	{
		printf("system has infinitely many solutions\n");
		return 0;
	}
    printf("system is consistent and has a unique solution\n");
	return 1;
}

void equivalent(double **ab, int n, double **c, double *d)
{
	int i, j;
	for (i = 0; i < n; i++)
	{
        if (ab[i][i] == 0)
        {
            printf("division by zero\n");
            exit(3);
        }
		d[i] = ab[i][n] / ab[i][i];
		for (j = 0; j < n; j++)
		{
			if (i == j)
			{
				c[i][j] = 0;
			}
			else
			{
				c[i][j] = -ab[i][j] / ab[i][i];
			}
		}
	}
}
