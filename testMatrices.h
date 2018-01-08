/*
 * testMatrices.h
 *
 *  Created on: Jan 8, 2018
 *      Author: maxim koroteev
 */

#ifndef TESTMATRICES_H_
#define TESTMATRICES_H_

// Two functions below implement two test matrices used in WR (Wilkinson-Reinsch)
// Handbook to test their implementation of Jacobi rotations. Use them to reproduce
// results indicated in WR to convince yourself the algorithm is correct. Note,
// the "correctness" does not propagate further than WR tests, so it may be the algorithm
// would fail on some convoluted matrices.

// this matrix was used in WR, page 208. N=30
void testMatrix_1(double**& mymatr, int N)
{
	for (int i = 1; i<N + 1; i++)
		for (int j = 1; j<N + 1; j++)
			mymatr[i - 1][j - 1] = (i>j) ? i : j;

	// print the test matrix
	for (int i = 0; i<N; i++)
	{
		for (int j = 0; j<N; j++)
			std::cout << mymatr[i][j] << " ";
		std::cout << std::endl;
	}
	// -->

   return;
}

// this matrix was used in WR, page 209, N=44
void testMatrix_2(double**& mymatr, int N)
{
    // initialization
	for (int i = 0; i<N; i++)
		for (int j = 0; j<N; j++)
			mymatr[i][j] = 0.0;

	// main diagonal and two sub-diagonals
	for (int i = 1; i < N; i++)
	{
		mymatr[i][i]     = 6.0;
		mymatr[i][i + 1] = 3.0;
		mymatr[i][i-1]   = 3.0;
	}
	mymatr[0][0] = 5.0;     mymatr[N-1][N-1] = 5.0;
	mymatr[0][1] = 2.0;     mymatr[1][0] = 2.0;
	mymatr[N-1][N - 2] = 2; mymatr[N - 2][N-1] = 2.0;

	// further sub-diagonals
	for (int i = 0; i < N - 2; i++)
	{
		mymatr[i][i + 2] = 1.0;
		mymatr[i + 2][i] = 1.0;
	}
	for (int i = 0; i < N - 3; i++)
	{
		mymatr[i][i + 3] = 1.0;
		mymatr[i + 3][i] = 1.0;
	}

	// print the test matrix
	for (int i = 0; i<N; i++)
	{
		for (int j = 0; j<N; j++)
			std::cout << mymatr[i][j] << " ";
		std::cout << std::endl;
	}
	// -->

   return;
}

#endif /* TESTMATRICES_H_ */
