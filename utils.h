/*
 * utils.h
 *
 *  Created on: Jan 8, 2018
 *      Author: maxim koroteev
 */

#ifndef UTILS_H_
#define UTILS_H_

#include <algorithm>
#include <numeric>
#include <iterator>

template <typename T>
void sort_indexes(std::vector<T> &lam, std::vector<int>& idx)
{
	std::iota(idx.begin(), idx.end(), 0);

	// sort indexes based on comparing values in v
	std::sort(idx.begin(), idx.end(),
	     [&lam](size_t i1, size_t i2) { return lam[i1] < lam[i2]; } );
}

template <typename T>
void sort_indexes_backward(std::vector<T> &lam, std::vector<int>& idx)
{
	std::iota(idx.begin(), idx.end(), 0);

	// sort indexes based on comparing values in v
	std::sort(idx.begin(), idx.end(),
	     [&lam](size_t i1, size_t i2) { return lam[i1] > lam[i2]; } );
}

void MemAlloc(double**& z, int rows, int columns)
{
	if (!z)
	{
		z = new double*[rows];
		for (int i = 0; i< rows; i++)
			z[i] = new double[columns];
	}
}

void MemRelease(double**& z, int rows)
{
	for (int i = 0; i<rows; i++)
		delete[] z[i];
	delete[] z;
	z = NULL;
}

#endif /* UTILS_H_ */
