// Numerically solves the equation -u''(x) = f(x) for u(x) with f(x) known for x = [0, 1].
// The equation is rewritten as a set of linear equations, Av = w, which will be solved using Gaussian elimination.
// This code considers the problem when A is a tridiagonal matrix where the lower-, main- and upper diagonal entries are -1, 2 and -1 respectively.

#include "stdafx.h"
#include <iostream>
#include <armadillo>
#include "time.h"

using namespace std;
using namespace arma;

int main()
{	
	// Define number of grid points n and steplength h
	float n = 10;
	float h = 1 / (1 + n);

	// Start taking the time
	clock_t start, finish;
	start = clock();

	// Define data points
	vec x(n);
	for(int i = 0; i < n; i = i + 1)
	{
		x(i) = (i + 1)*h;
	}

	// Define function f(x) and vector w
	vec f = 100 * exp(-10 * x);
	vec w = h*h*f;

	// New vector elements from forward substitution (SPECIAL CASE)
	vec wf(n); wf(0) = w(0);
	for (float i = 1; i < n; i = i + 1)
	{
		wf(i) = w(i) + (i / (i + 1))*wf(i - 1);
	}
	// Solution vector elements from backward substitution (SPEACIAL CASE)
	vec v(n); v(n - 1) = wf(n - 1)*(n/(n + 1));
	float i = n - 1;
	while (i > 0)
	{
		v(i - 1) = (i / (i + 1))*(wf(i - 1) + v(i));
		i = i - 1;
	}

	// Define analytical solution and save
	vec u = 1 - (1 - exp(-10))*x - exp(-10 * x);

	// Stop taking the time
	finish = clock();
	cout << "Takes a total of " << ((finish - start) / (double CLOCKS_PER_SEC)) << " seconds." << endl;

	// Save calculated vecotrs to later be graphically displayed using different software
	x.save("data_points.txt", raw_ascii);
	v.save("numerical.txt", raw_ascii);
	u.save("analytical.txt", raw_ascii);
	return 0;
}

