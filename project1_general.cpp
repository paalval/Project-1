// Numerically solves the equation -u''(x) = f(x) for u(x) with f(x) known for x = [0, 1].
// The equation is rewritten as a set of linear equations, Av = w, which will be solved using Gaussian elimination.
// This code considers the problem when A is a tridiagonal matrix where the lower-, main- and upper diagonal entries are elements in the vectors c, a and b respectively.

#include "stdafx.h"
#include <iostream>
#include <armadillo>
#include "time.h"

using namespace std;
using namespace arma;

int main()
{
	// Define number of grid points n and steplength h
	float n = 100000;
	float h = 1 / (1 + n);

	// Start taking the time
	clock_t start, finish;
	start = clock();

	// Define data points
	vec x(n);
	for (int i = 0; i < n; i = i + 1)
	{
		x(i) = (i + 1)*h;
	}

	// Define function f(x) and vector w
	vec f = 100 * exp(-10 * x);
	vec w = h*h*f;

	// Define lower-, main- and upper diagonal vectors c, a, and b
	vec c(n); c.fill(-1);
	vec a(n); a.fill(2);
	vec b(n); b.fill(-1);

	// New vector elements from forward substitution (GENERAL CASE)
	vec af(n); af(0) = a(0);
	vec wf(n); wf(0) = w(0);
	for (int i = 1; i < n; i = i + 1)
	{
		af(i) = a(i) - b(i - 1)*(c(i - 1) / af(i - 1));
		wf(i) = w(i) - wf(i - 1)*(c(i - 1) / af(i - 1));
	}

	// Solution vector elements from backward substitution (GENERAL CASE)
	vec v(n); v(n - 1) = wf(n - 1) / af(n - 1);
	int i = n - 1;
	while (i > 0)
	{
		v(i - 1) = (wf(i - 1) - b(i - 1)*v(i)) / af(i - 1);
		i = i - 1;
	}

	// Define analytical solution and save
	vec u = 1 - (1 - exp(-10))*x - exp(-10 * x);

	// Stop taking the time
	finish = clock();
	cout << "Takes a total of " << ((finish - start) / (double CLOCKS_PER_SEC)) << " seconds." << endl;

	// Save calculated vectors to later be graphically displayed using different software
	x.save("data_points.txt", raw_ascii);
	v.save("numerical.txt", raw_ascii);
	u.save("analytical.txt", raw_ascii);
	return 0;
}

