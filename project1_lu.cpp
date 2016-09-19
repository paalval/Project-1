// Numerically solves the equation -u''(x) = f(x) for u(x) with f(x) known for x = [0, 1].
// The equation is rewritten as a set of linear equations, Av = w, which will be solved by finding the inverse of A using LU-decomposition then doing the multiplication v = inv(A)w.
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
	float n = 6;
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
	w.print("w = ");

	// Define lower-, main- and upper diagonal vectors c, a, and b
	vec c(n - 1); c.fill(-1);
	vec a(n);     a.fill(2);
	vec b(n - 1); b.fill(-1);

	// Define matrix A
	mat A = zeros(n, n);
	A.diag()   = a; // Main diagonal
	A.diag(-1) = c; // Lower diagonal
	A.diag(1)  = b; // Upper diagonal
	A.print("A = ");

	// Use A.i to compute inverse of A, and calculte solution vector v
	//mat B = A.i();
	//vec v = B*w;

	// Use solve() to find the solution vector v
	//vec v = solve(A, w);

	// Decompose matrix A into lower- and upper triangular matrices L and U.
	mat L, U;
	lu(L, U, A);
	//L.print("L = ");
	//U.print("U = ");
	vec y = solve(U, x);
	vec v = solve(L, y);

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

