#include <iostream>
#include "DirectSolvers.h"
#include <fstream>
#include <chrono>

template<class T>
double Relative_Error(Vector<T> X, Vector<T> x)
{
	double s = 0.0;

	for (int i = 0; i < X.size(); i++)
		s += std::pow(X[i] - x[i], 2);

	return std::sqrt(s) / x.Norma();
}

int main()
{
	try
	{
		int n = 500;
		Matrix<double> A(n, n);
		Vector<double> X_true(n);

		for (int i = 0; i < n; ++i)
		{
			for (int j = 0; j < n; ++j)
			{
				A[i][j] = 1.0 / (i + j + 1.0);
			}
			X_true[i] = 1;
		}

		auto F = A * X_true;

		auto start = std::chrono::steady_clock::now();
		
		QRDecomposition<double> solver(A, QRAlgorithm::Householder);
		auto X = solver.Solve(F);

		auto end = std::chrono::steady_clock::now();
		auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
		
		X.ConsolePrint();

		std::cout << '\n' << "Error: " << Relative_Error(X, X_true) << '\n';
		std::cout << "Runtime: " << elapsed.count() << '\n';

		int x;
		std::cin >> x;
	}
	catch (std::exception& ex)
	{
		std::cout << ex.what() << '\n';
	}

	return 0;
}