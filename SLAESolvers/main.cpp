#include <iostream>
#include "DirectSolvers.h"
#include <fstream>

int main()
{
	try
	{
		std::ifstream file("A.txt");

		int n, m;

		file >> n >> m;

		Matrix<double> A(n, m);

		for (int i = 0; i < n; ++i)
		{
			for (int j = 0; j < m; ++j)
			{
				file >> A[i][j];
			}
		}

		file.close();

		/*file.open("F.txt");

		file >> n;

		Vector<double> F(n);

		for (int i = 0; i < n; ++i)
		{
			file >> F[i];
		}

		file.close();*/

		Matrix<double> Q(n, m);
		Matrix<double> R(n, m);

		for (int i = 0; i < Q.getRowCount(); ++i)
			Q[i][i] = 1;

		GivensTransformation::GivensOrthogonalization(A, Q, R);
		Q.ConsolePrint();
		std::cout << '\n';
		R.ConsolePrint();

		/*Vector<double> res = GaussMethod::Solve(A, F);
		res.ConsolePrint();

		std::cout << '\n';

		LUDecomposition<double> LUSolver(A);
		res = LUSolver.Solve(F);
		res.ConsolePrint();\*/

		int k;
		std::cout << "enter1 : ";
		std::cin >> k;
		std::cout << '\n';
	}
	catch (std::exception& ex)
	{
		std::cout << ex.what() << '\n';
	}

	int k;
	std::cout << "enter2 : ";
	std::cin >> k;

	return 0;
}