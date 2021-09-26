#pragma once
#include "Matrix.h"
#include "CONST.h"

class GivensTransformation
{
public:
	template<class T>
	static void GivensOrthogonalization(const Matrix<T>&, Matrix<T>&, Matrix<T>&);
};

template<class T>
void GivensTransformation::GivensOrthogonalization(const Matrix<T>& A, Matrix<T>& Q, Matrix<T>& R)
{
	double temp1, temp2;
	R = A;

	for (int j = 0; j < R.getColumnCount() - 1; ++j)
	{
		for (int i = j + 1; i < R.getRowCount(); ++i)
		{
			if (std::abs(R[i][j]) > CONST::EPS)
			{
				double del = std::sqrt(std::pow(R[i][j], 2) + std::pow(R[j][j], 2));
				double c = R[j][j] / del;
				double s = R[i][j] / del;

				for (int k = j; k < R.getColumnCount(); ++k)
				{
					temp1 = c * R[j][k] + s * R[i][k];
					temp2 = c * R[i][k] - s * R[j][k];
					R[j][k] = temp1;
					R[i][k] = temp2;
				}

				for (int k = 0; k < Q.getRowCount(); ++k)
				{
					temp1 = c * Q[k][j] + s * Q[k][i];
					temp2 = c * Q[k][i] - s * Q[k][j];
					Q[k][j] = temp1;
					Q[k][i] = temp2;
				}
			}
		}
	}
}

class HouseholderTransformation
{
public:
	template<class T>
	static void HouseholderOrthogonalization(const Matrix<T>&, Matrix<T>&, Matrix<T>&);
};