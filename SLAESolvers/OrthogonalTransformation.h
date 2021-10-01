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

template<class T>
void HouseholderTransformation::HouseholderOrthogonalization(const Matrix<T>& A, Matrix<T>& Q, Matrix<T>& R)
{
	R = A;
	double tmp, beta, mu;
	Vector<T> p(R.getRowCount());

	for (int i = 0; i < R.getColumnCount() - 1; ++i)
	{
		tmp = 0.0;

		for (int k = i; k < R.getRowCount(); ++k)
			tmp += std::pow(R[k][i], 2);

		if (std::sqrt(std::abs(tmp - R[i][i] * R[i][i])) > CONST::EPS)
		{
			beta = std::sqrt(tmp);

			if (R[i][i] >= 0)
				beta = -beta;

			mu = 1.0 / beta / (beta - R[i][i]);

			for (int k = 0; k < R.getRowCount(); ++k)
			{
				p[k] = 0;
				
				if(k >= i)
					p[k] = R[k][i];
			}

			p[i] -= beta;

			for (int m = i; m < R.getColumnCount(); ++m)
			{
				tmp = 0;

				for (int n = i; n < R.getRowCount(); ++n)
					tmp += R[n][m] * p[n];

				tmp *= mu;

				for (int n = i; n < R.getRowCount(); ++n)
					R[n][m] -= tmp * p[n];
			}

			for (int m = 0; m < Q.getRowCount(); ++m)
			{
				tmp = 0;

				for (int n = i; n < Q.getRowCount(); ++n)
					tmp += Q[m][n] * p[n];

				tmp *= mu;

				for (int n = i; n < Q.getRowCount(); ++n)
					Q[m][n] -= tmp * p[n];
			}
		}
	}
}
