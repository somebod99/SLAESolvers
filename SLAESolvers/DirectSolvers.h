#pragma once
#include "Substitution.h"
#include "OrthogonalTransformation.h"

class GaussMethod
{
public:
	template<class T>
	static int FindMainElement(const Matrix<T>& A, const int& j);

	template<class T>
	static void DirectWay(Matrix<T>& A, Vector<T>& F);

	template<class T>
	static void DirectWay(Matrix<T>& A);

	template<class T>
	static Vector<T> Solve(Matrix<T> A, Vector<T> F);
};

template<class T>
int GaussMethod::FindMainElement(const Matrix<T>& A, const int& j)
{
	int index = j;
	for (int i = j + 1; i < A.getRowCount(); ++i)
		if (std::abs(A[i][j]) > std::abs(A[index][j]))
			index = i;

	if (std::abs(A[index][j]) < CONST::EPS)
		throw std::exception("Degenerate matrix");

	return index;
}

template<class T>
void GaussMethod::DirectWay(Matrix<T>& A, Vector<T>& F)
{
	for (int i = 0; i < A.getRowCount() - 1; ++i)
	{
		int ind = FindMainElement(A, i);

		if (ind != i)
		{
			auto temp = A[i];
			A[i] = A[ind];
			A[ind] = temp;

			auto temp_f = F[i];
			F[i] = F[ind];
			F[ind] = temp_f;
		}

		for (int j = i + 1; j < A.getRowCount(); ++j)
		{
			auto del = A[j][i] / A[i][i];

			A[j][i] = 0;

			for (int k = i + 1; k < A.getColumnCount(); ++k)
			{
				A[j][k] -= del * A[i][k];
			}

			F[j] -= del * F[i];
		}
	}
}

template<class T>
void GaussMethod::DirectWay(Matrix<T>& A)
{
	for (int i = 0; i < A.getRowCount() - 1; ++i)
	{
		int ind = FindMainElement(A, i);

		if (ind != i)
		{
			auto temp = A[i];
			A[i] = A[ind];
			A[ind] = temp;
		}

		for (int j = i + 1; j < A.getRowCount(); ++j)
		{
			auto del = A[j][i] / A[i][i];

			A[j][i] = 0;

			for (int k = i + 1; k < A.getColumnCount(); ++k)
			{
				A[j][k] -= del * A[i][k];
			}
		}
	}
}

template<class T>
Vector<T> GaussMethod::Solve(Matrix<T> A, Vector<T> F)
{
	DirectWay(A, F);

	return Substitution::BackRowSubstitution(A, F);
}

template<class T>
class LUDecomposition
{
private:
	Matrix<T> LU;

public:
	LUDecomposition(const Matrix<T>&);
	Vector<T> DirectWay(const Matrix<T>&, const Vector<T>&);
	Vector<T> Solve(const Vector<T>&);
};

template<class T>
LUDecomposition<T>::LUDecomposition(const Matrix<T>& A)
{
	LU = A;

	GaussMethod::DirectWay(LU);

	for (int i = 0; i < LU.getRowCount(); ++i)
	{
		for (int j = 0; j < i; ++j)
		{
			double sum = 0;

			for (int k = 0; k < j; ++k)
			{
				sum += LU[i][k] * LU[k][j];
			}

			LU[i][j] = (A[i][j] - sum) / LU[j][j];
		}
	}
}

template<class T>
Vector<T> LUDecomposition<T>::DirectWay(const Matrix<T>& A, const Vector<T>& F)
{
	Vector<T> res(F);

	for (int i = 1; i < F.size(); ++i)
	{
		for (int j = 0; j < i; ++j)
		{
			res[i] -= A[i][j] * res[j];
		}
	}

	return res;
}

template<class T>
Vector<T> LUDecomposition<T>::Solve(const Vector<T>& F)
{
	auto res = DirectWay(LU, F);
	res.ConsolePrint();
	return Substitution::BackRowSubstitution(LU, res);
}

enum class QRAlgorithm
{
	Householder,
	Givens
};

template<class T>
class QRDecomposition
{
private:
	Matrix<T> Q;
	Matrix<T> R;

public:
	QRDecomposition(const Matrix<T>& A, const QRAlgorithm& method);
	Vector<T> Solve(const Vector<T>& F);
};

template<class T>
QRDecomposition<T>::QRDecomposition(const Matrix<T>& A, const QRAlgorithm& method)
{
	Q.resize(A.getRowCount(), A.getColumnCount());

	switch (method)
	{
	case QRAlgorithm::Givens:
		for (int i = 0; i < A.getRowCount(); ++i)
			Q[i][i] = 1;

		GivensTransformation::GivensOrthogonalization(A, Q, R);
		break;
		
	case QRAlgorithm::Householder:
		for (int i = 0; i < A.getRowCount(); ++i)
			Q[i][i] = 1;

		HouseholderTransformation::HouseholderOrthogonalization(A, Q, R);
		break;
	}
}

template<class T>
Vector<T> QRDecomposition<T>::Solve(const Vector<T>& F)
{
	auto res = Q.MultiplicationTransMatrixVector(F);
	return Substitution::BackRowSubstitution(R, res);
}
