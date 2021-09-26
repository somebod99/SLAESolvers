#pragma once
#include "Matrix.h"
#include "CONST.h"

class Substitution
{
public:
	template<class T>
	static Vector<T> DirectRowSubstitution(const Matrix<T>& A, const Vector<T>& F);

	template<class T>
	static Vector<T> DirectColumnSubstitution(const Matrix<T>& A, const Vector<T>& F);

	template<class T>
	static Vector<T> BackRowSubstitution(const Matrix<T>& A, const Vector<T>& F);

	template<class T>
	static Vector<T> BackColumnSubstitution(const Matrix<T>& A, const Vector<T>& F);
};

template<class T>
Vector<T> Substitution::DirectRowSubstitution(const Matrix<T>& A, const Vector<T>& F)
{
	Vector<T> res(F);

	for (int i = 0; i < F.size(); ++i)
	{
		if (std::abs(A[i][i]) < CONST::EPS)
			throw std::exception("Division by 0");

		for (int j = 0; j < i; ++j)
		{
			res[i] -= A[i][j] * res[j];
		}

		res[i] /= A[i][i];
	}

	return res;
}

template<class T>
Vector<T> Substitution::DirectColumnSubstitution(const Matrix<T>& A, const Vector<T>& F)
{
	Vector<T> res(F);

	for (int j = 0; j < F.size(); ++j)
	{
		if (std::abs(A[j][j]) < CONST::EPS)
			throw std::exception("Division by 0");

		res[j] /= A[j][j];

		for (int i = j + 1; i < F.size(); ++i)
		{
			res[i] -= A[i][j] * res[j];
		}
	}

	return res;
}

template<class T>
Vector<T> Substitution::BackRowSubstitution(const Matrix<T>& A, const Vector<T>& F)
{
	Vector<T> res(F);

	for (int i = F.size() - 1; i >= 0; --i)
	{
		if (std::abs(A[i][i]) < CONST::EPS)
			throw std::exception("Division by 0");

		for (int j = i + 1; j < F.size(); ++j)
		{
			res[i] -= A[i][j] * res[j];
		}

		res[i] /= A[i][i];
	}

	return res;
}

template<class T>
Vector<T> Substitution::BackColumnSubstitution(const Matrix<T>& A, const Vector<T>& F)
{
	Vector<T> res(F);

	for (int j = F.size() - 1; j >= 0; --j)
	{
		if (std::abs(A[j][j]) < CONST::EPS)
			throw std::exception("Division by 0");

		res[j] /= A[j][j];

		for (int i = j - 1; i >= 0; --i)
		{
			res[i] -= A[i][j] * res[j];
		}
	}

	return res;
}