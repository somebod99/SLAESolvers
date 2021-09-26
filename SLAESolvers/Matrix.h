#pragma once
#include "Vector.h"

template<class T>
class Matrix
{
private:
	int N;
	int M;

public:
	T** elem;

	Matrix();
	Matrix(const int&, const int&);
	~Matrix();

	Matrix(const Matrix<T>&);
	Matrix<T>& operator=(const Matrix<T>&);

	int getRowCount() const;
	int getColumnCount() const;

	T*& operator[](const int&);
	const T* operator[](const int&) const;

	Matrix<T> operator+(const Matrix<T>&) const;
	Matrix<T> operator-(const Matrix<T>&) const;
	void operator-();
	Matrix<T> operator*(const Matrix<T>&) const;
	Matrix<T> operator*(const T&) const;
	Vector<T> operator*(const Vector<T>&) const;
	void operator+=(const Matrix<T>&);
	void operator-=(const Matrix<T>&);
	void operator*=(const T&);

	bool operator==(const Matrix<T>&) const;

	Vector<T> getColumn(const int&);
	Matrix<T> TransposeMatrix();
	Vector<T> MultiplicationTransMatrixVector(const Vector<T>&);
	void ConsolePrint() const;
};

template<class T>
Matrix<T>::Matrix() : N(0), M(0), elem(nullptr) {}

template<class T>
Matrix<T>::Matrix(const int& n, const int& m) : N(n), M(m)
{
	elem = new T * [N];

	for (int i = 0; i < N; ++i)
	{
		elem[i] = new T[M];

		for (int j = 0; j < M; ++j)
			elem[i][j] = 0;
	}
}

template<class T>
Matrix<T>::~Matrix()
{
	for (int i = 0; i < N; ++i)
		delete[] elem[i];

	delete[] elem;
}

template<class T>
Matrix<T>::Matrix(const Matrix<T>& copy) : N(copy.N), M(copy.M)
{
	elem = new T * [N];

	for (int i = 0; i < N; ++i)
	{
		elem[i] = new T[M];

		for (int j = 0; j < M; ++j)
		{
			elem[i][j] = copy.elem[i][j];
		}
	}
}

template<class T>
Matrix<T>& Matrix<T>::operator=(const Matrix<T>& copy)
{
	if (elem == nullptr)
	{
		N = copy.N;
		M = copy.M;

		elem = new T * [N];

		for (int i = 0; i < N; ++i)
			elem[i] = new T[M];
	}
	else if (N != copy.N || M != copy.M)
	{
		throw std::exception("Dimensional mismatch");
	}

	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < M; ++j)
		{
			elem[i][j] = copy.elem[i][j];
		}
	}

	return *this;
}

template<class T>
int Matrix<T>::getRowCount() const
{
	return N;
}

template<class T>
int Matrix<T>::getColumnCount() const
{
	return M;
}

template<class T>
T*& Matrix<T>::operator[](const int& i)
{
	return elem[i];
}

template<class T>
const T* Matrix<T>::operator[](const int& i) const
{
	return elem[i];
}

template<class T>
Matrix<T> Matrix<T>::operator+(const Matrix<T>& rhs) const
{
	if (N != rhs.N || M != rhs.M)
		throw std::exception("Dimensional mismatch");

	Matrix<T> res(N, M);

	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < M; ++j)
		{
			res.elem[i][j] = elem[i][j] + rhs.elem[i][j];
		}
	}

	return res;
}

template<class T>
Matrix<T> Matrix<T>::operator-(const Matrix<T>& rhs) const
{
	if (N != rhs.N || M != rhs.M)
		throw std::exception("Dimensional mismatch");

	Matrix<T> res(M, N);

	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < M; ++j)
		{
			res[i][j] = elem[i][j] - rhs.elem[i][j];
		}
	}

	return res;
}

template<class T>
void Matrix<T>::operator-()
{
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < M; ++j)
			elem[i][j] = -elem[i][j];
	}
}

template<class T>
Matrix<T> Matrix<T>::operator*(const Matrix<T>& rhs) const
{
	if (M != rhs.N)
		throw std::exception("Dimensional mismatch");

	Matrix<T> res(N, rhs.M);

	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < rhs.M; ++j)
		{
			for (int k = 0; k < M; ++k)
			{
				res[i][j] += elem[i][k] * rhs.elem[k][j];
			}
		}
	}

	return res;
}

template<class T>
Matrix<T> Matrix<T>::operator*(const T& x) const
{
	Matrix<T> res(N, M);

	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < M; ++j)
		{
			res[i][j] = elem[i][j] * x;
		}
	}

	return res;
}

template<class T>
Vector<T> Matrix<T>::operator*(const Vector<T>& rhs) const
{
	if (M != rhs.N)
		throw std::exception("Dimensional mismatch");

	Vector<T> res(N);

	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < M; ++j)
		{
			res[i] += elem[i][j] * rhs.elem[j];
		}
	}

	return res;
}

template<class T>
void Matrix<T>::operator+=(const Matrix<T>& rhs)
{
	if (N != rhs.N || M != rhs.M)
		throw std::exception("Dimensional mismatch");

	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < M; ++j)
		{
			elem[i][j] += rhs.elem[i][j];
		}
	}
}

template<class T>
void Matrix<T>::operator-=(const Matrix<T>& rhs)
{
	if (N != rhs.N || M != rhs.M)
		throw std::exception("Dimensional mismatch");

	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < M; ++j)
		{
			elem[i][j] -= rhs.elem[i][j];
		}
	}
}

template<class T>
void Matrix<T>::operator*=(const T& x)
{
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < M; ++j)
		{
			elem[i][j] *= x;
		}
	}
}

template<class T>
bool Matrix<T>::operator==(const Matrix<T>& rhs) const
{
	if (N != rhs.N || M != rhs.M)
		return false;

	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < M; ++j)
		{
			if (elem[i][j] != rhs.elem[i][j])
				return false;
		}
	}

	return true;
}

template<class T>
Vector<T> Matrix<T>::getColumn(const int& num)
{
	Vector<T> res(N);

	for (int i = 0; i < N; ++i)
	{
		res[i] = elem[i][num];
	}

	return res;
}

template<class T>
Matrix<T> Matrix<T>::TransposeMatrix()
{
	Matrix<T> res(N, M);

	for (int i = 0; i < N; ++i)
	{
		for (int j = i; j < M; ++j)
		{
			res.elem[i][j] = elem[j][i];
		}
	}

	return res;
}

template<class T>
Vector<T> Matrix<T>::MultiplicationTransMatrixVector(const Vector<T>& rhs)
{
	if (M != rhs.size())
		throw std::exception("Dimensional mismatch");

	Vector<T> res(N);

	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < M; ++j)
		{
			res.elem[i][j] += elem[j][i] * rhs[j];
		}
	}

	return res;
}

template<class T>
void Matrix<T>::ConsolePrint() const
{
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < M; ++j)
		{
			std::cout << elem[i][j] << '\t';
		}

		std::cout << '\n';
	}
}
