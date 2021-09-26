#pragma once
#include <exception>
#include <cmath>

template<class T>
class Vector
{
private:
	int N;

public:
	T* elem;

	Vector();
	Vector(const int&);
	~Vector();

	Vector(const Vector<T>&);
	Vector<T>& operator=(const Vector<T>&);

	int size() const;

	T& operator[](const int&);
	const T& operator[](const int&) const;

	Vector<T> operator+(const Vector<T>&) const;
	Vector<T> operator-(const Vector<T>&) const;
	void operator-();
	T operator*(const Vector<T>&) const;
	Vector<T> operator*(const T&) const;
	void operator+=(const Vector<T>&);
	void operator-=(const Vector<T>&);
	void operator*=(const T&);

	bool operator==(const Vector<T>&) const;

	double Norma() const;
	void Normalizing();
	void ConsolePrint() const;
};

template<class T>
Vector<T>::Vector() : N(0), elem(nullptr) {}

template<class T>
Vector<T>::Vector(const int& n) : N(n)
{
	elem = new T[N];

	for (int i = 0; i < N; ++i)
		elem[i] = 0;
}

template<class T>
Vector<T>::~Vector()
{
	delete[] elem;
}

template<class T>
Vector<T>::Vector(const Vector<T>& copy) : N(copy.N)
{
	elem = new T[N];

	for (int i = 0; i < N; ++i)
		elem[i] = copy.elem[i];
}

template<class T>
Vector<T>& Vector<T>::operator=(const Vector<T>& copy)
{
	if (elem == nullptr)
	{
		N = copy.N;
		elem = new T[N];
	}
	else if (N != copy.N)
	{
		throw std::exception("Dimensional mismatch");
	}

	for (int i = 0; i < N; ++i)
		elem[i] = copy.elem[i];

	return *this;
}

template<class T>
int Vector<T>::size() const
{
	return N;
}

template<class T>
T& Vector<T>::operator[](const int& i)
{
	if (i >= N || i < 0)
		throw std::exception("Array out of bounds");

	return elem[i];
}

template<class T>
const T& Vector<T>::operator[](const int& i) const
{
	if (i >= N || i < 0)
		throw std::exception("Array out of bounds");

	return elem[i];
}

template<class T>
Vector<T> Vector<T>::operator+(const Vector<T>& rhs) const
{
	if (rhs.N != N)
		throw std::exception("Dimensional mismatch");

	Vector<T> res(N);

	for (int i = 0; i < N; ++i)
		res[i] = elem[i] + rhs.elem[i];

	return res;
}

template<class T>
Vector<T> Vector<T>::operator-(const Vector<T>& rhs) const
{
	if (rhs.N != N)
		throw std::exception("Dimensional mismatch");

	Vector<T> res(N);

	for (int i = 0; i < N; ++i)
		res[i] = elem[i] - rhs.elem[i];

	return res;
}

template<class T>
void Vector<T>::operator-()
{
	for (int i = 0; i < N; ++i)
		elem[i] = -elem[i];
}

template<class T>
T Vector<T>::operator*(const Vector<T>& rhs) const
{
	if (rhs.N != N)
		throw std::exception("Dimensional mismatch");

	T res = 0;

	for (int i = 0; i < N; ++i)
		res += elem[i] * rhs.elem[i];

	return res;
}

template<class T>
Vector<T> Vector<T>::operator*(const T& x) const
{
	Vector<T> res(N);

	for (int i = 0; i < N; ++i)
		res[i] = elem[i] * x;

	return res;
}

template<class T>
void Vector<T>::operator+=(const Vector<T>& rhs)
{
	if (rhs.N != N)
		throw std::exception("Dimensional mismatch");

	for (int i = 0; i < N; ++i)
		elem[i] += rhs.elem[i];
}

template<class T>
void Vector<T>::operator-=(const Vector<T>& rhs)
{
	if (rhs.N != N)
		throw std::exception("Dimensional mismatch");

	for (int i = 0; i < N; ++i)
		elem[i] -= rhs.elem[i];
}

template<class T>
void Vector<T>::operator*=(const T& x)
{
	for (int i = 0; i < N; ++i)
		elem[i] *= x;
}

template<class T>
bool Vector<T>::operator==(const Vector<T>& rhs) const
{
	if (rhs.N != N)
		throw std::exception("Dimensional mismatch");

	for (int i = 0; i < N; ++i)
		if (elem[i] != rhs.elem[i])
			return false;

	return true;
}

template<class T>
double Vector<T>::Norma() const
{
	return std::sqrt(this->operator*(*this));
}

template<class T>
void Vector<T>::Normalizing()
{
	double norm = Norma();

	for (int i = 0; i < N; ++i)
		elem[i] /= norm;
}

template<class T>
void Vector<T>::ConsolePrint() const
{
	for (int i = 0; i < N; ++i)
		std::cout << elem[i] << '\n';
}
