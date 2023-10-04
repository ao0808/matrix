#ifndef __S21MATRIX_H__
#define __S21MATRIX_H__

#include <cmath>
#include <cstring>
#include <iostream>

class S21Matrix {
 private:
  // attributes
  int _rows, _cols;  // rows and columns attributes
  double** _matrix;       // pointer to the memory where the matrix will be allocated
  void ZeroMatrix() noexcept;
  void CopyMatrix(const S21Matrix& other) noexcept;
  int Sdvig(int i, int j) noexcept;
  void SumStr(int row, int ro, double tmp) noexcept;
  void Minor(int i, int j, S21Matrix& other) noexcept;
  void Remove() noexcept;
  int Triangulate(S21Matrix& other) noexcept;

 public:
  S21Matrix() noexcept;                   // default constructor
  S21Matrix(int rows, int cols);          // parameterized constructor
  S21Matrix(const S21Matrix& other);      // copy constructor
  S21Matrix(S21Matrix&& other) noexcept;  // move constructor
  ~S21Matrix() noexcept;                           // destructor

  bool operator==(const S21Matrix& other) const;
  S21Matrix operator+(const S21Matrix& other) const;
  S21Matrix operator-(const S21Matrix& other) const;
  S21Matrix operator*(const S21Matrix& other) const;
  S21Matrix operator*(const double num) const;
  S21Matrix& operator=(const S21Matrix& other);
  S21Matrix& operator-=(const S21Matrix& other);
  S21Matrix& operator+=(const S21Matrix& other);
  S21Matrix& operator*=(const S21Matrix& other);
  S21Matrix& operator*=(const double num);
  double& operator()(int row, int col) const;
  void SumMatrix(const S21Matrix& other);
  void SubMatrix(const S21Matrix& other);
  void MulNumber(const double num);
  void MulMatrix(const S21Matrix& other);
  bool EqMatrix(const S21Matrix& other) const noexcept;
  S21Matrix Transpose() const;
  double Determinant() const;
  S21Matrix InverseMatrix() const;
  S21Matrix CalcComplements() const;
  int GetRows() const noexcept;
  int GetCols() const noexcept;
  void SetRows(const int rows);
  void SetCols(const int cols);
};

#endif
