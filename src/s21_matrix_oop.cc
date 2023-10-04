#include "s21_matrix_oop.h"

S21Matrix::S21Matrix() noexcept {
  _rows = 0;
  _cols = 0;
  _matrix = nullptr;
}

S21Matrix::S21Matrix(int rows, int cols) {
  if (rows < 1 || cols < 1) {
    throw std::out_of_range("ERROR: incorrect matrix");
  }
  _rows = rows;
  _cols = cols;
  _matrix = new double*[rows]();
  for (int i = 0; i < rows; i++) {
    _matrix[i] = new double[cols]();
  }
 ZeroMatrix();
}

S21Matrix::S21Matrix(const S21Matrix& other)
    : _rows(other._rows), _cols(other._cols) {
  _matrix = new double*[_rows]();
  for (int i = 0; i < _rows; i++) {
    _matrix[i] = new double[_cols]();
  }
  for (int i = 0; i < other._rows; i++) {
    for (int j = 0; j < other._cols; j++) {
      _matrix[i][j] = other._matrix[i][j];
    }
  }
}

void S21Matrix::ZeroMatrix() noexcept {
  for (int i = 0; i < _rows; i++) {
    for (int j = 0; j < _cols; j++) {
      _matrix[i][j] = 0;
    }
  }
}

S21Matrix S21Matrix::Transpose() const{
  if (_rows < 1 || _cols < 1) {
    throw std::out_of_range("ERROR: incorrect matrix");
  }
  S21Matrix tmp(this->_cols, this->_rows);
  for (int i = 0; i < tmp._rows; i++) {
    for (int j = 0; j < tmp._cols; j++) {
      tmp._matrix[i][j] = this->_matrix[j][i];
    }
  }
  return tmp;
}

bool S21Matrix::EqMatrix(const S21Matrix& other) const noexcept {
  bool flag = 1;
  if (this->_cols != other._cols || this->_rows != other._rows) {
    flag = 0;
  }
  for (int i = 0; i < other._rows && flag != 0; i++) {
    for (int j = 0; j < other._cols; j++) {
      if (fabs(_matrix[i][j] - other._matrix[i][j]) > 1e-6) {
        flag = 0;
      }
    }
  }
  return flag;
}

void S21Matrix::SumMatrix(const S21Matrix& other) {
  if (_rows != other._rows || _cols != other._cols) {
    throw std::out_of_range(
        "Incorrect input, matrices should have the same size");
  }
  if (_rows < 1 || _cols < 1 || other._rows < 1 || other._cols < 1) {
    throw std::out_of_range("ERROR: incorrect matrix");
  }
  for (int i = 0; i < other._rows; i++) {
    for (int j = 0; j < other._cols; j++) {
      _matrix[i][j] += other._matrix[i][j];
    }
  }
}

void S21Matrix::SubMatrix(const S21Matrix& other) {
  if (_rows != other._rows || _cols != other._cols) {
    throw std::out_of_range(
        "Incorrect input, matrices should have the same size");
  }
  if (_rows < 1 || _cols < 1 || other._rows < 1 || other._cols < 1) {
    throw std::out_of_range("ERROR: incorrect matrix");
  }
  for (int i = 0; i < other._rows; i++) {
    for (int j = 0; j < other._cols; j++) {
      _matrix[i][j] -= other._matrix[i][j];
    }
  }
}

void S21Matrix::MulNumber(const double num) {
  if (_rows < 1 || _cols < 1) {
    throw std::out_of_range("ERROR: incorrect matrix");
  }
  for (int i = 0; i < this->_rows; i++) {
    for (int j = 0; j < this->_cols; j++) {
      this->_matrix[i][j] = this->_matrix[i][j] * num;
    }
  }
}

void S21Matrix::MulMatrix(const S21Matrix& other) {
  if (_rows < 1 || _cols < 1 || other._rows < 1 || other._cols < 1) {
    throw std::out_of_range("ERROR: incorrect matrix");
  }
  if (_cols != other._rows) {
    throw std::out_of_range(
        "Incorrect input, matrices should have the same size");
  }
  S21Matrix tmp(_rows, _cols);
  for (int i = 0; i < this->_rows; i++) {
    for (int j = 0; j < other._cols; j++) {
      double res = 0;
      for (int k = 0; k < this->_cols; k++) {
        res += this->_matrix[i][k] * other._matrix[k][j];
      }
      tmp._matrix[i][j] = res;
    }
  }
  *this = tmp;
  tmp.Remove();
}

// operator

S21Matrix S21Matrix::operator+(const S21Matrix& other) const {
  S21Matrix res(*this);
  res.SumMatrix(other);
  return res;
}

S21Matrix S21Matrix::operator-(const S21Matrix& other) const {
  S21Matrix res(*this);
  res.SubMatrix(other);
  return res;
}

S21Matrix S21Matrix::operator*(const S21Matrix& other) const {
  S21Matrix res(*this);
  res.MulMatrix(other);
  return res;
}

S21Matrix S21Matrix::operator*(const double num) const {
  S21Matrix res(*this);
  res.MulNumber(num);
  return res;
}

S21Matrix& S21Matrix::operator=(const S21Matrix& other) {
  this->CopyMatrix(other);
  return *this;
}

S21Matrix& S21Matrix::operator+=(const S21Matrix& other) {
  this->SumMatrix(other);
  return *this;
}

S21Matrix& S21Matrix::operator-=(const S21Matrix& other) {
  this->SubMatrix(other);
  return *this;
}

S21Matrix& S21Matrix::operator*=(const S21Matrix& other) {
  this->MulMatrix(other);
  return *this;
}

S21Matrix& S21Matrix::operator*=(const double num) {
  this->MulNumber(num);
  return *this;
}

bool S21Matrix::operator==(const S21Matrix& other) const {
  bool res = EqMatrix(other);
  return res;
}

S21Matrix S21Matrix::CalcComplements() const {
  if (_rows < 1 || _cols < 1) {
    throw std::out_of_range("ERROR: incorrect matrix");
  }
  if (_cols != _rows) {
    throw std::out_of_range("ERROR: matrix is not square");
  }
  S21Matrix result(*this);
  S21Matrix minor(*this);
  for (int i = 0; i < _rows; i++) {
    for (int j = 0; j < _cols; j++) {
      S21Matrix tmp(_rows - 1, _cols - 1);
      minor.Minor(i, j, tmp);
      double det = 0;
      det = tmp.Determinant();
      result._matrix[i][j] = pow(-1.0, i + j) * det;
      tmp.Remove();
    }
  }
  minor.Remove();
  return result;
}

S21Matrix S21Matrix::InverseMatrix() const {
  if (_rows < 1 || _cols < 1) {
    throw std::out_of_range("ERROR: incorrect matrix");
  }
  double det = 0;
  S21Matrix tmp(*this);
  det = this->Determinant();
  if (fabs(det) < 1e-7) {
    tmp.Remove();
    throw std::out_of_range("ERROR: calculation impossible: Determinant = 0");
  } else {
    if (_rows > 1 && _cols > 1) {
      S21Matrix tmp2 = tmp.CalcComplements();
      S21Matrix transpon = tmp2.Transpose();
      transpon.MulNumber(1.0 / det);
      tmp.CopyMatrix(transpon);
      transpon.Remove();
      tmp2.Remove();
    } else {
      S21Matrix tmp2(1, 1);
      tmp._matrix[0][0] = 1.0 / _matrix[0][0];
    }
  }
  return tmp;
}

double S21Matrix::Determinant() const {
  if (_rows < 1 || _cols < 1) {
    throw std::out_of_range("ERROR: incorrect matrix");
  }
  if (_cols != _rows) {
    throw std::out_of_range("ERROR: matrix is not square");
  }
  double res = 0;
  S21Matrix tmp;
  S21Matrix copy(*this);
  res = (-copy.Triangulate(tmp) % 2) ? -1 : 1;
  for (int x = 0; x < tmp._cols; x++) {
    res *= tmp._matrix[x][x];
  }
  tmp.Remove();
  copy.Remove();
  return -res;
}

int S21Matrix::Triangulate(S21Matrix& other) noexcept {
  int znak = 1;
  other.CopyMatrix(*this);
  for (int i = 0; i < other._cols; i++) {
    if (!other._matrix[i][i]) {
      for (int j = i + 1; j < other._rows; j++) {
        if (other._matrix[j][i]) {
          znak -= other.Sdvig(j, i);
        }
      }
    }
    for (int j = i + 1; j < _rows; j++) {
      other.SumStr(j, i, -other._matrix[j][i] / other._matrix[i][i]);
    }
  }
  return -znak;
}

int S21Matrix::Sdvig(int i, int j) noexcept {
  double* tmp = _matrix[i];
  _matrix[i] = _matrix[j];
  _matrix[j] = tmp;
  return -1;
}

void S21Matrix::SumStr(int row, int ro, double tmp) noexcept {
  for (int x = 0; x < _cols; x++) _matrix[row][x] += _matrix[ro][x] * tmp;
}

void S21Matrix::Minor(int x, int y, S21Matrix& other) noexcept {
  int i1 = 0, i2 = 0, j1 = 0, j2 = 0;
  for (i1 = 0; i1 < _rows - 1; i1++) {
    if (i1 == x) {
      i2 = 1;
    }
    j2 = 0;
    for (j1 = 0; j1 < _rows - 1; j1++) {
      if (j1 == y) {
        j2 = 1;
      }
      other._matrix[i1][j1] = _matrix[i1 + i2][j1 + j2];
    }
  }
}

S21Matrix::S21Matrix(S21Matrix&& other) noexcept
    : _rows(other._rows), _cols(other._cols) {
  if (this != &other) {
    _matrix = other._matrix;
    other._matrix = nullptr;
    other._rows = 0;
    other._cols = 0;
  }
}

S21Matrix::~S21Matrix() noexcept {
  if (_matrix) {
    for (int i = 0; i < _rows; i++) {
      delete[] _matrix[i];
    }
    delete[] _matrix;
  }
}

double& S21Matrix::operator()(int row, int col) const {
  if (row >= _rows || col >= _cols || row < 0 || col < 0) {
    throw std::out_of_range("Incorrect input, index is out of range");
  } else {
    return _matrix[row][col];
  }
}

void S21Matrix::CopyMatrix(const S21Matrix& other) noexcept {
  if (this->_matrix) {
    for (int i = 0; i < this->_rows; i++) {
      delete[] this->_matrix[i];
    }
    delete[] this->_matrix;
  }
  this->_cols = other._cols;
  this->_rows = other._rows;
  this->_matrix = new double*[this->_rows]();
  for (int i = 0; i < this->_rows; i++) {
    _matrix[i] = new double[this->_cols]();
  }
  for (int i = 0; i < other._rows; i++) {
    for (int j = 0; j < other._cols; j++) {
      this->_matrix[i][j] = other._matrix[i][j];
    }
  }
}

void S21Matrix::Remove() noexcept {
  if (_matrix) {
    for (int i = 0; i < _rows; i++) {
      if (_matrix[i]) {
        delete[] _matrix[i];
      }
    }
    delete[] _matrix;
    _matrix = nullptr;
  }
  _rows = 0;
  _cols = 0;
}

int S21Matrix::GetRows() const noexcept { return _rows; }

int S21Matrix::GetCols() const noexcept { return _cols; }

void S21Matrix::SetRows(const int rows) {
  if (rows < 1) {
    throw std::out_of_range("ERROR: incorrect matrix");
  }
  S21Matrix tmp(rows, _cols);
  tmp.ZeroMatrix();
  int a = 0;
  if (tmp._rows < _rows) {
    a = tmp._rows;
  } else {
    a = _rows;
  }
  for (int i = 0; i < a; i++) {
    for (int j = 0; j < tmp._cols; j++) {
      tmp._matrix[i][j] = _matrix[i][j];
    }
  }
  CopyMatrix(tmp);
  tmp.Remove();
}

void S21Matrix::SetCols(const int cols) {
  if (cols < 1) {
    throw std::out_of_range("ERROR: incorrect matrix");
  }
  S21Matrix tmp(_rows, cols);
  tmp.ZeroMatrix();
  int a = 0;
  if (tmp._cols < _cols) {
    a = tmp._cols;
  } else {
    a = _cols;
  }
  for (int i = 0; i < tmp._rows; i++) {
    for (int j = 0; j < a; j++) {
      tmp._matrix[i][j] = _matrix[i][j];
    }
  }
  CopyMatrix(tmp);
  tmp.Remove();
}
