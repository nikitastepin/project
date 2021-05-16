#include <iostream>
#include <cstdlib>
#include <stdexcept>
 
using MatrixAllocationError = std::bad_alloc;
 
class MatrixWrongSizeError : public std::domain_error {
public:
    MatrixWrongSizeError() : std::domain_error("Matrices sizes mismatch") {
    }
};
 
class MatrixIndexError : public std::out_of_range {
public:
    MatrixIndexError() : std::out_of_range("Matrix element does not exist") {
    }
};
 
template <typename T>
T GetZero() {
    return T(0);
}
 
template <class T>
T** Alloc(int rows, int cols) {
    T* bulk = rows * cols > 0 ? new T[rows * cols] : nullptr;
    T** arr = new T*[rows];
    for (int i = 0; i < rows; ++i) {
        arr[i] = bulk + cols * i;
    }
    return arr;
}
 
template <typename T>
class Matrix {
private:
    int rows_count_ = 0;
    int cols_count_ = 0;
    T** matrix_ = nullptr;
 
    void Copy(const Matrix& matr) {
        Clear();
        rows_count_ = matr.rows_count_;
        cols_count_ = matr.cols_count_;
        matrix_ = Alloc<T>(rows_count_, cols_count_);
        for (int i = 0; i < rows_count_; ++i) {
            for (int j = 0; j < cols_count_; ++j) {
                matrix_[i][j] = matr.matrix_[i][j];
            }
        }
    }
 
public:
    Matrix(int rows, int cols, const T& value = GetZero<T>()) : rows_count_(rows), cols_count_(cols) {
        matrix_ = Alloc<T>(rows, cols);
        for (int i = 0; i < rows_count_; ++i) {
            for (int j = 0; j < cols_count_; ++j) {
                matrix_[i][j] = value;
            }
        }
    }
 
    void Clear() {
        if (matrix_ != nullptr) {
            delete[] matrix_[0];
            delete[] matrix_;
            matrix_ = nullptr;
        }
        rows_count_ = 0;
        cols_count_ = 0;
    }
 
    ~Matrix() {
        Clear();
    }
 
    Matrix(const Matrix& matr) {
        Copy(matr);
    }
 
    Matrix& operator=(const Matrix& matr) {
        if (this != &matr) {
            Copy(matr);
        }
        return *this;
    }
 
    Matrix(Matrix&& matr) {
        rows_count_ = matr.rows_count_;
        cols_count_ = matr.cols_count_;
        matrix_ = matr.matrix_;
        matr.matrix_ = nullptr;
    }
 
    Matrix& operator=(Matrix&& matr) {
        if (this != &matr) {
            Clear();
            rows_count_ = matr.rows_count_;
            cols_count_ = matr.cols_count_;
            matrix_ = matr.matrix_;
            matr.matrix_ = nullptr;
        }
        return *this;
    }
 
    int GetRowsNumber() const {
        return rows_count_;
    }
 
    int GetColumnsNumber() const {
        return cols_count_;
    }
 
    T& operator()(int rowNumber, int colNumber) {
        if (rowNumber < 0 || rowNumber >= rows_count_ || colNumber < 0 || colNumber >= cols_count_) {
            throw MatrixIndexError();
        }
        return matrix_[rowNumber][colNumber];
    }
 
    const T& operator()(int rowNumber, int colNumber) const {
        if (rowNumber < 0 || rowNumber >= rows_count_ || colNumber < 0 || colNumber >= cols_count_) {
            throw MatrixIndexError();
        }
        return matrix_[rowNumber][colNumber];
    }
 
    Matrix GetTransposed() const {
        Matrix res(cols_count_, rows_count_);
 
        for (int i = 0; i < cols_count_; ++i) {
            for (int j = 0; j < rows_count_; ++j) {
                res.matrix_[i][j] = matrix_[j][i];
            }
        }
        return res;
    }
 
    Matrix& Transpose() {
        *this = GetTransposed();
        return *this;
    }
 
    Matrix& operator+=(const Matrix& matr) {
        if (rows_count_ != matr.rows_count_ || cols_count_ != matr.cols_count_) {
            throw MatrixWrongSizeError();
        }
        for (int i = 0; i < rows_count_; ++i) {
            for (int j = 0; j < cols_count_; ++j) {
                matrix_[i][j] += matr.matrix_[i][j];
            }
        }
        return *this;
    }
 
    Matrix& operator-=(const Matrix& matr) {
        if (rows_count_ != matr.rows_count_ || cols_count_ != matr.cols_count_) {
            throw MatrixWrongSizeError();
        }
        for (int i = 0; i < rows_count_; ++i) {
            for (int j = 0; j < cols_count_; ++j) {
                matrix_[i][j] -= matr.matrix_[i][j];
            }
        }
        return *this;
    }
 
    Matrix& operator*=(const T& value) {
        for (int i = 0; i < rows_count_; ++i) {
            for (int j = 0; j < cols_count_; ++j) {
                matrix_[i][j] *= value;
            }
        }
        return *this;
    }
 
    Matrix& operator*=(const Matrix& matr) {
        *this = *this * matr;
        return *this;
    }
 
    template <class Other>
    friend Matrix<Other> operator*(const Matrix<Other>& lhs, const Matrix<Other>& rhs);
};
 
template <class T>
Matrix<T> operator*(const Matrix<T>& matr, const T& scalar) {
    Matrix<T> res = matr;
    res *= scalar;
    return res;
}
 
template <class T>
Matrix<T> operator*(const T& scalar, const Matrix<T>& matr) {
    Matrix<T> res = matr;
    res *= scalar;
    return res;
}
 
template <class T>
Matrix<T> operator+(const Matrix<T>& lhs, const Matrix<T>& rhs) {
    Matrix<T> res = lhs;
    res += rhs;
    return res;
}
 
template <class T>
Matrix<T> operator-(const Matrix<T>& lhs, const Matrix<T>& rhs) {
    Matrix<T> res = lhs;
    res -= rhs;
    return res;
}
 
template <class T>
Matrix<T> operator*(const Matrix<T>& lhs, const Matrix<T>& rhs) {
    if (lhs.cols_count_ != rhs.rows_count_) {
        throw MatrixWrongSizeError();
    }
    Matrix<T> res(lhs.rows_count_, rhs.cols_count_);
    for (int i = 0; i < lhs.rows_count_; ++i) {
        for (int j = 0; j < rhs.cols_count_; ++j) {
            for (int k = 0; k < lhs.cols_count_; ++k) {
                res(i, j) += lhs(i, k) * rhs(k, j);
            }
        }
    }
    return res;
}
 
template <class T>
std::istream& operator>>(std::istream& in, Matrix<T>& matr) {
    for (int i = 0; i < matr.GetRowsNumber(); ++i) {
        for (int j = 0; j < matr.GetColumnsNumber(); ++j) {
            in >> matr(i, j);
        }
    }
    return in;
}
 
template <class T>
std::ostream& operator<<(std::ostream& out, const Matrix<T>& matr) {
    for (int i = 0; i < matr.GetRowsNumber(); ++i) {
        for (int j = 0; j < matr.GetColumnsNumber(); ++j) {
            out << matr(i, j) << ' ';
        }
        out << std::endl;
    }
    return out;
}
 
 
class RationalDivisionByZero : std::range_error {
public:
    RationalDivisionByZero() : std::range_error("Dividing Rational by zero") {
    }
};
 
template <class T>
T GCD(const T& a, const T& b) {
    if (b == 0) {
        return a;
    }
    return GCD(b, a % b);
}
 
class Rational {
private:
    int numerator_;
    int denominator_; // positive
 
    void Reduce() {
        if (denominator_ < 0) {
            denominator_ = -denominator_;
            numerator_ = -numerator_;
        }
        int gcd = GCD(abs(numerator_), denominator_);
        numerator_ /= gcd;
        denominator_ /= gcd;
    }
 
public:
    Rational(int numerator = 0, int denominator = 1) :
        numerator_(numerator), denominator_(denominator) {
        Reduce();
    }
 
    int GetNumerator() const {
        return numerator_;
    }
 
    int GetDenominator() const {
        return denominator_;
    }
 
    Rational& operator+=(const Rational& rhs) {
        numerator_ = numerator_ * rhs.denominator_ + rhs.numerator_ * denominator_;
        denominator_ *= rhs.denominator_;
        Reduce();
        return *this;
    }
 
    Rational& operator-=(const Rational& rhs) {
        numerator_ = numerator_ * rhs.denominator_ - rhs.numerator_ * denominator_;
        denominator_ *= rhs.denominator_;
        Reduce();
        return *this;
    }
 
    Rational& operator*=(const Rational& rhs) {
        numerator_ *= rhs.numerator_;
        denominator_ *= rhs.denominator_;
        Reduce();
        return *this;
    }
 
    Rational& operator/=(const Rational& rhs) {
        if (rhs.numerator_ == 0) {
            throw RationalDivisionByZero();
        }
        numerator_ *= rhs.denominator_;
        denominator_ *= rhs.numerator_;
        Reduce();
        return *this;
    }
 
    Rational& operator++() {
        numerator_ += denominator_;
        return *this;
    }
 
    Rational operator++(int) {
        Rational copy = *this;
        numerator_ += denominator_;
        return copy;
    }
 
    Rational& operator--() {
        numerator_ -= denominator_;
        return *this;
    }
 
    Rational operator--(int) {
        Rational copy = *this;
        numerator_ -= denominator_;
        return copy;
    }
 
    friend Rational operator-(const Rational& number);
    friend bool operator<(const Rational& lhs, const Rational& rhs);
    friend bool operator==(const Rational& lhs, const Rational& rhs);
    friend std::istream& operator>>(std::istream& in, Rational& number);
    friend std::ostream& operator<<(std::ostream& out, const Rational& number);
};
 
Rational operator+(const Rational& number) {
    return number;
}
 
Rational operator-(const Rational& number) {
    return Rational(-number.numerator_, number.denominator_);
}
 
Rational operator+(const Rational& lhs, const Rational& rhs) {
    Rational res = lhs;
    res += rhs;
    return res;
}
 
Rational operator-(const Rational& lhs, const Rational& rhs) {
    Rational res = lhs;
    res -= rhs;
    return res;
}
 
Rational operator*(const Rational& lhs, const Rational& rhs) {
    Rational res = lhs;
    res *= rhs;
    return res;
}
 
Rational operator/(const Rational& lhs, const Rational& rhs) {
    Rational res = lhs;
    res /= rhs;
    return res;
}
 
bool operator<(const Rational& lhs, const Rational& rhs) {
    return lhs.numerator_ * rhs.denominator_ < rhs.numerator_ * lhs.denominator_;
}
 
bool operator<=(const Rational& lhs, const Rational& rhs) {
    return !(rhs < lhs);
}
 
bool operator>(const Rational& lhs, const Rational& rhs) {
    return rhs < lhs;
}
 
bool operator>=(const Rational& lhs, const Rational& rhs) {
    return !(lhs < rhs);
}
 
bool operator==(const Rational& lhs, const Rational& rhs) {
    return lhs.numerator_ == rhs.numerator_ && lhs.denominator_ == rhs.denominator_;
}
 
bool operator!=(const Rational& lhs, const Rational& rhs) {
    return !(lhs == rhs);
}
 
std::istream& operator>>(std::istream& in, Rational& number) {
    const int MAX_SIZE = 20;
    char str[MAX_SIZE];
    in >> str;
 
    if (sscanf(str, "%d/%d", &number.numerator_, &number.denominator_) < 2 || number.denominator_ == 0) {
        number.denominator_ = 1;
    }
    number.Reduce();
 
    return in;
}
 
std::ostream& operator<<(std::ostream& out, const Rational& number) {
    out << number.numerator_;
    if (number.denominator_ != 1) {
        out << "/" << number.denominator_;
    }
    return out;
}
 
//=================== main() ===============//
 
using namespace std;
 
int main() {
    int m, n, p, q;
    cin >> m >> n >> p >> q;
 
    Matrix<int> A(m, n), B(p, q);
    cin >> A >> B;
 
    A = A;
    try {
        cout << A + B * 2 - m * A << endl;
        cout << (A -= B += A *= 2) << endl;
        cout << (((A -= B) += A) *= 2) << endl;
    } catch (const MatrixWrongSizeError&) {
        cout << "A and B are of different size." << endl;
    }
    B = A;
 
    {
        Matrix<int> AA(A);
        Matrix<int> AAA(1, 1);
        AAA = A;
        cout << AA << endl;
        cout << (AAA += Matrix<int>(m, n)) + B << endl;
    }
 
    Rational r;
    cin >> r;
    Matrix<Rational> C(m, n), D(p, q);
    cin >> C >> D;
    try {
        cout << C * D << endl;
        cout << (C *= D) << endl;
        cout << C << endl;
    } catch (const MatrixWrongSizeError&) {
        cout << "C and D have not appropriate sizes for multiplication." << endl;
    }
    cout << C.GetTransposed() * (r * C) << endl;
    cout << C.Transpose() << endl;
    try {
        (C(0, 0) *= 6) /= 3;
        cout << C(0, 0) << endl;
        cout << C(m, m) << endl;
    } catch (const MatrixIndexError&) {
        cout << "Index out of range." << endl;
    }
 
    {
        const Matrix<Rational>& rC = C;
        cout << rC << endl;
        cout << rC.GetRowsNumber() << ' ' << rC.GetColumnsNumber() << ' ' << rC(0, 0) << endl;
        cout << (C = C) * (Rational(1, 2) * rC).GetTransposed() << endl;
    }
    return 0;
}

