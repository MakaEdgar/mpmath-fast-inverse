#include "mpinv.h"


#include <cmath>
#include <iostream>
#include <fstream>
#include <cassert>
#include <utility>
#include <new>
#include <string>

#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/eigen.hpp>
#include <Eigen/Dense>


const int MP_PRECISION = 50;
typedef boost::multiprecision::
number<boost::multiprecision::backends::cpp_dec_float<MP_PRECISION>> doubleMP;
typedef Eigen::Matrix<doubleMP, Eigen::Dynamic, Eigen::Dynamic> MatrixMP;


template <typename MatrixType>
inline typename MatrixType::Scalar _logdet(const MatrixType& M, bool use_cholesky = false) {
	using namespace Eigen;
	using std::log;
	typedef typename MatrixType::Scalar Scalar;
	Scalar ld = 0;
	if (use_cholesky) {
		LLT<Matrix<Scalar, Dynamic, Dynamic>> chol(M);
		auto& U = chol.matrixL();
		for (unsigned i = 0; i < M.rows(); ++i)
			ld += log(U(i, i));
		ld *= 2;
	}
	else {
		PartialPivLU<Matrix<Scalar, Dynamic, Dynamic>> lu(M);
		auto& LU = lu.matrixLU();
		Scalar c = lu.permutationP().determinant(); // -1 or 1
		for (unsigned i = 0; i < LU.rows(); ++i) {
			const auto& lii = LU(i, i);
			if (lii < Scalar(0)) c *= -1;
			ld += log(abs(lii));
		}
		ld += log(c);
	}
	return ld;
}

doubleMP _logdet_LLT(const MatrixMP& m, bool is_symmetric) {
	assert(is_symmetric && "Only symmetric matrix allowed in _logdet_LLT");
	doubleMP logdet_val = 0;

	Eigen::LLT<MatrixMP> llt(m);
	auto& U = llt.matrixU();

	for (int i = 0; i < m.rows(); ++i)
		logdet_val += log(U(i,i));

	return 2*logdet_val;

	//return D.log().trace();
}

doubleMP _logdet_LDLT(const MatrixMP& m, bool is_symmetric) {
	assert(is_symmetric && "Only symmetric matrix allowed in _logdet_LDLT");
	doubleMP logdet_val = 0;
	Eigen::LDLT<MatrixMP> ldlt(m);

	const auto& D = ldlt.vectorD();
	for (int i = 0; i < m.rows(); ++i)
		logdet_val += log(D(i));

	return logdet_val;
	//return D.log().trace();
}







MatrixMP& _get_matrix(void* ptr_mpmat) {
	assert(ptr_mpmat && "matrix is not initialized");
	return *((MatrixMP*)ptr_mpmat);
}




mpmat::mpmat() {
	ptr_mpmat = nullptr;
	ptr_mpmat_inv = nullptr;
	int cols = 0;
	int rows = 0;
}
mpmat::mpmat(const char* file_in) {
	std::ifstream ifs(file_in);
	if (ifs.is_open()) {
		ifs >> rows >> cols;

		MatrixMP* ptr_matrix = new MatrixMP(rows, cols);
		ptr_mpmat = (void*)ptr_matrix;

		for (int row = 0; row != rows; ++row)
			for (int col = 0; col != cols; ++col)
			{
				ifs >> (*ptr_matrix)(row, col);

			}
	}
	ifs.close();

}
mpmat::mpmat(int cols_, int rows_) : cols(cols_), rows(rows_) {
	MatrixMP* ptr_matrix = new MatrixMP(cols, rows);
	ptr_mpmat = (void*)ptr_matrix;
}


mpmat::~mpmat() {
	clear_ptr_mpmat(this->ptr_mpmat);
	clear_ptr_mpmat(this->ptr_mpmat_inv);

	// this->mpmat::mpmat();    // deprecated in g++

	new (this) mpmat();
}



void mpmat::set_matrix_coeff(int i, int j, const char* doubleMP_str) {
	_get_matrix(ptr_mpmat)(i, j) = doubleMP(doubleMP_str);
}

void _get_matrix_coeff(int i, int j, void* ptr_mpmat, char* coeff_chars) {
	doubleMP& d = _get_matrix(ptr_mpmat)(i, j);
	std::string d_str = d.str();
	//const char* d_chars = &d_str[0];
	strcpy(coeff_chars, &d_str[0]);
}

const char* mpmat::get_matrix_coeff(int i, int j) {
	_get_matrix_coeff(i, j, ptr_mpmat, _curr_coeff_chars);
	return _curr_coeff_chars;
}
const char* mpmat::get_matrix_inv_coeff(int i, int j) {
	_get_matrix_coeff(i, j, ptr_mpmat_inv, _curr_coeff_chars);
	return _curr_coeff_chars;
}




void mpmat::clear_ptr_mpmat(void* ptr_mpmat) {
	if (ptr_mpmat) {
		delete (MatrixMP*)ptr_mpmat;
	}
}

void mpmat::load_mpmat(const char* file_in) {
	this->~mpmat();
	//this->mpmat::mpmat(file_in);    // deprecated in g++
	new (this) mpmat(file_in);
}

void mpmat::save_mpmat(const char* file_out) {
	std::ofstream ofs(file_out);
	ofs << cols << " " << rows << std::endl;
	ofs << std::setprecision(std::numeric_limits<doubleMP>::digits10);
	ofs << _get_matrix(ptr_mpmat) << std::endl;
}

void mpmat::save_mpmat_inv(const char* file_out) {
	std::ofstream ofs(file_out);
	ofs << cols << " " << rows << std::endl;
	ofs << std::setprecision(std::numeric_limits<doubleMP>::digits10);
	ofs << _get_matrix(ptr_mpmat_inv) << std::endl;
}


void mpmat::calc_logdet1(bool is_symmetric) {
	logdet = _logdet(_get_matrix(ptr_mpmat), is_symmetric).convert_to<double>();
}
void mpmat::calc_logdet2(bool is_symmetric) {
	logdet = _logdet_LLT(_get_matrix(ptr_mpmat), is_symmetric).convert_to<double>();
}
void mpmat::calc_logdet3(bool is_symmetric) {
	logdet = _logdet_LDLT(_get_matrix(ptr_mpmat), is_symmetric).convert_to<double>();
}




void mpmat::calc_logdet(bool is_symmetric) {
	if (is_symmetric) {
		logdet = _logdet_LDLT(_get_matrix(ptr_mpmat), is_symmetric).convert_to<double>();
	}
	else {
		logdet = _logdet(_get_matrix(ptr_mpmat), false).convert_to<double>();
	}

}
void mpmat::calc_invert() {
	clear_ptr_mpmat(ptr_mpmat_inv);

	MatrixMP* ptr_matrix_inv = new MatrixMP(_get_matrix(ptr_mpmat).inverse());
	ptr_mpmat_inv = (void*)ptr_matrix_inv;
}

std::pair<doubleMP, MatrixMP> _calc_invert_with_logdet(MatrixMP& m, bool is_symmetric = false) {
	assert(is_symmetric && "Only symmetric matrix allowed in _calc_invert_with_logdet");

	Eigen::LDLT<MatrixMP> ldlt(m);
	MatrixMP inv1 = ldlt.solve(MatrixMP::Identity(m.cols(), m.rows()));

	doubleMP logdet_val = 0;
	const auto& D = ldlt.vectorD();
	for (int i = 0; i < m.rows(); ++i)
		logdet_val += log(D(i));

	return std::make_pair(logdet_val, inv1);

}

void mpmat::calc_invert_with_logdet(bool is_symmetric) {

	clear_ptr_mpmat(ptr_mpmat_inv);
	
	auto logdet_and_invmat_pair = _calc_invert_with_logdet(_get_matrix(ptr_mpmat), is_symmetric);

	logdet = logdet_and_invmat_pair.first.convert_to<double>();
	
	MatrixMP* ptr_matrix_inv = new MatrixMP(logdet_and_invmat_pair.second);
	ptr_mpmat_inv = (void*)ptr_matrix_inv;
}

double mpmat::get_logdet() {
	return logdet;
}

void mpmat::get_data(const char* data) {
	//std::cout << typeid(data).name() << std::endl;
	std::cout << data << std::endl;
	doubleMP d = doubleMP(data);
	std::cout << d << std::endl;
	std::cout << 2*d << std::endl;
}












//
//MatrixMP& _get_matrix(mpmat& m) {
//	assert(m.ptr_mpmat && "matrix is not initialized");
//	return *((MatrixMP*)m.ptr_mpmat);
//}
//
//
//void _set_matrix(MatrixMP& m_matrix, mpmat& m_mpmat) {
//	m_mpmat.~mpmat();
//
//	MatrixMP* ptr_matrix = new MatrixMP(m_matrix);
//	m_mpmat.ptr_mpmat = (void*)ptr_matrix;
//	m_mpmat.ptr_mpmat_inv = nullptr;
//	m_mpmat.cols = m_matrix.cols();
//	m_mpmat.rows = m_matrix.rows();
//}


//
//mpmat& mpmat::operator=(mpmat& m) {
//	ptr_mpmat = m.ptr_mpmat;
//	ptr_mpmat_inv = m.ptr_mpmat_inv;
//
//	cols = m.cols;
//	rows = m.rows;
//	return *this;
//}
