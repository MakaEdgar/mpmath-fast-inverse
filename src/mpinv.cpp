// Turn off MSVC error when using strlen 
#define _CRT_SECURE_NO_WARNINGS

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


//const int MP_PRECISION defined in mpinv.h
typedef boost::multiprecision::
number<boost::multiprecision::backends::cpp_dec_float<MP_PRECISION>> doubleMP;
typedef Eigen::Matrix<doubleMP, Eigen::Dynamic, Eigen::Dynamic> MatrixMP;

// Private mpinv.cpp functions
doubleMP _logdet_eigen(const MatrixMP& m) {
	return log(m.determinant());
}
doubleMP _logdet_LU(const MatrixMP& m) {
	Eigen::PartialPivLU<MatrixMP> lu(m);
	auto& LU = lu.matrixLU();
	
	doubleMP logdet_val = doubleMP("0.0");
	doubleMP c = lu.permutationP().determinant(); // -1 or 1
	for (int i = 0; i < LU.rows(); ++i) {
		const auto& lii = LU(i, i);
		if (lii < doubleMP("0.0")) c *= doubleMP("-1.0");
		logdet_val += log(abs(lii));
	}
	logdet_val += log(c);

	return logdet_val;
}
doubleMP _logdet_LLT(const MatrixMP& m) {
	//assert(is_symmetric && "Only symmetric matrix allowed in _logdet_LLT");
	Eigen::LLT<MatrixMP> llt(m);
	auto& U = llt.matrixU();

	doubleMP logdet_val = doubleMP("0.0");
	for (int i = 0; i < m.rows(); ++i)
		logdet_val += log(U(i,i));
	logdet_val *= 2;

	return logdet_val;
}
doubleMP _logdet_LDLT(const MatrixMP& m) {
	//assert(is_symmetric && "Only symmetric matrix allowed in _logdet_LDLT");
	Eigen::LDLT<MatrixMP> ldlt(m);

	doubleMP logdet_val = doubleMP("0.0");
	const auto& D = ldlt.vectorD();
	for (int i = 0; i < m.rows(); ++i)
		logdet_val += log(D(i));

	return logdet_val;
}

doubleMP _det_eigen(const MatrixMP& m) {
	return m.determinant();
}
doubleMP _det_LU(const MatrixMP& m) {
	doubleMP logdet_val = 0;
	Eigen::PartialPivLU<MatrixMP> lu(m);
	auto& LU = lu.matrixLU();

	doubleMP det_val = doubleMP("1.0");
	doubleMP c = lu.permutationP().determinant(); // -1 or 1
	for (int i = 0; i < LU.rows(); ++i) {
		det_val *= LU(i, i);
		//const auto& lii = LU(i, i);
		//if (lii < doubleMP("0.0")) c *= doubleMP("-1.0");
		//det_val *= abs(lii);
	}
	det_val *= c;

	return det_val;
}
doubleMP _det_LLT(const MatrixMP& m) {
	//assert(is_symmetric && "Only symmetric matrix allowed in _logdet_LLT");
	Eigen::LLT<MatrixMP> llt(m);
	auto& U = llt.matrixU();

	doubleMP det_val = doubleMP("1.0");
	for (int i = 0; i < m.rows(); ++i)
		det_val *= U(i, i);
	det_val *= det_val;

	return det_val;
}
doubleMP _det_LDLT(const MatrixMP& m) {
	//assert(is_symmetric && "Only symmetric matrix allowed in _logdet_LDLT");
	Eigen::LDLT<MatrixMP> ldlt(m);

	doubleMP det_val = doubleMP("1.0");
	const auto& D = ldlt.vectorD();
	for (int i = 0; i < m.rows(); ++i)
		det_val *= D(i);

	return det_val;
}

MatrixMP _inverse_eigen(const MatrixMP& m) {
	return m.inverse();
}
MatrixMP _inverse_LU(const MatrixMP& m) {
	Eigen::PartialPivLU<MatrixMP> lu(m);
	return lu.solve(MatrixMP::Identity(m.cols(), m.rows()));
}
MatrixMP _inverse_LLT(const MatrixMP& m) {
	Eigen::LLT<MatrixMP> llt(m);
	return llt.solve(MatrixMP::Identity(m.cols(), m.rows()));
}
MatrixMP _inverse_LDLT(const MatrixMP& m) {
	Eigen::LDLT<MatrixMP> ldlt(m);
	return ldlt.solve(MatrixMP::Identity(m.cols(), m.rows()));
}

std::pair<doubleMP, MatrixMP> _calc_inverse_with_logdet_LU(MatrixMP& m) {
	Eigen::PartialPivLU<MatrixMP> lu(m);
	auto& LU = lu.matrixLU();

	doubleMP logdet_val = doubleMP("0.0");
	doubleMP c = lu.permutationP().determinant(); // -1 or 1
	for (int i = 0; i < LU.rows(); ++i) {
		const auto& lii = LU(i, i);
		if (lii < doubleMP("0.0")) c *= doubleMP("-1.0");
		logdet_val += log(abs(lii));
	}
	logdet_val += log(c);

	MatrixMP m_inv = lu.solve(MatrixMP::Identity(m.cols(), m.rows()));

	return std::make_pair(logdet_val, m_inv);

}
std::pair<doubleMP, MatrixMP> _calc_inverse_with_logdet_LLT(MatrixMP& m) {
	Eigen::LLT<MatrixMP> llt(m);
	auto& U = llt.matrixU();

	doubleMP logdet_val = doubleMP("0.0");
	for (int i = 0; i < m.rows(); ++i)
		logdet_val += log(U(i, i));
	logdet_val *= 2;

	MatrixMP m_inv = llt.solve(MatrixMP::Identity(m.cols(), m.rows()));

	return std::make_pair(logdet_val, m_inv);

}
std::pair<doubleMP, MatrixMP> _calc_inverse_with_logdet_LDLT(MatrixMP& m) {
	Eigen::LDLT<MatrixMP> ldlt(m);

	doubleMP logdet_val = doubleMP("0.0");
	const auto& D = ldlt.vectorD();
	for (int i = 0; i < m.rows(); ++i)
		logdet_val += log(D(i));

	MatrixMP m_inv = ldlt.solve(MatrixMP::Identity(m.cols(), m.rows()));

	return std::make_pair(logdet_val, m_inv);

}


MatrixMP& _get_matrix(void* ptr_mpmat) {
	assert(ptr_mpmat && "matrix is not initialized");
	return *((MatrixMP*)ptr_mpmat);
}
void _save_mpmat(const char* file_out, void* ptr_mpmat) {
	MatrixMP& m = _get_matrix(ptr_mpmat);
	std::ofstream ofs(file_out);
	ofs << m.cols() << " " << m.rows() << std::endl;
	ofs << std::setprecision(std::numeric_limits<doubleMP>::digits10);
	ofs << _get_matrix(ptr_mpmat) << std::endl;
}
void _doubleMP_to_str(doubleMP d, char* d_str) {
	std::string d_string = d.str();
	std::copy(d_string.begin(), d_string.end(), d_str);
	//d_str = d_string.c_str();
	//strcpy(d_str, &d_string[0]);
}
void _get_matrix_coeff(int i, int j, void* ptr_mpmat, char* coeff_str) {
	doubleMP& d = _get_matrix(ptr_mpmat)(i, j);
	_doubleMP_to_str(d, coeff_str);
}
void _clear_ptr_mpmat(void* ptr_mpmat) {
	if (ptr_mpmat) {
		delete (MatrixMP*)ptr_mpmat;
	}
}




// mpinv class public functions
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
	_clear_ptr_mpmat(this->ptr_mpmat);
	_clear_ptr_mpmat(this->ptr_mpmat_inv);

	// this->mpmat::mpmat();    // deprecated in g++

	new (this) mpmat();
}

void mpmat::set_matrix_coeff(int i, int j, const char* doubleMP_str) {
	_get_matrix(ptr_mpmat)(i, j) = doubleMP(doubleMP_str);
}
const char* mpmat::get_matrix_coeff(int i, int j) {
	_get_matrix_coeff(i, j, ptr_mpmat, _curr_coeff_str);
	return _curr_coeff_str;
}
const char* mpmat::get_matrix_inv_coeff(int i, int j) {
	_get_matrix_coeff(i, j, ptr_mpmat_inv, _curr_coeff_str);
	return _curr_coeff_str;
}

void mpmat::load_mpmat(const char* file_in) {
	this->~mpmat();
	//this->mpmat::mpmat(file_in);    // deprecated in g++
	new (this) mpmat(file_in);
}
void mpmat::save_mpmat(const char* file_out) {
	_save_mpmat(file_out, ptr_mpmat);
}
void mpmat::save_mpmat_inv(const char* file_out) {
	_save_mpmat(file_out, ptr_mpmat_inv);
}


void mpmat::calc_logdet_LU() {
	logdet = _logdet_LU(_get_matrix(ptr_mpmat)).convert_to<double>();
}
void mpmat::calc_logdet_LLT() {
	logdet = _logdet_LLT(_get_matrix(ptr_mpmat)).convert_to<double>();
}
void mpmat::calc_logdet_LDLT() {
	logdet = _logdet_LDLT(_get_matrix(ptr_mpmat)).convert_to<double>();
}
void mpmat::calc_logdet() {
	logdet = _logdet_eigen(_get_matrix(ptr_mpmat)).convert_to<double>();
}

void mpmat::calc_det_LU() {
	doubleMP det_ = _det_LU(_get_matrix(ptr_mpmat));
	det = det_.convert_to<double>();
	_doubleMP_to_str(det_, _det_str);
}
void mpmat::calc_det_LLT() {
	doubleMP det_ = _det_LLT(_get_matrix(ptr_mpmat));
	det = det_.convert_to<double>();
	_doubleMP_to_str(det_, _det_str);
}
void mpmat::calc_det_LDLT() {
	doubleMP det_ = _det_LDLT(_get_matrix(ptr_mpmat));
	det = det_.convert_to<double>();
	_doubleMP_to_str(det_, _det_str);
}
void mpmat::calc_det() {
	doubleMP det_ = _det_eigen(_get_matrix(ptr_mpmat));
	det = det_.convert_to<double>();
	_doubleMP_to_str(det_, _det_str);
}

void mpmat::calc_inverse_LU() {
	_clear_ptr_mpmat(ptr_mpmat_inv);
	MatrixMP* ptr_matrix_inv = new MatrixMP(_inverse_LU(_get_matrix(ptr_mpmat)));
	ptr_mpmat_inv = (void*)ptr_matrix_inv;

}
void mpmat::calc_inverse_LLT() {
	_clear_ptr_mpmat(ptr_mpmat_inv);
	MatrixMP* ptr_matrix_inv = new MatrixMP(_inverse_LLT(_get_matrix(ptr_mpmat)));
	ptr_mpmat_inv = (void*)ptr_matrix_inv;

}
void mpmat::calc_inverse_LDLT() {
	_clear_ptr_mpmat(ptr_mpmat_inv);
	MatrixMP* ptr_matrix_inv = new MatrixMP(_inverse_LDLT(_get_matrix(ptr_mpmat)));
	ptr_mpmat_inv = (void*)ptr_matrix_inv;

}
void mpmat::calc_inverse() {
	_clear_ptr_mpmat(ptr_mpmat_inv);
	MatrixMP* ptr_matrix_inv = new MatrixMP(_inverse_eigen(_get_matrix(ptr_mpmat)));
	ptr_mpmat_inv = (void*)ptr_matrix_inv;
}

void mpmat::calc_inverse_with_logdet_LU() {
	_clear_ptr_mpmat(ptr_mpmat_inv);
	auto logdet_and_invmat_pair = _calc_inverse_with_logdet_LU(_get_matrix(ptr_mpmat));

	logdet = logdet_and_invmat_pair.first.convert_to<double>();

	MatrixMP* ptr_matrix_inv = new MatrixMP(logdet_and_invmat_pair.second);
	ptr_mpmat_inv = (void*)ptr_matrix_inv;
}
void mpmat::calc_inverse_with_logdet_LLT() {
	_clear_ptr_mpmat(ptr_mpmat_inv);
	auto logdet_and_invmat_pair = _calc_inverse_with_logdet_LLT(_get_matrix(ptr_mpmat));

	logdet = logdet_and_invmat_pair.first.convert_to<double>();

	MatrixMP* ptr_matrix_inv = new MatrixMP(logdet_and_invmat_pair.second);
	ptr_mpmat_inv = (void*)ptr_matrix_inv;
}
void mpmat::calc_inverse_with_logdet_LDLT() {
	_clear_ptr_mpmat(ptr_mpmat_inv);
	auto logdet_and_invmat_pair = _calc_inverse_with_logdet_LDLT(_get_matrix(ptr_mpmat));

	logdet = logdet_and_invmat_pair.first.convert_to<double>();

	MatrixMP* ptr_matrix_inv = new MatrixMP(logdet_and_invmat_pair.second);
	ptr_mpmat_inv = (void*)ptr_matrix_inv;
}
void mpmat::calc_inverse_with_logdet() {
	calc_inverse_with_logdet_LLT();
}


double mpmat::get_logdet() {
	return logdet;
}
double mpmat::get_det() {
	return det;
}
const char* mpmat::get_det_str() {
	return _det_str;
}

