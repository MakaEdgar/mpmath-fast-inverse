
#include "mpinv.h"

#include <fstream>
#include <iostream>

#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/eigen.hpp>
#include <Eigen/Dense>

const int MP_PRECISION = 50;
typedef boost::multiprecision::
	number<boost::multiprecision::backends::cpp_dec_float<MP_PRECISION>> doubleMP;
typedef Eigen::Matrix<doubleMP, Eigen::Dynamic, Eigen::Dynamic> MatrixMP;


// set use_cholesky if M is symmetric - it's faster and more stable
// for dep paring it won't be
template <typename MatrixType>
inline typename MatrixType::Scalar logdet_stable(const MatrixType& M, bool use_cholesky = false) {
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




mpinv::mpinv(char* file_in, char* file_out) {
    file_in_ = nullptr;
    file_out_ = nullptr;
    
    set_file_in(file_in);
    set_file_out(file_out);
}

mpinv::~mpinv() {
    delete[] file_in_;
    delete[] file_out_;
}

void mpinv::inv() {
    int N = 0, M = 0;
    MatrixMP A;

    std::ifstream ifs (file_in_);
    if (ifs.is_open()) {
        ifs >> N >> M;
        if (N != M) {std::cout<< "AAAA N != M" << std::endl; return;}
        
        A = MatrixMP::Zero(N,N);
        for (int row = 0; row != N; ++row)
            for (int col = 0; col != N; ++col)
            {
                doubleMP item = doubleMP("0.0");
                ifs >> item;
                A(row, col) = item;
            }
    } else { std::cout<< "AAAA file not open!" << std::endl; return;}
    ifs.close();
    
    MatrixMP A_inv = A.inverse();

    std::ofstream ofs(file_out_);
    ofs << N << " " << N << std::endl;
    ofs << std::setprecision(std::numeric_limits<doubleMP>::digits10);
    ofs << A_inv << std::endl;
    
    return;
}



double mpinv::logdet(bool is_symmetric) {
    int N = 0, M = 0;
    MatrixMP A;

    std::ifstream ifs (file_in_);
    if (ifs.is_open()) {
        ifs >> N >> M;
        if (N != M) {std::cout<< "AAAA N != M" << std::endl; return -1;}
        
        A = MatrixMP::Zero(N,N);
        for (int row = 0; row != N; ++row)
            for (int col = 0; col != N; ++col)
            {
                doubleMP item = doubleMP("0.0");
                ifs >> item;
                A(row, col) = item;
            }
    } else { std::cout<< "AAAA file not open!" << std::endl; return -1;}
    ifs.close();
    
    double logdet_A = logdet_stable(A, is_symmetric).convert_to<double>();
    
    return logdet_A;

}




char* mpinv::get_file_in() {
    return file_in_;
}

void mpinv::set_file_in(const char* val) {
	if (file_in_ != nullptr) {
		delete[] file_in_;
	}

	auto count = strlen(val) + 1;
	file_in_ = new char[count];
	strcpy(file_in_, val);
}

char* mpinv::get_file_out() {
    return file_out_;
}

void mpinv::set_file_out(const char* val) {
	if (file_out_ != nullptr) {
		delete[] file_out_;
	}

	auto count = strlen(val) + 1;
	file_out_ = new char[count];
	strcpy(file_out_, val);
}

    

