#ifndef MPINV_H
#define MPINV_H

const int MP_PRECISION = 50;
const int MAX_MP_PRECISION = 100;

class mpmat {
private:
	void* ptr_mpmat = nullptr;
	void* ptr_mpmat_inv = nullptr;

	double det = 0.0;
	double logdet = 0.0;

	char _det_str[MAX_MP_PRECISION + 1] = "";
	char _curr_coeff_str[MAX_MP_PRECISION + 1] = "";

public:
	const char* _VERSION = "0.1.2";
	int _MP_PRECISION = MP_PRECISION;
	int _MAX_MP_PRECISION = MAX_MP_PRECISION;

	int rows = 0;
	int cols = 0;

	mpmat();
	mpmat(const char* file_in);
	mpmat(int rows_, int cols_);

	virtual ~mpmat();

	void set_matrix_coeff(int i, int j, const char* doubleMP_str);
	const char* get_matrix_coeff(int i, int j);
	const char* get_matrix_inv_coeff(int i, int j);

	void load_mpmat(const char* file_in);
	void save_mpmat(const char* file_out="m.mpmat");
	void save_mpmat_inv(const char* file_out="m_inv.mpmat");

	
	void calc_logdet_LU();
	void calc_logdet_LLT();
	void calc_logdet_LDLT();
	void calc_logdet();

	void calc_det_LU();
	void calc_det_LLT();
	void calc_det_LDLT();
	void calc_det();

	void calc_inverse_LU();
	void calc_inverse_LLT();
	void calc_inverse_LDLT();
	void calc_inverse();

	void calc_inverse_with_logdet_LU();
	void calc_inverse_with_logdet_LLT();
	void calc_inverse_with_logdet_LDLT();
	void calc_inverse_with_logdet();

	double get_logdet();
	double get_det();
	const char* get_det_str();

};

#endif // MPINV_H


