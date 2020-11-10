#ifndef MPINV_H
#define MPINV_H

const int MP_PRECISION = 50;
const int MAX_MP_PRECISION = 100;

class mpmat {
private:
	void* ptr_mpmat = nullptr;
	void* ptr_mpmat_inv = nullptr;
	double logdet = 0.0;
	double det = 0.0;
		
	void _clear_ptr_mpmat(void* ptr_mpmat);

	char _curr_coeff_chars[MAX_MP_PRECISION + 1] = "";
	void _get_matrix_coeff(int i, int j, void* ptr_mpmat, char* coeff_chars);

	void _save_mpmat(const char* file_out, void* ptr_mpmat);

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

	double get_logdet(bool to_calculate=true);
	double get_det(bool to_calculate = true);

};

#endif // MPINV_H


