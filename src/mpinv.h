#ifndef MPINV_H
#define MPINV_H


class mpmat {
private:
	void* ptr_mpmat = nullptr;
	void* ptr_mpmat_inv = nullptr;
	double logdet = 0.0;
	
	//static const int MAX_POSSIBLE_PRECISION = 100;
    //char _curr_element_chars[MAX_POSSIBLE_PRECISION] = "";
	char _curr_coeff_chars[100] = "";

	void calc_logdet1(bool is_symmetric = false);
	void calc_logdet2(bool is_symmetric = false);
	void calc_logdet3(bool is_symmetric = false);

	void clear_ptr_mpmat(void* ptr_mpmat);

public:
	const char* version_ = "0.1.1";
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
	void save_mpmat(const char* file_out="M.mpmat");
	void save_mpmat_inv(const char* file_out="M_inv.mpmat");

	void calc_invert();
	void calc_logdet(bool is_symmetric = false);
	void calc_invert_with_logdet(bool is_symmetric = false);


	double get_logdet();
	
	void get_data(const char* data);


	//mpmat& operator=(mpmat& m);

};


#endif // MPINV_H


