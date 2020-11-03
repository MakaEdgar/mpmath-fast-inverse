#ifndef MPINV_H
#define MPINV_H


class mpmat {
//private:
	void* ptr_mpmat = nullptr;
	void* ptr_mpmat_inv = nullptr;
	double logdet = 0.0;

	void calc_logdet1(bool is_symmetric = false);
	void calc_logdet2(bool is_symmetric = false);
	void calc_logdet3(bool is_symmetric = false);

	void clear_ptr_mpmat(void* ptr_mpmat);


public:
	const char* version_ = "0.1.0";
	int cols = 0;
	int rows = 0;

	mpmat();
	mpmat(const char* file_in);

	virtual ~mpmat();
	
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


