%Module(name=mpinv, language="C++")
%DefaultEncoding "UTF-8"

class mpmat {
    %TypeHeaderCode
    #include <mpinv.h>
    %End

    public:
		const char* _VERSION;
		int _MP_PRECISION;
		int _MAX_MP_PRECISION;

        int rows;
		int cols;

		mpmat();
		mpmat(const char*);
		mpmat(int rows_, int cols_);

		void set_matrix_coeff(int i, int j, const char* doubleMP_str);
		const char* get_matrix_coeff(int i, int j);
		const char* get_matrix_inv_coeff(int i, int j);

        void load_mpmat(const char*);
		void save_mpmat(const char*);
        void save_mpmat_inv(const char*);


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
		double get_det(bool to_calculate=true);


};
