%Module(name=mpinv2, language="C++")
%DefaultEncoding "UTF-8"

class mpmat {
    %TypeHeaderCode
    #include <mpinv.h>
    %End

    public:
		const char* version_;
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
        
		void calc_invert();
		void calc_logdet(bool);
		void calc_invert_with_logdet(bool);
			   		 
        double get_logdet();

		void get_data(const char*);

};
