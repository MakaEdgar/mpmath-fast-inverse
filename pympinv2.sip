%Module(name=mpinv2, language="C++")
%DefaultEncoding "UTF-8"

class mpmat {
    %TypeHeaderCode
    #include <mpinv.h>
    %End

    public:
        int cols;
		int rows;

		mpmat();
		mpmat(const char*);
        
        void load_mpmat(const char*);
		void save_mpmat(const char*);
        void save_mpmat_inv(const char*);
        
		void calc_invert();
		void calc_logdet(bool);
		void calc_invert_with_logdet(bool);
			   		 
        double get_logdet();

		void get_data(const char*);

};
