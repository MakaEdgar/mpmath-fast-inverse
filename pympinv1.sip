%Module(name=mpinv1, language="C++")
%DefaultEncoding "UTF-8"

class mpinv {
    %TypeHeaderCode
    #include <mpinv.h>
    %End

    public:
        mpinv(char*, char*);
        
        void inv();
        
        double logdet(bool);
        
        char* get_file_in();
        void set_file_in(const char*);

        char* get_file_out();
        void set_file_out(const char*);


};