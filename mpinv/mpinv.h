#ifndef FOO_LIB
#define FOO_LIB

class mpinv {
private:
    char* file_in_;
    char* file_out_;
    
public:
    mpinv(char* file_in, char* file_out);
    virtual ~mpinv();
    
    void inv();
    
    double logdet(bool is_symmetric=false);
    
    char* get_file_in();
    void set_file_in(const char* val);

    char* get_file_out();
    void set_file_out(const char* val);


};

#endif
