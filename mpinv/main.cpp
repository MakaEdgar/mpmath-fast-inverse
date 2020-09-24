
#include "mpinv.h"
#include <iostream>

int main(int, char**)
{
    mpinv foo("in.mpmat", "in_inverse.mpmat");
    foo.inv();
    std::cout << "Logdet is " << foo.logdet() << std::endl;
    std::cout << foo.get_file_in() << std::endl;
    std::cout << foo.get_file_out() << std::endl;
    
    foo.set_file_in("bazaza");
    foo.set_file_out("bazaza_out");
    
    std::cout << foo.get_file_in() << std::endl;
    std::cout << foo.get_file_out() << std::endl;

    return 0;
}



