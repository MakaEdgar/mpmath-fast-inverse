import mpmath as mp
from .mpinv import mpmat as mpmat_cpp


class mpmat:
    _mat = None
    
    def __init__(self, m):
        self._mat = mpmat_cpp(m.rows, m.cols)
        for i in range(m.rows):
            for j in range(m.cols):
                self._mat.set_matrix_coeff(i,j,str(m[i,j])) 
                
    def get_mat(self):
        return mp.matrix([[self._mat.get_matrix_coeff(i,j) for j in range(m.cols)] for i in range(m.rows)])
    
    def set_mat(self, m):
        self.__init__(m)
        
    def get_inv(self, to_calculate=True, method="base"):
        if to_calculate:
            if method is "base":
                self._mat.calc_inverse()
            elif method is "LU":
                self._mat.calc_inverse_LU()
            elif method is "LLT":
                self._mat.calc_inverse_LLT()
            elif method is "LDLT":
                self._mat.calc_inverse_LDLT()
        return mp.matrix([[self._mat.get_matrix_inv_coeff(i,j) for j in range(m.cols)] for i in range(m.rows)])
    
    def get_det(self, to_calculate=True, method="base", return_float=False):
        if to_calculate:
            if method is "base":
                self._mat.calc_det()
            elif method is "LU":
                self._mat.calc_det_LU()
            elif method is "LLT":
                self._mat.calc_det_LLT()
            elif method is "LDLT":
                self._mat.calc_det_LDLT()
        if return_float:   
            return self._mat.get_det()
        return mp.mpf(self._mat.get_det_str())
     
    def get_logdet(self, to_calculate=True, method="base"):
        if to_calculate:
            if method is "base":
                self._mat.calc_logdet()
            elif method is "LU":
                self._mat.calc_logdet_LU()
            elif method is "LLT":
                self._mat.calc_logdet_LLT()
            elif method is "LDLT":
                self._mat.calc_logdet_LDLT()        
        return self._mat.get_logdet() 
    
    def get_logdet_with_inverse(self, to_calculate=True, method="base"):
        if to_calculate:
            if method is "base":
                self._mat.calc_inverse_with_logdet()
            elif method is "LU":
                self._mat.calc_inverse_with_logdet_LU()
            elif method is "LLT":
                self._mat.calc_inverse_with_logdet_LLT()
            elif method is "LDLT":
                self._mat.calc_inverse_with_logdet_LDLT()        
        return (self.get_logdet(False), self.get_inv(False))
    
    
    
def fast_mp_matrix_inverse(m):
    return mpmat(m).get_inv()

def fast_mp_matrix_det(m):
    return mpmat(m).get_det()

def fast_mp_matrix_logdet(m):
    return mpmat(m).get_logdet()

def fast_mp_matrix_inverse_with_logdet(m):
    return mpmat(m).get_logdet_with_inverse()


def fast_mp_matrix_inverse_symm(m):
    return mpmat(m).get_inv(method="LLT")

def fast_mp_matrix_det_symm(m):
    return mpmat(m).get_logdet(method="LLT")

def fast_mp_matrix_logdet_symm(m):
    return mpmat(m).get_logdet(method="LLT")

def fast_mp_matrix_inverse_with_logdet_symm(m):
    return mpmat(m).get_logdet_with_inverse(method="LLT")

