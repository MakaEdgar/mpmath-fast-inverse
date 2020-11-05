import mpmath as mp; mp.mp.dps = 50
from mpinv2 import mpmat



def slow_mp_matrix_inverse(m):
    return m**-1

def slow_mp_matrix_logdet(m):
    return mp.log(mp.det(m))

def slow_mp_matrix_inverse_with_logdet(m):
    logdet = slow_mp_matrix_logdet(m)
    m_inv = slow_mp_matrix_inverse(m)
    return (logdet, m_inv)





def fast_mp_matrix_inverse(m):
    a = mpmat(m.rows, m.cols)
    for i in range(m.rows):
        for j in range(m.cols):
            a.set_matrix_coeff(i,j,str(m[i,j])) 
    a.calc_invert()
    m_inv = mp.matrix([[a.get_matrix_inv_coeff(i,j) for j in range(m.cols)] for i in range(m.rows)])
    return m_inv


def fast_mp_matrix_logdet(m, is_symmetric=False):
    a = mpmat(m.rows, m.cols)
    for i in range(m.rows):
        for j in range(m.cols):
            a.set_matrix_coeff(i,j,str(m[i,j])) 
    a.calc_logdet(is_symmetric)
    logdet = a.get_logdet()
    return logdet

def fast_mp_matrix_inverse_with_logdet(m, is_symmetric=False):
    a = mpmat(m.rows, m.cols)
    for i in range(m.rows):
        for j in range(m.cols):
            a.set_matrix_coeff(i,j,str(m[i,j])) 
    a.calc_invert_with_logdet(is_symmetric)
    logdet = a.get_logdet()
    m_inv = mp.matrix([[a.get_matrix_inv_coeff(i,j) for j in range(m.cols)] for i in range(m.rows)])
    return(logdet, m_inv)


def fast_mp_matrix_inverse_with_logdet_symm(m):
    return fast_mp_matrix_inverse_with_logdet(m,True)

def fast_mp_matrix_inverse_symm(m):
    return fast_mp_matrix_inverse(m)

def fast_mp_matrix_logdet_symm(m):
    return fast_mp_matrix_logdet(m,True)

