import mpmath as mp

def slow_mp_matrix_inverse(m):
    return m**-1

def slow_mp_matrix_logdet(m):
    return mp.log(mp.det(m))

def slow_mp_matrix_inverse_with_logdet(m):
    logdet = slow_mp_matrix_logdet(m)
    m_inv = slow_mp_matrix_inverse(m)
    return (logdet, m_inv)



