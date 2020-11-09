print("""\
##################
### TEST mpinv ###  
##################\
""")

import mpmath as mp; mp.mp.dps = 50;
m = mp.randmatrix(10);
m = m.T*m;

import mpinv
print("\ttest mpinv.fast_mp_matrix_inverse_symm(m): ", end="")
m_inv = mpinv.fast_mp_matrix_inverse_symm(m)
print("m_inv calculated", end="")
print("\r  pass!");

print("\ttest mpinv.fast_mp_matrix_logdet_symm(m): ", end="");
m_ld = mpinv.fast_mp_matrix_logdet_symm(m)
m_inv_ld = mpinv.fast_mp_matrix_logdet_symm(m_inv)
print("logdet m=", m_ld, ", m_inv=",m_inv_ld,sep="",end="")
print("\r  pass!");


print("\ttest mpinv.fast_mp_matrix_inverse_with_logdet_symm(m): ", end="")
m_ld, m_inv = mpinv.fast_mp_matrix_inverse_with_logdet_symm(m)
print("tuple (logdet,m_inv) calculated", end="")
print("\r  pass!");

print()