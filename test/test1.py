print("""\
##################
### TEST mpinv ###   !!! current version works only with 50 mpmath digits
##################   !!! include \"mp.mp.dps = 50\" in your code""")


print("test mpinv(base)");
print("\ttest import: ", end="")
import mpinv2
a=mpinv2.mpmat()
print("lib version=", a.version_, end="")
print("\r  pass!")



print("test mpinv.py")
import sys;sys.path.append("bin/py_module")
print("\ttest import",end="\t")
import mpinv
print("\r  pass!");



import mpmath as mp; mp.mp.dps = 50;
m = mp.randmatrix(10);
m = m.T*m;

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