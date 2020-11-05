mpinv2: python library to calculate mp.matrix invert and logdet.
(c) Edgar Makarov, 2020 


How to use:
1. clone git repo
2. execute ./run_ruild.sh
3. after successful building py_mpinv2*.whl file will be created in bin/py_module/, together with python wrapper mpinv.py
4. execute "pip3 install --user py_mpinv2-blabla.whl"
5. import mpinv.py and use functions:
    m_inv = fast_mp_matrix_inverse_symm(m)                                  # returns mpmath.matrix
    m_logdet = fast_mp_matrix_logdet_symm(m)                                # returns float64
    m_logdet, m_inv = fast_mp_matrix_inverse_with_logdet_symm(m)            # returns tuple (float64, mpmath.matrix)

    where m is symmetric (and SPD) mpmath.matrix

Other functions and capabilities are under development.
In case of any bugs send the log and other info to e.makarov@skoltech.ru
