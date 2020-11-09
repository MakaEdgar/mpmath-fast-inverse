mpinv: python library to calculate mpmath.matrix invert and logdet.
(c) Edgar Makarov, 2020 


How to use:
1. clone git repo
2. execute `./run_build.sh`
3. after successful building `mpinv*.tar.gz and .whl` files will be created in bin/py_module/
4. execute `pip3 install --user mpinv*.tar.gz(or .whl)`
5. Enjoy mpinv library by using functions:
```
    import mpinv
    import mpmath as mp; mp.mp.dps = 50
    
    m_inv = mpinv.fast_mp_matrix_inverse_symm(m)                                  # returns mp.matrix
    m_logdet = mpinv.fast_mp_matrix_logdet_symm(m)                                # returns float64
    m_logdet, m_inv = mpinv.fast_mp_matrix_inverse_with_logdet_symm(m)            # returns tuple (float64, mp.matrix)

    where m is symmetric (and SPD) mp.matrix
```


Other functions and capabilities are under development.
In case of any bugs send the log and other info to e.makarov@skoltech.ru
