import mpmath as mp

def slow_mp_matrix_inverse(m):
    return m**-1

def slow_mp_matrix_logdet(m):
    return mp.log(mp.det(m))

def slow_mp_matrix_inverse_with_logdet(m):
    logdet = slow_mp_matrix_logdet(m)
    m_inv = slow_mp_matrix_inverse(m)
    return (logdet, m_inv)


def mp_matrix_to_mpmatfile(m, outfile="matrix.mpmat"):
    with open(outfile, "w") as f:
        f.write(str(m.rows) + " " + str(m.cols) + "\n")
        for row in m.tolist():
            for m_ij in row:
                f.write(str(m_ij) + " ")
            f.write("\n")

def mp_matrix_from_mpmatfile(infile="matrix.mpmat"):
    with open(infile, "r") as f:
        file = [x.strip().split() for x in f.read().split("\n") if x.strip() != ""]
    matrix_list = [list(map(mp.mpf, x)) for x in file[1:]]
    return mp.matrix(matrix_list)

