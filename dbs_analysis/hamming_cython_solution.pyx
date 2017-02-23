def hamming_loop(char *s0, char *s1):
    if len(s0) != len(s1):
        raise ValueError()
    cdef:
        int N = len(s0)
        int i, count = 0
    for i in range(N):
        count += (s0[i] != s1[i])
    return count

# adapted from https://github.com/kwmsmith/scipy-2015-cython-tutorial
