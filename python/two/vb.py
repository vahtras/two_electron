from util.full import matrix
from .core import list_integrals

def vb_transform(dens, delta, **kwargs):
    filename = kwargs.get('filename', '/tmp/AOTWOINT')
    a, m = dens[0].shape
    H_uvmn = matrix((a, a, m, m))
    for ig, g in list_integrals(filename):
        p, q, r, s = ig
        s, r, q, p = p-1, q-1, r-1, s-1
        if p == q: g *= 0.5
        if r == s: g *= 0.5
        if (p, q) == (r, s): g *= 0.5

        for d1, D1 in zip(dens, delta):
            for d2, D2 in zip(dens, delta):
                H_uvmn += D1[:, q].x(D2[:, s].x(d1[p, :].x(d2[r, :]*g)))
                H_uvmn += D1[:, p].x(D2[:, s].x(d1[q, :].x(d2[r, :]*g)))
                H_uvmn += D1[:, q].x(D2[:, r].x(d1[p, :].x(d2[s, :]*g)))
                H_uvmn += D1[:, p].x(D2[:, r].x(d1[q, :].x(d2[s, :]*g)))
                H_uvmn += D1[:, s].x(D2[:, q].x(d1[r, :].x(d2[p, :]*g)))
                H_uvmn += D1[:, r].x(D2[:, q].x(d1[s, :].x(d2[p, :]*g)))
                H_uvmn += D1[:, s].x(D2[:, p].x(d1[r, :].x(d2[q, :]*g)))
                H_uvmn += D1[:, r].x(D2[:, p].x(d1[s, :].x(d2[q, :]*g)))

            H_uvmn -= D1[:, s].x(D1[:, q].x(d1[p, :].x(d1[r, :]*g)))
            H_uvmn -= D1[:, s].x(D1[:, p].x(d1[q, :].x(d1[r, :]*g)))
            H_uvmn -= D1[:, r].x(D1[:, q].x(d1[p, :].x(d1[s, :]*g)))
            H_uvmn -= D1[:, r].x(D1[:, p].x(d1[q, :].x(d1[s, :]*g)))
            H_uvmn -= D1[:, q].x(D1[:, s].x(d1[r, :].x(d1[p, :]*g)))
            H_uvmn -= D1[:, q].x(D1[:, r].x(d1[s, :].x(d1[p, :]*g)))
            H_uvmn -= D1[:, p].x(D1[:, s].x(d1[r, :].x(d1[q, :]*g)))
            H_uvmn -= D1[:, p].x(D1[:, r].x(d1[s, :].x(d1[q, :]*g)))
        
    return H_uvmn.transpose(0, 2, 1, 3)

def vb_transform2(Dma, Dam, Delta1, Delta2, **kwargs):
    filename = kwargs.get('filename', '/tmp/AOTWOINT')
    H_umvn = matrix(Dam[0].shape + Dam[0].shape)
    for ig, g in list_integrals(filename):
        p, q, r, s = ig
        s, r, q, p = p-1, q-1, r-1, s-1
        if p == q: g *= 0.5
        if r == s: g *= 0.5
        if (p, q) == (r, s): g *= 0.5

        for d1a, d2a, D1a, D2a in zip(Dma, Dam, Delta1, Delta2):
            for d1b, d2b, D1b, D2b in zip(Dma, Dam, Delta1, Delta2):
               H_umvn += D1a[p,:].x(d1a[:, q].x(D2b[:, s].x(d2b[r, :]*g)))
               H_umvn += D1a[q,:].x(d1a[:, p].x(D2b[:, s].x(d2b[r, :]*g)))
               H_umvn += D1a[p,:].x(d1a[:, q].x(D2b[:, r].x(d2b[s, :]*g)))
               H_umvn += D1a[q,:].x(d1a[:, p].x(D2b[:, r].x(d2b[s, :]*g)))

               H_umvn += D1a[r,:].x(d1a[:, s].x(D2b[:, q].x(d2b[p, :]*g)))
               H_umvn += D1a[r,:].x(d1a[:, s].x(D2b[:, p].x(d2b[q, :]*g)))
               H_umvn += D1a[s,:].x(d1a[:, r].x(D2b[:, q].x(d2b[p, :]*g)))
               H_umvn += D1a[s,:].x(d1a[:, r].x(D2b[:, p].x(d2b[q, :]*g)))

            H_umvn -= D1a[p, :].x(d1a[:, s].x(D2a[:, q].x(d2a[r, :]*g)))
            H_umvn -= D1a[q, :].x(d1a[:, s].x(D2a[:, p].x(d2a[r, :]*g)))
            H_umvn -= D1a[p, :].x(d1a[:, r].x(D2a[:, q].x(d2a[s, :]*g)))
            H_umvn -= D1a[q, :].x(d1a[:, r].x(D2a[:, p].x(d2a[s, :]*g)))

            H_umvn -= D1a[r, :].x(d1a[:, q].x(D2a[:, s].x(d2a[p, :]*g)))
            H_umvn -= D1a[r, :].x(d1a[:, p].x(D2a[:, s].x(d2a[q, :]*g)))
            H_umvn -= D1a[s, :].x(d1a[:, q].x(D2a[:, r].x(d2a[p, :]*g)))
            H_umvn -= D1a[s, :].x(d1a[:, p].x(D2a[:, r].x(d2a[q, :]*g)))

    return H_umvn
