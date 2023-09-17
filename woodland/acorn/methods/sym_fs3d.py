from sympy import *
import re

# Column-major R. Opposite of what we decided to use in util.hpp, so RARt is
# actually RtAR in util.hpp.
def rotate_sym_tensor_3x3_RARt():
    var('t11 t12 t13 t22 t23 t33 x1 x2 x3 y1 y2 y3 z1 z2 z3')
    R = Matrix([[x1, y1, z1], [x2, y2, z2], [x3, y3, z3]])
    T = Matrix([[t11, t12, t13], [t12, t22, t23], [t13, t23, t33]])
    assert(Eq(T - T.T, Matrix.zeros(3)))
    S = R*T*R.T
    pprint(R)
    pprint(T)
    for i in range(3):
        for j in range(3):
            if j < i: continue
            print(simplify(S[i,j]))

def rotate_sym_tensor_3x3_RtAR():
    var('t11 t12 t13 t22 t23 t33 x1 x2 x3 y1 y2 y3 z1 z2 z3')
    R = Matrix([[x1, y1, z1], [x2, y2, z2], [x3, y3, z3]])
    T = Matrix([[t11, t12, t13], [t12, t22, t23], [t13, t23, t33]])
    assert(Eq(T - T.T, Matrix.zeros(3)))
    S = R.T*T*R
    pprint(R)
    pprint(T)
    for i in range(3):
        for j in range(3):
            if j < i: continue
            print(simplify(S[i,j]))

# refs:
# [GGBT12] doi: 10.1785/0120120127

# GGBT12 eqs 1b, 2b. W.r.t. GGBT12, y = xi, lam = lambda, a = alpha.
# U[i,j] = U_i^j
# T[i,j,k] = T_{ik}^j
def make_U_and_T_K(x, y):
    var('mu lam', real=True)
    a = (lam + mu)/(lam + 2*mu)
    R = 0
    for i in range(3):
        R += (x[i] - y[i])**2
    R = sqrt(R)
    U = Matrix.zeros(3)
    for i in range(3):
        for j in range(3):
            # N.B. this disagrees in two ways with the eq before GGBT12 eq 20:
            #   1. sign
            #   2. compensating for the 1/R term when i = j.
            # Segall ch. 3 eq 3.90 -> 3.91 handles these details correctly.
            U[i,j] = -a*R.diff(y[j]).diff(y[i])
            if i == j:
                U[i,j] += 2/R
            U[i,j] *= 1/(8*pi*mu)
            U[i,j] = simplify(U[i,j])
    T = array.MutableDenseNDimArray.zeros(3,3,3)
    for i in range(3):
        Uidiv = 0
        for n in range(3):
            Uidiv += U[i,n].diff(y[n])
        for j in range(3):
            for k in range(3):
                T[i,j,k] = mu*(U[i,j].diff(y[k]) + U[i,k].diff(y[j]))
                if j == k:
                    T[i,j,k] += lam*Uidiv
    return (mu, lam, R, U, T)

def calc_Tvd(T, v, d):
    u = Matrix.zeros(3,1)
    for i in range(3):
        for k in range(3):
            for j in range(3):
                u[i] += T[i,j,k]*v[k]*d[j]
    return u

def calc_sigma(mu, lam, x, u):
    sigma = Matrix.zeros(3)
    udiv = 0
    for i in range(3):
        udiv += u[i].diff(x[i])
    for i in range(3):
        for j in range(3):
            sigma[i,j] = mu*(u[i].diff(x[j]) + u[j].diff(x[i]))
            if i == j:
                sigma[i,j] += lam*udiv
    return sigma

def sub_dict(d, e, verbose=False):
    if verbose: print('e:', e)
    for k, v in d.items():
        if e == k:
            if verbose: print('  found:', v)
            return v
    if len(e.args) == 0:
        return e
    new_args = []
    for i in range(len(e.args)):
        new_args.append(sub_dict(d, e.args[i], verbose))
    return e.func(*new_args)

def write_c(d):
    def sub(c):
        c = re.sub(r'pow\(([^\(,]+),\s*2\)',
                   r'square(\1)', c)
        c = re.sub(r'pow\(([^,]+),\s*2\)',
                   r'square(\1)', c)
        c = re.sub(r'pow\(([^,]+),\s*3\)',
                   r'cube(\1)', c)
        # I can't get the math sub to work, so sub it at the code level.
        c = c.replace('pow(R2, 7.0/2.0)', 'R7')
        # More text substitutions.
        for (s, o) in [('z11', 'square(z1)'),
                       ('z111', 'cube(z1)'),
                       ('z12', 'z1*z2'),
                       ('z13', 'z1*z3'),
                       ('z22', 'square(z2)'),
                       ('z222', 'cube(z2)'),
                       ('z23', 'z2*z3'),
                       ('z33', 'square(z3)'),
                       ('z333', 'cube(z3)'),
                       ('z11pz22', '(z11 + z22)'),
                       ('z11pz33', '(z11 + z33)'),
                       ('z22pz33', '(z22 + z33)'),
                       ('R4', 'square(R2)'),
                       ('den1', '(M_PI*R7*l2m*mu)'),
                       ('den2', '(M_PI*R7*l2m)')]:
            c = c.replace(o, s)
        return c
    for k, v in d.items():
        if type(v) == Symbol:
            lhs = ccode(v)
        elif type(v) == str:
            lhs = v
        else:
            continue
        print('{} = {};'.format(lhs, sub(ccode(k))), flush=True)

def cpp_sigma():
    # Monolithic code. Giant one-liner for each sigma.
    var('x1 x2 x3 y1 y2 y3 v1 v2 v3 d1 d2 d3 ' +
        'R2 z1 z2 z3 lm l2m l3m ' +
        'sigma11 sigma12 sigma13 sigma22 sigma23 sigma33', real=True)
    x = Matrix([x1, x2, x3])
    y = Matrix([y1, y2, y3])
    v = Matrix([v1, v2, v3])
    d = Matrix([d1, d2, d3])
    sigma_sym = Matrix([[sigma11, sigma12, sigma13],
                        [sigma12, sigma22, sigma23],
                        [sigma13, sigma23, sigma33]])
    mu, lam, R, U, T = make_U_and_T_K(x, y)
    u = calc_Tvd(T, v, d)
    sigma = calc_sigma(mu, lam, x, u)
    print(U[1,1])
    print(U[1,2])
    print(T[1,2,0])
    subs = {R**2: R2,
            x1-y1: z1, x2-y2: z2, x3-y3: z3,
            -x1+y1: -z1, -x2+y2: -z2, -x3+y3: -z3,
            z1+x1-y1: 2*z1, z2+x2-y2: 2*z2, z3+x3-y3: 2*z3,
            z1+2*x1-2*y1: 3*z1, z2+2*x2-2*y2: 3*z2, z3+2*x3-2*y3: 3*z3,
            lam+mu: lm, lam+2*mu: l2m, lam+3*mu: l3m}
    write_c(subs)
    for i in range(3):
        for j in range(3):
            if j < i: continue
            s = sigma[i,j]
            s = sub_dict(subs, s)
            s = sub_dict(subs, simplify(s))
            write_c({s: sigma_sym[i,j]})

def cpp_sigma_alt():
    # u(i) = T(i,j,k) v(k) d(j) is linear in T, so take derivs of T rather than
    # u to produce 3^4 = 81 separate T lines to compute.
    #   With simplification, this turns out to be better than the monolithic
    # code:
    #   $ wc fs3d_sigma*cpp
    #    34  3274 29577 fs3d_sigma_monolithic.cpp
    #   117  2030 19709 fs3d_sigma_Tx81.cpp
    # There's way more line-level parallelism and possibly fewer calcs if
    # character count is a reasonable proxy.
    var('x1 x2 x3 y1 y2 y3 v1 v2 v3 d1 d2 d3 ' +
        'R2 R7 z1 z2 z3 lm l2m l3m ' +
        'sigma11 sigma12 sigma13 sigma22 sigma23 sigma33', real=True)
    x = Matrix([x1, x2, x3])
    y = Matrix([y1, y2, y3])
    v = Matrix([v1, v2, v3])
    d = Matrix([d1, d2, d3])
    mu, lam, R, U, T = make_U_and_T_K(x, y)
    subs = {R**2: R2,
            #R2**(Integer(7)/2): R7,  # not getting this to sub successfully
            x1-y1: z1, x2-y2: z2, x3-y3: z3,
            -x1+y1: -z1, -x2+y2: -z2, -x3+y3: -z3,
            z1+x1-y1: 2*z1, z2+x2-y2: 2*z2, z3+x3-y3: 2*z3,
            z1+2*x1-2*y1: 3*z1, z2+2*x2-2*y2: 3*z2, z3+2*x3-2*y3: 3*z3,
            lam+mu: lm, lam+2*mu: l2m, lam+3*mu: l3m}
    #write_c(subs)
    print("""  const Real
    lm = lam + mu,
    l2m = lam + 2*mu,
    l3m = lam + 3*mu;
  const RealT
    x1 = x[0], x2 = x[1], x3 = x[2],
    y1 = y[0], y2 = y[1], y3 = y[2],
    z1 = x1 - y1,
    z2 = x2 - y2,
    z3 = x3 - y3,
    z11 = square(z1),
    z111 = z1*z11,
    z12 = z1*z2,
    z13 = z1*z3,
    z22 = square(z2),
    z222 = z2*z22,
    z23 = z2*z3,
    z33 = square(z3),
    z333 = z3*z33,
    z11pz22 = z11 + z22,
    z11pz33 = z11 + z33,
    z22pz33 = z22 + z33,
    R2 = z11 + z22 + z33,
    R4 = square(R2),
    R7 = acorn::pow(R2, 7.0/2.0),
    den1 = M_PI*R7*l2m*mu,
    den2 = M_PI*R7*l2m;
  RealT T_x[81];""")
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for d in range(3):
                    T_x = T[i,j,k].diff(x[d])
                    T_x = sub_dict(subs, T_x)
                    T_x = sub_dict(subs, simplify(T_x))
                    write_c({T_x: 'T_x[{}]'.format(3*(3*(3*i+j)+k)+d)})
    print('Real u_x[9];')
    for i in range(3):
        for d in range(3):
            s = ''
            for j in range(3):
                for k in range(3):
                    s += (' T_x[{}]*d[{}]*v[{}]'.
                          format(3*(3*(3*i+j)+k)+d, j, k))
                    if j == 2 and k == 2:
                        s += ';'
                    else:
                        s += ' +'
            print('u_x[{}] ={}'.format(3*i+d, s))
    print('const Real udiv = u_x[0] + u_x[4] + u_x[8];')
    k = 0
    for i in range(3):
        for j in range(3):
            if j < i: continue
            s = ''
            if i == j: s += 'lam*udiv + ';
            s += 'mu*(u_x[{}] + u_x[{}])'.format(3*i+j, 3*j+i)
            print('sigma[{}] = {};'.format(k, s))
            k = k + 1

def cmp_fsps():
    # Compare with sym_fsps.c1().
    var('x1 x2 x3 y1 y2 y3 v1 v2 v3 d1 d2 d3 ' +
        'R2 z1 z2 z3 lm l2m l3m ' +
        'sigma11 sigma12 sigma13 sigma22 sigma23 sigma33', real=True)
    x = Matrix([x1, x2, x3])
    y = Matrix([y1, y2, y3])
    v = Matrix([0, 1, 0])
    d = Matrix([1, 0, 0])
    mu, lam, R, U, T = make_U_and_T_K(x, y)
    u = calc_Tvd(T, v, d)
    sigma = calc_sigma(mu, lam, x, u)
    f = (sigma.subs(y1, 0).subs(y2, 0)
         .subs(x3, 0) # not needed by symmetry, but make the problem easier
         .subs(lam, 1).subs(mu, 1))
    # infinite out-of-plane dislocation
    sigma = simplify(integrate(f, (y3, -oo, oo)))
    # point dislocation for a quick look at signs
    #sigma = simplify(f.subs(y3, 0))
    print(sigma[0,0])
    print(sigma[0,1])
    print(sigma[1,1])

cpp_sigma_alt()
