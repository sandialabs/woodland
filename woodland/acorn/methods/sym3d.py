from sympy import *
import sympy, re

# refs:
# [GGBT12] doi: 10.1785/0120120127
# Okada, 1992
# Segall book, 2010
# GGBT12 eqs 1b, 2b. W.r.t. GGBT12, y = xi, lam = lambda, a = alpha.
# U[i,j] = U_i^j
# T[i,j,k] = T_{ik}^j

def make_T(x, y, U):
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
    return T

def calc_sigma(lam, mu, x, u):
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

def lam_mu_to_alpha(lam, mu):
    return (lam + mu)/(lam + 2*mu)

# Fullspace U and derivatives.
def calc_T_K(x, y):
    var('lam mu', real=True)
    a = lam_mu_to_alpha(lam, mu)
    R = 0
    for i in range(3):
        R += (x[i] - y[i])**2
    R = sqrt(R)
    U = Matrix.zeros(3)
    for i in range(3):
        for j in range(3):
            # N.B. the following disagrees in two ways with the eq before GGBT12
            # eq 20:            
            #   1. sign
            #   2. compensating for the 1/R term when i = j.
            # It agrees with other sources, e.g., Segall ch. 3 eq 3.90 -> 3.91.
            U[i,j] = -a*R.diff(y[j]).diff(y[i])
            if i == j:
                U[i,j] += 2/R
            U[i,j] *= 1/(8*pi*mu)
    T = make_T(x, y, U)
    return lam, mu, R, T

# Halfspace terms.
def calc_halfspace_coords(x, y):
    r = [x[0] - y[0], x[1] - y[1], -(x[2] + y[2])]
    R = 0
    for i in range(3):
        R += r[i]**2
    R = sqrt(R)
    return r, R

def calc_T_A(lam, mu, x, y):
    a = lam_mu_to_alpha(lam, mu)
    r, R = calc_halfspace_coords(x, y)
    U = Matrix.zeros(3)
    for i in range(3):
        for j in range(3):
            U[i,j] = a*r[i]*r[j]/R**3
            if i == j:
                U[i,j] += (2-a)/R
            U[i,j] *= -1/(8*pi*mu)
    T = make_T(x, y, U)
    return R, T

def calc_T_B(lam, mu, x, y):
    a = lam_mu_to_alpha(lam, mu)
    r, R = calc_halfspace_coords(x, y)
    f = (1 - a)/a
    rr3 = R + r[2]
    U = Matrix.zeros(3)
    for i in range(3):
        for j in range(3):
            U[i,j] = 0
            if i == j:
                U[i,j] += 1/R
            U[i,j] += r[i]*r[j]/R**3
            if i == j:
                U[i,j] += f/rr3
            if j == 2:
                U[i,j] += f*r[i]/(R*rr3)
            if i == 2 and j != 2:
                U[i,j] += -f*r[j]/(R*rr3)
            if i != 2 and j != 2:
                U[i,j] += -f*r[i]*r[j]/(R*rr3**2)
            U[i,j] /= 4*pi*mu
    T = make_T(x, y, U)
    return R, T

def calc_T_C(lam, mu, x, y):
    a = lam_mu_to_alpha(lam, mu)
    r, R = calc_halfspace_coords(x, y)
    U = Matrix.zeros(3)
    for i in range(3):
        for j in range(3):
            dij = 1 if i == j else 0
            di  = 1 if i == 2 else 0
            dj  = 1 if j == 2 else 0
            U[i,j] = (x[2]/(4*pi*mu)*(1 - 2*di)*
                      ((2 - a)*(r[i]*dj - r[j]*di)/R**3 +
                       a*y[2]*(dij/R**3 - 3*r[i]*r[j]/R**5)))
    T = make_T(x, y, U)
    return R, T

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

def make_c(d, text_subs=None):
    if text_subs is None: text_subs = []
    def sub(c):
        c = re.sub(r'pow\(([^\(,]+),\s*2\)',
                   r'square(\1)', c)
        c = re.sub(r'pow\(([^,]+),\s*2\)',
                   r'square(\1)', c)
        c = re.sub(r'pow\(([^,]+),\s*3\)',
                   r'cube(\1)', c)
        c = c.replace('pow(R2, 7.0/2.0)', 'R7')
        for (s, o) in [('z11', 'square(z1)'),
                       ('z111', 'cube(z1)'),
                       ('z12', 'z1*z2'),
                       ('z13', 'z1*z3'),
                       ('z22', 'square(z2)'),
                       ('z222', 'cube(z2)'),
                       ('z23', 'z2*z3'),
                       ('z33', 'square(z3)'),
                       ('z333', 'cube(z3)')]:
            c = c.replace(o, s)
        for (s, o) in text_subs:
            c = c.replace(o, s)
        return c
    s = ''
    for k, v in d.items():
        if type(v) == Symbol:
            lhs = ccode(v)
        elif type(v) == str:
            lhs = v
        else:
            continue
        s += '{} = {};\n'.format(lhs, sub(ccode(k)))
    return s

def write_c(d, text_subs=None):
    s = make_c(d, text_subs)
    print(s, flush=True, end='')

def make_vars():
    var('x1 x2 x3 y1 y2 y3 v1 v2 v3 d1 d2 d3 ' +
        'R2 z1 z2 z3 lm l2m l3m ' +
        'sigma11 sigma12 sigma13 sigma22 sigma23 sigma33', real=True)
    x = Matrix([x1, x2, x3])
    y = Matrix([y1, y2, y3])
    v = Matrix([v1, v2, v3])
    d = Matrix([d1, d2, d3])
    return x, y, v, d

def make_subs_fullspace(lam, mu, R):
    subs = {R**2: R2,
            #R2**(Integer(7)/2): R7,  # not getting this to sub successfully
            x1-y1: z1, x2-y2: z2, x3-y3: z3,
            -x1+y1: -z1, -x2+y2: -z2, -x3+y3: -z3,
            z1+x1-y1: 2*z1, z2+x2-y2: 2*z2, z3+x3-y3: 2*z3,
            z1+2*x1-2*y1: 3*z1, z2+2*x2-2*y2: 3*z2, z3+2*x3-2*y3: 3*z3,
            lam+mu: lm, lam+2*mu: l2m, lam+3*mu: l3m}
    return subs

def print_tmps_fullspace():
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
    den2 = M_PI*R7*l2m;""")

def print_calcs(x, T, subs, text_subs):
    print('RealT T_x[81];')
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for d in range(3):
                    T_x = T[i,j,k].diff(x[d])
                    T_x = sub_dict(subs, T_x)
                    T_x = sub_dict(subs, simplify(T_x))
                    write_c({T_x: 'T_x[{}]'.format(3*(3*(3*i+j)+k)+d)},
                            text_subs)
    print('RealT u_x[9];')
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
    print('const RealT udiv = u_x[0] + u_x[4] + u_x[8];')
    k = 0
    for i in range(3):
        for j in range(3):
            if j < i: continue
            s = ''
            if i == j: s += 'lam*udiv + ';
            s += 'mu*(u_x[{}] + u_x[{}])'.format(3*i+j, 3*j+i)
            print('sigma[{}] = {};'.format(k, s))
            k = k + 1

def print_calcs_compressed(x, nml, disloc, T, subs, text_subs, simplify=True):
    zero = set()
    previous = {}
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for d in range(3):
                    T_x = T[i,j,k].diff(x[d])
                    T_x = sub_dict(subs, T_x)
                    if simplify: T_x = sub_dict(subs, sympy.simplify(T_x))
                    if T_x == 0:
                        zero.add((i,j,k,d))
                        continue
                    if disloc[j] == 0 or nml[k] == 0:
                        continue
                    s = make_c({T_x: ''}, text_subs)
                    lhs = 'T_x_{:<2d}'.format(3*(3*(3*i+j)+k)+d)
                    if s in previous:
                        rhs = ' = ' + previous[s] + ';\n'
                    else:
                        rhs = s
                        previous[s] = lhs
                    print('const auto {}{}'.format(lhs, rhs), flush=True, end='')
    print('RealT u_x[9];')
    for i in range(3):
        for d in range(3):
            s = ''
            first = True
            for j in range(3):
                for k in range(3):
                    if ((i,j,k,d) not in zero and
                        disloc[j] != 0 and nml[k] != 0):
                        if not first: s += ' +'
                        s += (' T_x_{:<2d}*d{}*v{}'.
                              format(3*(3*(3*i+j)+k)+d, j+1, k+1))
                        first = False
                    if j == 2 and k == 2:
                        s += ';'
            print('u_x[{}] ={}'.format(3*i+d, s))
    print('const RealT udiv = u_x[0] + u_x[4] + u_x[8];')
    k = 0
    for i in range(3):
        for j in range(3):
            if j < i: continue
            s = ''
            if i == j: s += 'lam*udiv + ';
            s += 'mu*(u_x[{}] + u_x[{}])'.format(3*i+j, 3*j+i)
            print('sigma[{}] = {};'.format(k, s))
            k = k + 1

def cpp_sigma_fullspace():
    x, y, v, d = make_vars()
    lam, mu, R, T = calc_T_K(x, y)
    subs = make_subs_fullspace(lam, mu, R)
    print_tmps_fullspace()
    text_subs = [('z11pz22', '(z11 + z22)'),
                 ('z11pz33', '(z11 + z33)'),
                 ('z22pz33', '(z22 + z33)'),
                 # The next line is unsafe b/c of possible pre- or post-fixes,
                 # but it works here.
                 ('R2', 'z11 + z22 + z33'),
                 ('R7', 'pow(R2, 7.0/2.0)'),
                 ('R4', 'square(R2)'),
                 ('den1', '(M_PI*R7*l2m*mu)'),
                 ('den2', '(M_PI*R7*l2m)')]
    print_calcs(x, T, subs, text_subs)

def cpp_sigma_fullspace_canonical_zv():
    # z = (x1,0,0), n = (n1, 0, n3)
    x, y, v, d = make_vars()
    v = Matrix([v1, 0, v3])
    lam, mu, R, T = calc_T_K(x, y)
    subs = {x2: 0, x3: 0, y1: 0, y2: 0, y3: 0,
            lam+mu: lm, lam+2*mu: l2m, lam+3*mu: l3m}
    text_subs = [('R4', 'pow(x1, 4)'),
                 ('R', 'fabs(x1)'),
                 ('den1', '(M_PI*l2m*mu*R4)'),
                 ('den2', '(M_PI*l2m*R4)'),
                 ('fac1', 'R/den1'),
                 ('fac2', 'R/den2')]
    print("""  const Real
    lm = lam + mu,
    l2m = lam + 2*mu;
  const RealT
    R3 = std::abs(x1*x1*x1),
    fac1 = 1.0/(M_PI*l2m*mu*R3),
    fac2 = 1.0/(M_PI*l2m*R3);""")
    print_calcs_compressed(x, v, d, T, subs, text_subs)

def cpp_sigma_halfspace_extra_terms(term):
    x, y, v, d = make_vars()
    var('lam mu')
    text_subs = [('pow(R2,', 'pow(z11 + z22 + z33,'),
                 ('(R2)', '(z11 + z22 + z33)'),
                 ('R', 'sqrt(R2)'),
                 ('R3', 'pow(R2, 3.0/2.0)')]
    print("""  const RealT
    d1 = d[0], d2 = d[1], d3 = d[2],
    v1 = v[0], v2 = v[1], v3 = v[2],
    z1 = x[0] - y[0], z2 = x[1] - y[1], z3 = -(x[2] + y[2]),
    z11 = z1*z1, z22 = z2*z2, z33 = z3*z3,
    z12 = z1*z2, z13 = z1*z3, z23 = z2*z3,
    R2 = z11 + z22 + z33,
    R = std::sqrt(R2);""")
    if term == 'A':
        R, T = calc_T_A(lam, mu, x, y)
        text_subs += [('R4', 'square(R2)'),
                      ('R7', 'pow(R2, 7.0/2.0)'),
                      ('R5', 'pow(R2, 5.0/2.0)'),
                      ('R2', '(R2)'),
                      ('R5', 'R*R4'),
                      ('den1', '(M_PI*l2m*R7)'),
                      ('den2', '(M_PI*l2m*mu*R7)')]
        print("""  const RealT
    R4 = R2*R2, R5 = R4*R, R7 = R5*R2,
    den1 = M_PI*l2m*R7,
    den2 = M_PI*l2m*mu*R7;""")
    elif term == 'B':
        R, T = calc_T_B(lam, mu, x, y)
        text_subs += [('R4', 'square(R2)'),
                      ('R5', 'pow(R2, 5.0/2.0)'),
                      ('R2', '(z11 + z22 + z33)'),
                      ('R5', 'pow(R2, 5.0/2.0)'),
                      ('R6', 'cube(R2)'),
                      ('R2', '(R2)'),
                      ('z1111', 'pow(z1, 4)'),
                      ('z2222', 'pow(z2, 4)'),
                      ('(1.0/R)', 'pow(R2, -1.0/2.0)'),
                      ('(1.0/R3)', 'pow(R2, -3.0/2.0)'),
                      ('z3pR', '-x3 - y3 + R'),
                      ('z3pRpow4', 'pow(z3pR, 4)'),
                      ('(-z3pR)', '-(z3 + R)'),
                      ('z3pRpow2', 'square(z3pR)'),
                      ('z3pRpow3', 'cube(z3pR)'),
                      ('z3pRpow2', 'square((-z3pR))'),
                      ('(-z3pRpow3)', 'cube((-z3pR))'),
                      ('(-z3pR)', '((-z3pR))'),
                      ('c1', 'M_PI*lm*mu'),
                      ('c2', 'M_PI*lm')]
        print("""  const RealT
    R3 = R2*R, R4 = R2*R2, R5 = R4*R,
    z111 = z11*z1, z222 = z22*z2, z333 = z33*z3,
    z1111 = z11*z11, z2222 = z22*z22,
    z3pR = z3 + R,
    z3pRpow2 = z3pR*z3pR,
    z3pRpow3 = z3pRpow2*z3pR,
    z3pRpow4 = z3pRpow2*z3pRpow2""")
    elif term == 'C':
        R, T = calc_T_C(lam, mu, x, y)
        text_subs += [('R5', 'pow(R2, 5.0/2.0)'),
                      ('R7', 'pow(R2, 7.0/2.0)'),
                      ('R9', 'pow(R2, 9.0/2.0)'),
                      ('R2', '(R2)'),
                      ('(1.0/R3)', 'pow(R2, -3.0/2.0)'),
                      ('z1111', 'pow(z1, 4)'),
                      ('z2222', 'pow(z2, 4)'),
                      ('z3333', 'pow(z3, 4)'),
                      ('c1', 'M_PI*l2m*mu'),
                      ('c2', 'M_PI*l2m')]
        print("""  const RealT
    x3 = x[2], y3 = y[2],
    R3 = R2*R, R4 = R2*R2, R5 = R4*R, R7 = R5*R2, R9 = R7*R2,
    z111 = z11*z1, z222 = z22*z2, z333 = z33*z3,
    z1111 = z11*z11, z2222 = z22*z22, z3333 = z33*z33,
    c2 = M_PI*l2m,
    c1 = c2*mu;""")
    else:
        raise BaseException('Not a term: {}'.format(term))
    simplify = term in ('A')
    subs = {y1: 0, y2: 0, x1: z1, x2: z2,
            -(x3+y3): z3, x3+y3: -z3, -x3-y3: z3,
            2*x3+2*y3: -2*z3,
            -2*x3-2*y3: 2*z3, -3*x3-3*y3: 3*z3, -4*x3-4*y3: 4*z3,
            -5*x3-5*y3: 5*z3, -7*x3-7*y3: 7*z3,
            (x3+y3)**2: z3**2,
            lam+mu: lm, lam+2*mu: l2m, lam+3*mu: l3m}
    print_calcs_compressed(x, v, d, T, subs, text_subs, simplify=simplify)

def cpp_sigma_halfspace_extra_terms_canonical(term):
    # y = (0,0,y3), x = (x1,0,x3) => z = (x1,0,-(x3+y3))
    x, y, v, d = make_vars()
    var('lam mu')
    text_subs = [('z33', 'square(x3) + 2*x3*y3 + square(y3)'),
                 ('pow(R2,', 'pow(z11 + z33,'),
                 ('(R2)', '(z11 + z33)'),
                 ('R', 'sqrt(R2)'),
                 ('R3', 'pow(R2, 3.0/2.0)')]
    print("""  const Real
    lm = lam + mu,
    l2m = lam + 2*mu;
  const RealT
    d1 = d[0], d2 = d[1], d3 = d[2],
    v1 = v[0], v2 = v[1], v3 = v[2],
    z1 = x1, z3 = -(x3 + y3),
    z11 = z1*z1, z33 = z3*z3, z13 = z1*z3,
    z111 = z11*z1,
    R2 = z11 + z33,
    R = std::sqrt(R2),
    R3 = R2*R;""")
    if term == 'A':
        R, T = calc_T_A(lam, mu, x, y)
        text_subs += [('R4', 'square(R2)'),
                      ('R7', 'pow(R2, 7.0/2.0)'),
                      ('R5', 'pow(R2, 5.0/2.0)'),
                      ('R2', '(R2)'),
                      ('R4', '(pow(z1, 4) + 2*z11*z33 + pow(z3, 4))'),
                      ('R5', 'R*R4'),
                      ('z1111', 'pow(z1, 4)'),
                      ('z3333', 'pow(z3, 4)'),
                      ('R7', 'R*(pow(z1, 6) + 3*z1111*z33 + 3*z11*z3333 + pow(z3, 6))'),
                      ('den1', '(M_PI*l2m*R7)'),
                      ('den2', '(M_PI*l2m*mu*R7)'),
                      ('den3', '(M_PI*l2m*R5)'),
                      ('den4', '(M_PI*l2m*mu*R5)')]
        print("""  const RealT
    R4 = R2*R2, R5 = R4*R, R7 = R5*R2,
    z333 = z33*z3,
    z1111 = z11*z11, z3333 = z33*z33,
    den1 = M_PI*l2m*R7,
    den2 = M_PI*l2m*mu*R7,
    den3 = M_PI*l2m*R5,
    den4 = M_PI*l2m*mu*R5;""")
    elif term == 'B':
        R, T = calc_T_B(lam, mu, x, y)
        text_subs += [('R4', 'square(R2)'),
                      ('R5', 'pow(R2, 5.0/2.0)'),
                      ('R7', 'pow(R2, 7.0/2.0)'),
                      ('R2', '(R2)'),
                      ('(1.0/R)', 'pow(R2, -1.0/2.0)'),
                      ('(1.0/R3)', 'pow(R2, -3.0/2.0)'),
                      ('z1111', 'pow(z1, 4)'),
                      ('z3333', 'pow(z3, 4)'),
                      ('z3pR', '-x3 - y3 + R'),
                      ('z3pRpow4', 'pow(z3pR, 4)'),
                      ('(-z3pR)', '-(z3 + R)'),
                      ('z3pRpow2', 'square(z3pR)'),
                      ('z3pRpow3', 'cube(z3pR)'),
                      ('z3pRpow2', 'square((-z3pR))'),
                      ('(-z3pRpow3)', 'cube((-z3pR))'),
                      ('(-z3pR)', '((-z3pR))')]
        print("""  const RealT
    R4 = R2*R2, R5 = R4*R, R7 = R5*R2,
    z333 = z33*z3,
    z1111 = z11*z11, z3333 = z33*z33,
    z3pR = z3 + R,
    z3pRpow2 = z3pR*z3pR,
    z3pRpow3 = z3pRpow2*z3pR,
    z3pRpow4 = z3pRpow2*z3pRpow2;""")
    elif term == 'C':
        R, T = calc_T_C(lam, mu, x, y)
        text_subs += [('R5', 'pow(R2, 5.0/2.0)'),
                      ('R7', 'pow(R2, 7.0/2.0)'),
                      ('R9', 'pow(R2, 9.0/2.0)'),
                      ('R2', '(R2)'),
                      ('(1.0/R3)', 'pow(R2, -3.0/2.0)'),
                      ('z1111', 'pow(z1, 4)'),
                      ('z3333', 'pow(z3, 4)'),
                      ('c1', 'M_PI*l2m*mu'),
                      ('c2', 'M_PI*l2m')]
        print("""  const RealT
    R4 = R2*R2, R5 = R4*R, R7 = R5*R2, R9 = R7*R2,
    z1111 = z11*z11,
    z333 = z33*z3,
    z3333 = z33*z33,
    c2 = M_PI*l2m,
    c1 = c2*mu;""")
    else:
        raise BaseException('Not a term: {}'.format(term))
    simplify = term in ('A')
    subs = {y1: 0, y2: 0, x2: 0, x1: z1,
            -(x3+y3): z3, x3+y3: -z3, -2*x3-2*y3: 2*z3, -3*x3-3*y3: 3*z3,
            -5*x3-5*y3: 5*z3, -7*x3-7*y3: 7*z3,
            (x3+y3)**2: z3**2,
            x3**2 + 2*x3*y3 + y3**2: z3**2,
            lam+mu: lm, lam+2*mu: l2m, lam+3*mu: l3m}
    print_calcs_compressed(x, v, d, T, subs, text_subs, simplify=simplify)


#cpp_sigma_fullspace()
#cpp_sigma_fullspace_canonical_zv()
#cpp_sigma_halfspace_extra_terms('C')
cpp_sigma_halfspace_extra_terms_canonical('C')
