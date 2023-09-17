from sympy import *

def matrix_example():
    var('a b c d x y')
    M = Matrix([[a, b], [c, d]])
    v = Matrix([x, y])
    pretty_print(M*v)
    pretty_print(M.T*v)

def examples():
    x, t = symbols('x t')

    a = x**2 + t
    pretty_print(a)
    pretty_print(a.subs(x, 2*x))

    print()
    a = Integral(cos(x)*exp(x), x).doit()
    pretty_print(a)
    print('sexr:')
    print(srepr(a))
    print('preorder_traversal:')
    for e in preorder_traversal(a): print(e)
    print('postorder_traversal:')
    for e in postorder_traversal(a): print(e)

    print()
    a = integrate(pi*sin(x**2), (x, -oo, oo))
    pretty_print(a)

    print()
    y = Function('y')
    a = dsolve(Eq(y(t).diff(t,t) - y(t), exp(t)), y(t))
    pretty_print(a)

    print()
    a = Matrix([[1, 2], [2, 2]]).eigenvals()
    pretty_print(a)

    print()
    pretty_print(Eq(Integer(1)/Integer(3), Rational(1,3)))
    pretty_print(Integer(1)/Integer(3) == Rational(1,3))

    print()
    pretty_print(factor(x**2 + 4*x + 4))
    pretty_print(expand((x+2)**2))
    pretty_print(simplify(2*x + 4*x))

    matrix_example()

def numbers_only(e):
    for a in e.args:
        if type(a) not in (Integer, Rational):
            return False
    return True

def contains_symbol(e, x):
    for es in preorder_traversal(e):
        if es == x:
            return True
    return False

def leaf_terms_only(e):
    if numbers_only(e) or type(e) == Symbol: return True
    if type(e) != Add: return False
    for a in e.args:
        for es in preorder_traversal(a):
            if type(es) == Add:
                return False
    return True

def simplify_pieces(e, leaf_fn=factor, verbose=False):
    "Attempt to factor pieces of e in a reasonable way."
    if leaf_terms_only(e):
        el = leaf_fn(e)
        if verbose and not numbers_only(e):
            print(f'simplify_pieces:\n  {e}\nto\n  {el}')
        return el
    else:
        new_args = []
        for i in range(len(e.args)):
            new_args.append(simplify_pieces(e.args[i], leaf_fn, verbose))
        return e.func(*new_args)

def dev_simplify_pieces():
    var('x', real=True)
    e = 1 + (x**2 + 4*x + 4)/(x**2 - 2*x + 1) + (x - 1)**2 + 1/(x**2 - x - 12)
    print(srepr(e))
    print('simplify(factor(e):')
    pretty_print(simplify(factor(e)))
    print('simplify_pieces:')
    pretty_print(simplify_pieces(e, verbose=True))

def make_s(x, y):
    """Setting 1. Full space. Plane strain. slip = (1, 0). Sigma at receiver for
    one out-of-page long dislocation (oopld) at (0, 0). Corresponds to Segall
    book eqs 3.77-79 with r2 -> inf and omitting mu/2 pi (1-nu)."""
    r2 = x**2 + y**2
    sxx = -y*(y**2 + 3*x**2) / r2**2
    syy = -y*(y**2 -   x**2) / r2**2
    sxy = -x*(y**2 -   x**2) / r2**2
    return ('sxx', sxx), ('syy', syy), ('sxy', sxy)

def c1():
    "Calculate point Green's function for setting 1."
    var('x y z')
    s = make_s(x, y)
    p = [limit((s[i][1].subs(x, x-z) -
                s[i][1].subs(x, x+z)) /
               (2*z),
               z, 0, dir='+').factor()
         for i in range(len(s))]
    for i in range(len(s)):
        print('p' + s[i][0][1:] + ' =', p[i])
    # Equivalently, just take the derivative:
    print('Using partial deriv:')
    for i in range(len(s)):
        print('p' + s[i][0][1:] + ' =', -simplify(s[i][1].diff(x)))

def make_p(x, y):
    pxx = -2*x*y*(3*x**2 - y**2)/(x**2 + y**2)**3
    pyy =  2*x*y*(x**2 - 3*y**2)/(x**2 + y**2)**3
    pxy =  (x**2 - 2*x*y - y**2)*(x**2 + 2*x*y - y**2)/(x**2 + y**2)**3
    return ('pxx', pxx), ('pxy', pxy), ('pyy', pyy)

def fsps_integral_simplify(e):
    if not e.args: return e
    if type(e) == Pow:
        return Pow(factor(expand(e.args[0])), e.args[1])
    args = []
    for a in e.args:
        args.append(fsps_integral_simplify(a))
    return e.func(*args)

def z_leaf_fn(e):
    if not e.args: return e
    e = e.collect(y, evaluate=True)
    if type(e) != Add: return e
    #print('collect',e)
    args_no_y = []
    args_y = []
    for a in e.args:
        if (contains_symbol(a, y)):
            args_y.append(factor(a))
        else:
            args_no_y.append(factor(a))
    #print('args_y:',args_y)
    #print('args_no_y:',args_no_y)
    return Add(factor(Add(*args_no_y)), Add(*args_y))

def c2():
    """
    we want
        sigma(x, y) = f(d; x, y) - f(-d; x, y),
    where the antiderivative
        f(z; x, y) = int z^k p(x - z, y) dz.
    calculate f for each component and k."""
    var('x y z', real=True)
    ps = make_p(x, y)
    for i in range(len(ps)):
        p = ps[i][1]
        print('-> sigma' + ps[i][0][1:] + ':')
        for k in range(4):
            ip = integrate(p.subs(x, x-z)*z**k, z)
            ip = fsps_integral_simplify(ip)
            ip = simplify_pieces(ip, leaf_fn=z_leaf_fn, verbose=False)
            print('  s(z) = z^' + str(k) + ':')
            print(' '*4 + 'f(z) = ' + str(ip).replace('**', '^'))

def c3():
    # For lineflt_fsps_calc_sn_traction.
    var('x y z tc D s0 s1', real=True)
    dh = D/2
    x = tc*D-dh
    ss = Matrix([s0, s1])
    ps = make_p(x, y)
    p = ps[1][1].subs(y, 0).subs(x, x-z)
    print('p =', p)
    e = 0
    for k in (0,1):
        print(k)
        q = p*(z-x)**k
        print(q)
        f = integrate(q, z)
        f = fsps_integral_simplify(f)
        f = simplify_pieces(f, leaf_fn=z_leaf_fn, verbose=False)
        e += ss[k]*(f.subs(z, dh) - f.subs(z, -dh))
    print(e)

# We used essentially this code with sage. In pure sympy, I'm so far unable to
# get a solution.
def c7old():
    """Calculate f(th) in the normal traction f(th2) - f(th1) for a receiver on
    part of a circle. The circle is centered at (0, 0). th = 0 is at its bottom.
    """
    var('r x y sx sy rx ry th c s', real=True)

    # The receiver is at the bottom of the circle and, since it's on the circle,
    # oriented parallel to the x axis.
    rx = 0
    ry = -r

    c = cos(th)
    s = sin(th)
    # Given th, the angle location on the circle, the source's (x, y) is:
    sx =  r*s
    sy = -r*c
    # The source's orientation is, conveniently, th. Get the rotation matrix for
    # th.
    R = Matrix([[c, -s], [s, c]])
    # Use it to transform the receiver point into source-local coords. This next
    # line both makes the local origin at the source and does the rotation.
    rcv_local = R.T*Matrix([rx - sx, ry - sy])
    x = rcv_local[0]
    y = rcv_local[1]
    # Compute the point Green's functions in source-local coords.
    pxx, pxy, pyy = make_p(x, y)
    pxx = pxx[1]; pxy = pxy[1]; pyy = pyy[1]
    P = Matrix([[pxx, pxy], [pxy, pyy]])
    # Rotate the tensor into global coords.
    p = R*P*R.T
    # The receiver is at angle 0, so no further rotation is needed to get the
    # traction components on the receiver. At least for now, I want just the
    # normal traction.
    integrand_normal = simplify(p[1,1])
    print('integrand_normal =', integrand_normal)
    return
    arc_length = r
    print('attempting to integrate')
    ip = integrate(integrand_normal*arc_length, th)
    print('done')
    print('f_yy(th) =', simplify(ip))

def example_hfp():
    var('x y', real=True)
    integrand = -(x**2 - 2*x*y - y**2)*(x**2 + 2*x*y - y**2)/(x**2 + y**2)**3
    v = integrate(integrand, (x,-1,1))
    print(simplify(v))

def cquadratic():
    #   The fault is a quadratic f(x) = a x^2 centered at and tangent to rcv =
    # (0,0). The parameterization of this curve is p(x) = (x, f(x)). The global
    # coords of the source are (sx,f(sx)).
    var('x a', real=True)
    # src parabola data
    sx = x
    sy = a*x**2
    psx1 = 1      # x'
    psx2 = 2*a*x  # f'(sx) x'
    h = sqrt(psx1**2 + psx2**2)
    cth = psx1/h
    sth = psx2/h
    # rotation for glb to/from src-lcl coords
    R = Matrix([[cth, -sth], [sth, cth]])
    # local coords of rcv and src: R'(r_gbl - s_gbl)
    r_sel = R.T*Matrix([-sx, -sy])
    # Green's fn in src-local coords
    pxx, pxy, pyy = make_p(r_sel[0], r_sel[1])
    pxx = pxx[1]; pxy = pxy[1]; pyy = pyy[1]
    P = Matrix([[pxx, pxy], [pxy, pyy]])
    # Rotate the tensor into global coords.
    p = R*P*R.T
    # The receiver is at angle 0, so no further rotation is needed to get the
    # traction components on the receiver. At least for now, I want just the
    # normal traction.
    for k in (0, 1):
        print('k', k)
        integrand_normal = simplify(p[k,k])
        arc_length = h
        integrand = integrand_normal*arc_length
        print('attempting to integrate\n  ', integrand)
        continue
        v = integrate(integrand, x)
        print('done')
        print('integral =', simplify(v))

c1()
