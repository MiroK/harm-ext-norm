from weak_bcs.bc_apply import apply_bc
from dolfin import *
from xii import Trace, ii_assemble, EmbeddedMesh
import numpy as np


def get_system(surfaces, f):
    '''f: V -> expression that is f'''
    mesh = surfaces.mesh()
    bmesh = EmbeddedMesh(surfaces, 1)
    
    V = FunctionSpace(mesh, 'CG', 1)
    Q = FunctionSpace(bmesh, 'CG', 1)
    W = [V, Q]
    
    u, p = map(TrialFunction, W)
    v, q = map(TestFunction, W)
    Tu = Trace(u, bmesh)
    Tv = Trace(v, bmesh)

    # The line integral
    dx_ = Measure('dx', domain=bmesh)

    a00 = inner(grad(u), grad(v))*dx
    a01 = inner(Tv, p)*dx_
    a10 = inner(Tu, q)*dx_

    L0 = inner(Constant(0), v)*dx
    L1 = inner(Trace(f, bmesh), q)*dx_

    a = [[a00, a01], [a10, 0]]
    L = [L0, L1]

    A, b = map(ii_assemble, (a, L))

    V_bcs = [DirichletBC(V, Constant(0), surfaces, 2)]
    bcs = [V_bcs, []]
    # The system matrix with bcs
    A, b = apply_bc(A, b, bcs)

    return A, b, W


def make_lhs(r0, sine_series, cosine_series):
    '''
    V -> solution of -Delta u = 0 on a anulus (r0, 1) with 0 bcs on r=1
    and u(r0, theta) = f given as sine series, cosine series
    '''
    def make_u(V, sine_series=sine_series, cosine_series=cosine_series):
        x, y = V.tabulate_dof_coordinates().reshape((V.dim(), -1)).T
        
        theta = np.arctan2(y, x)
        r = np.sqrt(x**2 + y**2)

        # Sulution coefs are (sin(k theta), f(theta))_(0, 2*pi)
        #                    (cos(k theta), f(theta))_(0, 2*pi)
        #
        # so if f is fourier series the coefs simflify a lot
        values = np.zeros_like(x)
        for k, ak in enumerate(sine_series, 1):
            Ak = ak/(r0**k - 1./(r0**k))
            values += Ak*np.sin(k*theta)*(r**k - 1./(r**k))
            
        # Log term
        cosine_series = iter(cosine_series)
        try:
            b0 = next(cosine_series)
            F = b0/np.log(r0)

            values += F*np.log(r)
        except StopIteration:
            pass
            
        for k, bk in enumerate(cosine_series, 1):
            Bk = bk/(r0**k - 1./(r0**k))
            values += Bk*np.cos(k*theta)*(r**k - 1./(r**k))

        u = Function(V)
        u.vector().set_local(values)

        return u
    return make_u
            

def make_rhs(sine_series, cosine_series):
    '''
    V -> f in -Delta u = 0 on a anulus (r0, 1) with 0 bcs on r=1
    and u(r0, theta) = f given as sine series, cosine series
    '''
    def make_f(V, sine_series=sine_series, cosine_series=cosine_series):
        x, y = V.tabulate_dof_coordinates().reshape((V.dim(), -1)).T
        
        theta = np.arctan2(y, x)

        # so if f is fourier series the coefs simflify a lot
        values = np.zeros_like(x)
        for k, ak in enumerate(sine_series, 1):
            values += ak*np.sin(k*theta)
            
        for k, bk in enumerate(cosine_series, 0):
            values += bk*np.cos(k*theta)

        f = Function(V)
        f.vector().set_local(values)

        return f
    return make_f

# --------------------------------------------------------------------

if __name__ == '__main__':
    from xii import ii_convert, ii_Function
    from find_E import get_mesh    
    import argparse
    
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # Number of mesh refinements to use in convergence study
    parser.add_argument('-dx', type=float, default=0.25, help='Inner half width')
    parser.add_argument('-size0', type=float, default=0.5, help='Initial size of the mesh')
    parser.add_argument('-nrefs', type=int, default=3, help='Number of refinements')

    args = parser.parse_args()

    dx_, size = args.dx, args.size0
    table  = []

    sine_series = (1, 0.1, 0.3)
    cosine_series = (1, 0.2)

    make_u = make_lhs(args.dx, sine_series, cosine_series)
    make_f = make_rhs(sine_series, cosine_series)

    errors, hs = [], []
    for i in range(args.nrefs):
        surfaces = get_mesh('circle', dx_, size)
        size /= 2.

        mesh = surfaces.mesh()

        dsI = Measure('ds', subdomain_data=surfaces, subdomain_id=1)
        dsO = Measure('ds', subdomain_data=surfaces, subdomain_id=2)

        V = FunctionSpace(mesh, 'DG', 5)
        # Data in some higher order space
        uT = make_u(V)
        f = make_f(V)

        # We at least get the bcs right
        # print sqrt(abs(assemble(inner(uT-f, uT-f)*dsI)))
        # print sqrt(abs(assemble(inner(uT-Constant(0), uT-Constant(0))*dsO)))

        A, b, W = get_system(surfaces, f)

        # Solve
        A, b = map(ii_convert, (A, b))

        wh = ii_Function(W)

        solve(A, wh.vector(), b)
        uh = wh[0]
        
        e = errornorm(uT, uh, 'H10', degree_rise=2)
        h = uh.function_space().mesh().hmin()
        if errors:
            rate = ln(e/errors[-1])/ln(h/hs[-1])
        else:
            rate = np.nan
        errors.append(e)
        hs.append(h)

        print 'ndofs = %d h = %.2E e = %.4E r = %.2f' % (A.size(0), h, e, rate)

    e = interpolate(uh, V)
    e.vector().axpy(-1, uT.vector())
    File('plots/error_dx%g.pvd' % args.dx) << e
