from dolfin import *
from scipy.linalg import eigvalsh, eigh
from scipy.sparse import csr_matrix, csc_matrix
from slepc4py import SLEPc
from petsc4py import PETSc    
import numpy as np


# FIXME: SLEPc settings
def find_S(surfaces):
    '''
    Steklov eigenvalue problem for estimating the trace constant for the 
    L^2 space

      |Tu|_{0, \Gamma} \leq |\nabla u|_{0, \Omega} 

    -Delta u = 0 on \Omega
    grad(u).n = lambda*u on \Gamma.
    '''
    mesh = surfaces.mesh()
    
    V = FunctionSpace(mesh, 'CG', 1)
    u, v = TrialFunction(V), TestFunction(V)

    ds = Measure('ds', domain=mesh, subdomain_data=surfaces)

    a = inner(grad(u), grad(v))*dx
    m = inner(u, v)*ds(1) + Constant(0)*inner(u, v)*dx
    L = inner(Constant(0), v)*dx

    bcs = DirichletBC(V, Constant(0), surfaces, 2)

    A, _ = assemble_system(a, L, bcs)
    B, _ = assemble_system(m, L, bcs)

    bcs = DirichletBC(V, Constant(0), surfaces, 1)

    interior_dofs = set(range(V.dim())) - set(bcs.get_boundary_values().keys())

    return A, B, V, interior_dofs


def build_basis(V, dofs):
    '''Nodal basis for dofs of V'''
    v = Function(V).vector()
    values = np.zeros(V.dim())

    basis = []
    for d in dofs:
        values[d] = 1
        v_d = v.copy()
        v_d.set_local(values)
        basis.append(as_backend_type(v_d).vec())

        values *= 0
    return basis


def eig_min(A, B, V, idofs, opts=None):
    '''
    Smallest eigenvalue of A u = lmbda B u where B is positive semidefinite 
    with Z the nullspace of B
    '''
    # Defaults options
    defaults = {'-eps_max_it': 20000,
                '-eps_nev': 3,
                '-eps_monitor': None,
                '-eps_type': 'krylovschur',
                '-st_ksp_rtol': 1E-8,
                '-st_ksp_monitor_true_residual': None}
    
    if opts is not None:
        for key, value in opts.items():
            defaults[key] = value if value is not 'none' else None

    opts = PETSc.Options()
    [opts.setValue(k, v) for k, v in defaults.items()]
    
    A = as_backend_type(A).mat()
    B = as_backend_type(B).mat()

    Z = build_basis(V, idofs)

    E = SLEPc.EPS().create()
    E.setOperators(A, B)
    
    E.setProblemType(SLEPc.EPS.ProblemType.GHEP)
    E.setWhichEigenpairs(SLEPc.EPS.Which.SMALLEST_MAGNITUDE)

    ST = E.getST()
    ST.setType('sinvert')

    E.setDeflationSpace(Z)

    E.setFromOptions()

    print E.getTolerances()
    
    E.solve()

    its = E.getIterationNumber()
    nconv = E.getConverged()
    
    print 'Num converged %d of size(A) %g' % (nconv, A.size)
    
    eigen_pairs = []
    for i in range(nconv):
        eigv = Function(V)
        E.getEigenvector(i, as_backend_type(eigv.vector()).vec())

        eigw = E.getEigenvalue(i).real

        eigen_pairs.append((eigw, eigv))
        
    return eigen_pairs


def petsc_to_csr(A):
    '''PETScMatrix -> csr'''
    return csr_matrix(as_backend_type(A).mat().getValuesCSR()[::-1],
                      shape=(A.size(0), A.size(1)))

def identity(n, idx):
    '''n by len(idx) matrix that is the identity for index space'''
    indptr = np.arange(len(idx)+1)
    indices = np.fromiter(idx, dtype=int)
    data = np.ones_like(indices)
    
    return csc_matrix((data, indices, indptr), shape=(n, len(idx)))

# --------------------------------------------------------------------

if __name__ == '__main__':
    import matplotlib.pyplot as plt
    from find_E import get_mesh
    import argparse, os

    # NOTE: -eps_type 'lapack'
    # gives all the eigenvalues via LAPACK - and it work beutifully
    # with deflation

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # Number of mesh refinements to use in convergence study
    parser.add_argument('-dx', type=float, default=0.25, help='Inner half width')
    parser.add_argument('-size0', type=float, default=0.5, help='Initial size of the mesh')
    parser.add_argument('-nrefs', type=int, default=3, help='Number of refinements')

    args, unknown = parser.parse_known_args()

    opts = dict(zip(unknown[::2], unknown[1::2]))
    dx_, size = args.dx, args.size0

    table = []
    for i in range(args.nrefs):
        surfaces = get_mesh('circle', dx_, size)
        size /= 2.

        A, B, V, idofs = find_S(surfaces)

        pairs = eig_min(A, B, V, idofs, opts)
        spectrum = [p[0] for p in pairs]

        row = [A.size(0), np.nan, np.nan, spectrum[0], np.nan]
        table.append(row)
        print row
        
    File('plots/stekhlov_dx%g.pvd' % args.dx) << pairs[0][1]

    if not os.path.exists('results'): os.mkdir('results')
        
    np.savetxt('results/steklov_%s_dx_%g' % ('circle', dx_), np.array(table))
