import os, sys, shutil, subprocess
from weak_bcs.bc_apply import apply_bc
from weak_bcs.utils import matrix_fromHs
from dolfin import *
from xii import Trace, ii_assemble, EmbeddedMesh, block_diag_mat
from xii.linalg.convert import numpy_to_petsc
from block.algebraic.petsc import Cholesky
from scipy.linalg import sqrtm
from hsmg.hseig import HsNorm
from tqdm import trange
import numpy as np


def find_E_mixed(surfaces):
    '''Try via the mixed problem'''
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
    L1 = inner(Constant(0), q)*dx_

    a = [[a00, a01], [a10, 0]]
    L = [L0, L1]

    A, b = map(ii_assemble, (a, L))

    V_bcs = [DirichletBC(V, Constant(0), surfaces, 2)]
    bcs = [V_bcs, []]
    # The system matrix with bcs
    A, b = apply_bc(A, b, bcs)

    Hs = matrix_fromHs(HsNorm(Q, s=-0.5))
    # Inner product matrix for Riesz Map
    B = block_diag_mat([A[0][0], Hs])

    return A, B


def find_E_def(surfaces):
    '''From def |u(f)|/|f|'''
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
    L1 = inner(Constant(0), q)*dx_

    a = [[a00, a01], [a10, 0]]
    L = [L0, L1]

    A, b = map(ii_assemble, (a, L))

    V_bcs = [DirichletBC(V, Constant(0), surfaces, 2)]
    bcs = [V_bcs, []]
    # The system matrix with bcs
    A, b = apply_bc(A, b, bcs)

    Hs = matrix_fromHs(HsNorm(Q, s=0.5))

    # Get the Schur complement
    [[A, Bt],
     [B, _]] = A

    S = np.zeros((Q.dim(), Q.dim()))
    
    Ainv = Cholesky(A)
    # Build it column by column
    f = Function(Q)
    f_vec = f.vector()
    f_vals = np.zeros(Q.dim(), dtype=float)
    for col in trange(Q.dim(), desc='Build Schur complement'):
        f_vals[col] = 1.
        f_vec.set_local(f_vals)

        S[:, col] = (B*Ainv*Bt*f_vec).get_local()
        f_vals *= 0
    
    # S = sqrtm(M.dot(np.linalg.inv(S).dot(M.T))).real
    S = numpy_to_petsc(S)
    
    # Inner product matrix for Riesz Map
    Hs = matrix_fromHs(HsNorm(Q, s=-0.5))

    return S, Hs


def get_mesh(base, dx, size):
    '''Load or build'''
    if not os.path.exists('meshes'):
        os.mkdir('meshes')

    root = '%s_dx%g_size%g' % (base, dx, size)
    xml = 'meshes/%s.xml' % root
    surfaces = 'meshes/%s_facet_region.xml' % root
    # Load it
    if os.path.exists(xml):
        mesh = Mesh(xml)
        surfaces = MeshFunction('size_t', mesh, surfaces)

        volume = sum(cell.volume() for cell in cells(mesh))
                    
        if base == 'square':
            error = abs(volume - (1 - (2*dx)**2))
            tol = 1E-10
        elif base == 'circle':
            error = abs(volume - (np.pi - np.pi*(dx)**2))
            tol = 1E-2
        assert error < tol, error

        print '\t>>', volume, error
        
        return surfaces

    # Make it
    msh = '.'.join([root, 'msh'])
    args = (base, dx, size, msh)
    subprocess.call(['gmsh -2 %s.geo -setnumber dx %g -clscale %g -o %s' % args], shell=True)

    subprocess.call(['dolfin-convert %s %s' % (msh, '.'.join([root, 'xml']))], shell=True)
    
    shutil.move(os.path.basename(xml), xml)
    shutil.move(os.path.basename(surfaces), surfaces)
    
    os.remove(msh)
    map(os.remove, filter(lambda x: x.endswith('.xml'), os.listdir('.')))
    
    return get_mesh(base, dx, size)
    
# --------------------------------------------------------------------

if __name__ == '__main__':
    from scipy.linalg import eigvalsh
    from xii import ii_convert
    import argparse

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # Number of mesh refinements to use in convergence study
    parser.add_argument('-dx', type=float, default=0.25, help='Inner half width')
    parser.add_argument('-size0', type=float, default=0.5, help='Initial size of the mesh')
    parser.add_argument('-nrefs', type=int, default=3, help='Number of refinements')
    parser.add_argument('-domain', type=str, default='square', help='Geometry file')

    args = parser.parse_args()

    find_E = {'mixed': find_E_mixed,
              'def': find_E_def}['def']
    
    dx_, size = args.dx, args.size0
    table  = []
    msg = '%d eigs in [%g, %g] U [%g, %g]'
    for i in range(args.nrefs):
        surfaces = get_mesh(args.domain, dx_, size)
        size /= 2.

        A, B = map(ii_convert, find_E(surfaces))
        ndofs = A.size(0)
        
        eigw = eigvalsh(A.array(), B.array())

        pos = np.where(eigw > 0)[0]
        if len(pos):
            max_pos = max(eigw[pos])
            min_pos = min(eigw[pos])
        else:
            max_pos, min_pos = np.nan, np.nan

        neg = np.where(eigw < 0)[0]
        if len(neg):
            max_neg = max(eigw[neg])
            min_neg = min(eigw[neg])
        else:
            max_neg, min_neg = np.nan, np.nan

        row = [ndofs, min_neg, max_neg, min_pos, max_pos]
        table.append(row)
        # Always review
        print
        for row in table:
            print msg % tuple(row)

    if not os.path.exists('results'):
        os.mkdir('results')

    File('results/%s_mesh_dx_%g.pvd' % (args.domain, dx_)) << surfaces
        
    np.savetxt('results/%s_dx_%g' % (args.domain, dx_), np.array(table))
