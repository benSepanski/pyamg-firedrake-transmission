import firedrake.variational_solver as vs
import numpy as np

from firedrake import (
    ln, pi, Mesh, SpatialCoordinate, sqrt,
    FunctionSpace, Function, TrialFunction, TestFunction,
    FacetNormal, inner, grad, dx, ds, Constant,
    assemble
    )
from math import factorial
from two_D_helmholtz import AMGTransmissionPreconditioner


def hankel_function(expr, n):
    """
        Returns a :mod:`firedrake` expression approximation a hankel function
        of the first kind and order 0
        evaluated at :arg:`expr` by using the taylor
        series, expanded out to :arg:`n` terms.
    """
    j_0 = 0
    for i in range(n):
        j_0 += (-1)**i * (1 / 4 * expr**2)**i / factorial(i)**2

    g = Constant(0.57721566490153286)
    y_0 = (ln(expr / 2) + g) * j_0
    h_n = 0
    for i in range(n):
        h_n += 1 / (i + 1)
        y_0 += (-1)**(i) * h_n * (expr**2 / 4)**(i+1) / (factorial(i+1))**2
    y_0 *= Constant(2 / pi)

    imag_unit = Constant((np.zeros(1, dtype=np.complex128) + 1j)[0])
    h_0 = j_0 + imag_unit * y_0
    return h_0


mesh_file_name = 'meshes/circle_in_square/max0%25.msh'
h = mesh_file_name[mesh_file_name.find('%')-1:mesh_file_name.find('.')]
mesh = Mesh(mesh_file_name)
x, y = SpatialCoordinate(mesh)
if mesh.geometric_dimension() == 2:
    scatterer_bdy_id = 1
    outer_bdy_id = 2
else:
    scatterer_bdy_id = 1
    outer_bdy_id = 3


fspace = FunctionSpace(mesh, 'CG', 1)

u = TrialFunction(fspace)
v = TestFunction(fspace)

wave_number = 0.1  # one of [0.1, 1.0, 5.0, 10.0]
true_sol = Constant(1j / 4) * hankel_function(wave_number * sqrt(x**2 + y**2),
                                              n=80)

a = inner(grad(u), grad(v)) * dx \
    - Constant(wave_number**2) * inner(u, v) * dx \
    - Constant(1j * wave_number) * inner(u, v) * ds(outer_bdy_id)

n = FacetNormal(mesh)
L = inner(inner(grad(true_sol), n), v) * ds(scatterer_bdy_id)

solution = Function(fspace)

# Create a solver and return the KSP object with the solution so that can get
# PETSc information
# Create problem
problem = vs.LinearVariationalProblem(a, L, solution)

# Create solver and call solve
pyamg_tol = None
pyamg_maxiter = None
solver_params = {'pyamg_tol': pyamg_tol,
                 'pyamg_maxiter': pyamg_maxiter,
                 'ksp_monitor': None,
                 }
solver = vs.LinearVariationalSolver(problem,
                                    solver_parameters=solver_params)
# prepare to set up pyamg preconditioner if using it
A = assemble(a).M.handle
store_mat = True
if store_mat:
    dim = mesh.geometric_dimension()
    file_name = f"matrix_txt_files/helmholtz-{dim}D-h{h}.txt"
    from firedrake.petsc import PETSc
    myviewer = PETSc.Viewer().createASCII(
        file_name,
        mode=PETSc.Viewer.Format.ASCII_COMMON,
        comm=PETSc.COMM_WORLD)
    A.view(myviewer)

# solve
pc = solver.snes.getKSP().pc
pc.setType(pc.Type.PYTHON)
pc.setPythonContext(AMGTransmissionPreconditioner(wave_number,
                                                  fspace,
                                                  A,
                                                  tol=pyamg_tol,
                                                  maxiter=pyamg_maxiter,
                                                  use_plane_waves=True))

solver.solve()
