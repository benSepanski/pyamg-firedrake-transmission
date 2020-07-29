* transmission.py

    This holds the firedrake code that builds a true solution
    for the mesh and also the transmission weak form.
    It assembles the bilinear form into a matrix A.
    It also creates a preconditioner from the preconditioner
    made in two_D_helmholtz.py

* two_D_helmholtz.py

    Adapted from pyamg examples to be a firedrake preconditioner

* smoothed_aggregation_helmholtz_solver.py

    From pyamg examples with a wildly uneducated shot at planewaves3D
    implementing the 3D waves

* my_vis.py

    From pyamg examples
