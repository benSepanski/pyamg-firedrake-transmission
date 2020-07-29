* transmission.py

    - holds the firedrake code to build a bilinear form
      for transmission
    - assembles the bilinear form into a matrix A.
    - creates a preconditioner from the preconditioner class
        made in two_D_helmholtz.py

* two_D_helmholtz.py

    Adapted from pyamg examples to be a firedrake preconditioner

* smoothed_aggregation_helmholtz_solver.py

    From pyamg examples with a wildly uneducated shot at planewaves3D
    implementing the 3D waves

* my_vis.py

    From pyamg examples
