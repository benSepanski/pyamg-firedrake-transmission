# Setup

```bash
curl -O https://raw.githubusercontent.com/firedrakeproject/firedrake/wence/complex/scripts/firedrake-install
python3 firedrake-install --complex --package-branch firedrake wence/complex
source firedrake/bin/activate
pip install scipy pyamg
```
Install gmsh, then in the `pyamg-firedrake-transmission` directory run the command
```bash
make meshes
```

# Basic Usage

In `transmission.py`, set `wave_number`, `pyamg_tol`, `pyamg_maxiter`, and `mesh_file_name`
to your liking. If you want to save the matrix, make sure `store_mat = True`. Then,
run `python transmission.py`. This should solve the problem using a `pyamg`
solver as a preconditioner. It will save the matrix in `matrix_txt_files` if requested.

To further futz with `pyamg` settings, look at `two_D_helmholtz.py`.

# File Descriptions

## Matrix text files

Matrix text files are stored in `matrix_txt_files/`. The file name
format tells what dimension the helmholtz operator is in
followed by the maximum characteristic length with `'.'` replaced
by a `'%'`. For example, a 2D helmholtz matrix on a mesh
with maximum characteristic length `0.25` would have a matrix
filename of `helmholtz-2D-h0%25.txt`.

## Meshes

Meshes are stored in the `meshes/` directory. There is a 2D mesh
(a circle cut out of a square) and a 3D mesh (a ball cut out of a cube).
In both cases the sphere is the scatterer with a neumann BC
and the boundary of the square has the transmission boundary.
Running `make meshes` builds a bunch of `.msh` files in the
`meshes/circle_in_square` and `meshes/ball_in_cube` directories
(as well as `.options` files detailing the `gmsh` options used, you
can ignore those).
The filename will be `f'max{h}.msh'` where `h` is the maximum
characteristic length as a string, with `'.'` replaced by `'%'`.
For instance, the 2D mesh created with max characteristic length `0.25`
would be the file `meshes/circle_in_square/max0%25.msh`.

To change the number of refinements created of each mesh, change at
the global constants at the beginning of `bin/make_meshes`, then run
`make meshes` again.

## Python Files

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
