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

# Matrix text files

Matrix text files are stored in `matrix_txt_files/`. The file name
format tells what dimension the helmholtz operator is in
followed by the maximum characteristic length with `'.'` replaced
by a `'%'`. For example, a 2D helmholtz matrix on a mesh
with maximum characteristic length `0.25` would have a matrix
filename of `helmholtz-2D-h0%25.txt`.


# File Descriptions

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
