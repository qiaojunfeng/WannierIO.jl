## Sparse Binary Storage For `tb.dat`

`tb.dat` files can be converted to sparse-packed binary formats to reduce storage size.
The data container is `TbDat`.

Supported formats:

- `HDF5Format`
- `JLD2Format`
- `ZarrFormat`

Unified API:

```julia
tb = read_w90_tb("Si2_valence_tb.dat")
write_w90_tb("Si2_valence.h5", tb)
tb_h5 = read_w90_tb("Si2_valence.h5")
```

Or specify the format explicitly:

```julia
write_w90_tb("out.any", tb, JLD2Format(); atol=1e-10)
tb2 = read_w90_tb("out.any", JLD2Format())
```

To store reduced precision data, pass `value_type`:

```julia
write_w90_tb("out.h5", tb, HDF5Format(); value_type=Float32, index_type=Int16)
tb32 = read_w90_tb("out.h5", HDF5Format())
```

The shared sparse conversion API is:

```julia
pack = sparsify(tb; atol=1e-10, value_type=Float32)
tb2 = densify(pack)
```

For each operator (`H`, `r_x`, `r_y`, `r_z`), sparse matrices are stored in CSC form.
HDF5/Zarr backends persist the following arrays per operator:

- `*_is_complex`: whether nonzeros are complex-valued
- `*_nzptr`: offsets into concatenated nonzero arrays for each matrix
- `*_colptr`: column pointers for each matrix, stacked as a 2D array
- `*_rowval`: concatenated row indices
- `*_nzre`: concatenated real parts of nonzeros
- `*_nzim`: concatenated imaginary parts of nonzeros (stored only when `*_is_complex` is true)

Benchmark utility script (HDF5 vs JLD2 vs Zarr):

```bash
julia --project misc/benchmark_tb_storage.jl /path/to/tb_dir 200 1e-10
```

These features require optional packages:

```julia
using Pkg
Pkg.add(["HDF5", "JLD2", "Zarr"])

# Load the package to activate the corresponding extension method.
using HDF5
```
