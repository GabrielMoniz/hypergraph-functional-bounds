# hypergraph-functional-bounds

Julia code for computing the **local (classical)** and **quantum** bounds of the Bell functional associated with *k*-uniform complete hypergraph states, as constructed in:

> G. Moniz, *"k-Uniform complete hypergraph states stabilizers in terms of local operators"*, submitted to *Physical Review A* (arXiv preprint, 2025).

---

## Background

A *k-uniform complete hypergraph* on *n* vertices has vertex set {1, …, *n*} and hyperedge set *E* consisting of **all** *k*-element subsets of {1, …, *n*}. The associated quantum state is

$$
|G_{n,k}\rangle = \prod_{e \in E} \mathrm{CZ}_e\, |{+}\rangle^{\otimes n},
$$

where $\mathrm{CZ}_e$ denotes the generalized controlled-*Z* gate acting on the qubits in hyperedge *e*.

The stabilizer generators of $|G_{n,k}\rangle$ are

$$
g_j = X_j \prod_{\substack{e \in E \\ j \in e}} \mathrm{CZ}_{e \setminus \{j\}}, \qquad j = 1, \ldots, n.
$$

The central result of the manuscript is an explicit expression for each $g_j$ as a linear combination of **tensor products of local Pauli operators** (single-qubit *X* and *Z* gates). This Pauli expansion defines a multilinear Bell functional

$$
\beta_{n,k} = \sum_{j=1}^{n} \langle g_j \rangle,
$$

whose **quantum bound** is achieved by $|G_{n,k}\rangle$ and whose **local bound** is computed here.

---

## Contents

| File | Description |
|------|-------------|
| `bounds.jl` | Main script: constructs the Bell operator, builds the correlation tensor, and computes both bounds |
| `Project.toml` | Julia project specification (dependencies and compatibility) |
| `LICENSE` | MIT license |
| `CITATION.cff` | Machine-readable citation metadata |

---

## Requirements

* [Julia](https://julialang.org/) ≥ 1.9
* [Ket.jl](https://github.com/araujoms/ket) ≥ 0.9 — quantum-information toolkit (provides `local_bound` and `tsirelson_bound`)
* [Combinatorics.jl](https://github.com/JuliaMath/Combinatorics.jl)

---

## Installation

1. **Clone the repository**:

   ```bash
   git clone https://github.com/GabrielMoniz/hypergraph-functional-bounds.git
   cd hypergraph-functional-bounds
   ```

2. **Install dependencies** using Julia's built-in package manager. From the Julia REPL:

   ```julia
   using Pkg
   Pkg.activate(".")        # activate the project environment
   Pkg.instantiate()        # install all declared dependencies
   ```

   Or equivalently from the shell:

   ```bash
   julia --project=. -e 'using Pkg; Pkg.instantiate()'
   ```

---

## Usage

Run the main script from the command line:

```bash
julia --project=. bounds.jl
```

Or load it interactively from the Julia REPL:

```julia
using Pkg; Pkg.activate(".")
include("bounds.jl")
```

The script will print a table of local and quantum bounds for several (*n*, *k*) pairs.

### Calling individual functions

After `include("bounds.jl")` (or after copying the function definitions into the REPL), you can call the individual functions directly:

```julia
# Exact quantum bound (largest eigenvalue of the Bell operator)
qb = quantum_bound_exact(4, 3)

# Local (classical) bound via Ket.jl branch-and-bound
lb = local_bound_lhv(4, 3)

# NPA semidefinite upper bound at hierarchy level 1
npa = quantum_bound_npa(4, 3; level = 1)

# Construct the n-qubit k-uniform complete hypergraph state vector
ψ = hypergraph_state(4, 3)

# Obtain the Bell operator as a Hermitian matrix
B = bell_operator(4, 3)

# Correlation tensor in Ket.jl notation
FC = bell_functional_tensor(4, 3)
```

---

## Method

### Quantum bound

The quantum bound is computed as the **largest eigenvalue** of the Bell operator

$$
B_{n,k} = \sum_{j=1}^{n} g_j,
$$

which is a $2^n \times 2^n$ Hermitian matrix assembled from the stabilizer generators using `LinearAlgebra.eigvals`.

Because $|G_{n,k}\rangle$ is a common +1 eigenstate of all $g_j$, the state saturates the quantum bound $\beta_Q = n$.

### Local bound

The local (classical) bound is the maximum of $\beta_{n,k}$ over all **deterministic local strategies**: assignments of a fixed outcome $\pm 1$ to each measurement setting of each party. It is computed exactly by [`Ket.local_bound`](https://github.com/araujoms/ket), which implements the efficient branch-and-bound algorithm of

> Araújo, Hirsch, Quintino, *A quantum advantage for inferring causal structure*, arXiv:2005.13418.

The Bell functional is encoded as a correlation tensor `FC` with shape `(3, 3, …, 3)` (*n* indices, each of size 3), where

* index 1 → identity (marginal),
* index 2 → first dichotomic observable (*Z*),
* index 3 → second dichotomic observable (*X*).

### NPA bound

As an independent check on the quantum bound, the script also calls [`Ket.tsirelson_bound`](https://github.com/araujoms/ket), which implements the Navascués–Pironio–Acín (NPA) semidefinite-programming hierarchy

> Navascués, Pironio, Acín, *A convergent hierarchy of semidefinite programs characterizing the set of quantum correlations*, New J. Phys. **10**, 073013 (2008).

at NPA level 1 (the Tsirelson bound in two-dichotomic-observable scenarios).

---

## Expected output

```
k-Uniform Complete Hypergraph Bell Functional Bounds
=================================================
  n    k   Quantum (exact)   Local (LHV)
-------------------------------------------------
  3    2          3.000000      3.000000
  4    2          4.000000      4.000000
  5    2          5.000000      5.000000
  4    3          4.000000      4.000000
  5    3          5.000000      5.000000
  5    4          5.000000      5.000000
=================================================
```

---

## Citation

If you use this code, please cite both the manuscript and the software:

### Manuscript

```bibtex
@article{Moniz2025hypergraph,
  author  = {Moniz, Gabriel},
  title   = {{k-Uniform} complete hypergraph states stabilizers
             in terms of local operators},
  journal = {Physical Review A},
  year    = {2025},
  note    = {arXiv preprint}
}
```

### Software

```bibtex
@software{Moniz2025code,
  author    = {Moniz, Gabriel},
  title     = {hypergraph-functional-bounds},
  year      = {2025},
  url       = {https://github.com/GabrielMoniz/hypergraph-functional-bounds},
  license   = {MIT}
}
```

A machine-readable citation is also available in [`CITATION.cff`](CITATION.cff).

### Ket.jl

This code uses the [Ket.jl](https://github.com/araujoms/ket) package. Please also cite it if you use this work:

```bibtex
@misc{Araujo2024ket,
  author = {Araújo, Mateus},
  title  = {{Ket.jl}: a {Julia} package for quantum information},
  year   = {2024},
  url    = {https://github.com/araujoms/ket}
}
```

---

## License

This project is released under the [MIT License](LICENSE).
