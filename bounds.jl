# bounds.jl
#
# Computes the local (classical) and quantum bounds of the Bell functional
# associated with k-uniform complete hypergraph states on n qubits, as
# constructed in:
#
#   G. Moniz et al., "k-Uniform complete hypergraph states stabilizers in
#   terms of local operators" (arXiv preprint, 2025).
#
# The k-uniform complete hypergraph state |G_{n,k}⟩ is the n-qubit state
# defined by
#
#   |G_{n,k}⟩ = ∏_{e ∈ E} CZ_e |+⟩^⊗n ,
#
# where E is the set of all k-element subsets of {1,…,n} and CZ_e is the
# generalized controlled-Z gate on the qubits in e.  Its stabilizer
# generators are
#
#   g_j = X_j  ∏_{e ∈ E : j ∈ e}  CZ_{e \ {j}} ,     j = 1,…,n ,
#
# where each CZ_{e \ {j}} is a (k−1)-qubit gate.
#
# The Bell functional studied in the manuscript is defined through the
# Pauli expansion of these stabilizer generators.  Every generator g_j is
# expressed as a linear combination of tensor products of local Pauli
# operators; the resulting multi-linear Bell expression is evaluated both
# for quantum strategies (Tsirelson / quantum bound) and for local hidden-
# variable strategies (local / classical bound).
#
# Dependencies: Ket.jl, Combinatorics.jl, LinearAlgebra (stdlib).
#
# Usage:
#   julia bounds.jl
# or from the REPL:
#   include("bounds.jl")

using LinearAlgebra
using Combinatorics
using Printf
using Ket

# ─── Stabilizer-operator construction ─────────────────────────────────────────

"""
    stabilizer_matrix(j, n, k) → Matrix{ComplexF64}

Return the 2ⁿ × 2ⁿ matrix of the j-th stabilizer generator of the
k-uniform complete hypergraph state on n qubits:

    g_j = X_j ⊗ ∏_{e ∈ E : j ∈ e} CZ_{e \\ {j}} .

Qubits are ordered from most-significant (qubit 1) to least-significant
(qubit n) in the computational basis.
"""
function stabilizer_matrix(j::Int, n::Int, k::Int)
    dim     = 2^n
    others  = setdiff(1:n, [j])
    # (k−1)-element sub-edges: edges e\{j} for each e containing j
    sub_edges = collect(combinations(others, k - 1))

    mat = zeros(ComplexF64, dim, dim)
    for col in 0:(dim - 1)
        # Extract qubit values (MSB first: qubit i ↔ bit at position n-i)
        bits = [(col >> (n - i)) & 1 for i in 1:n]

        # Phase contributed by each CZ sub-edge: (-1) when all qubits are |1⟩
        phase = 1
        for edge in sub_edges
            if all(bits[v] == 1 for v in edge)
                phase *= -1
            end
        end

        # X_j flips qubit j
        new_bits    = copy(bits)
        new_bits[j] = 1 - bits[j]
        row = sum(new_bits[i] << (n - i) for i in 1:n)

        mat[row + 1, col + 1] += phase
    end
    return mat
end

"""
    bell_operator(n, k) → Hermitian{ComplexF64}

Return the Bell operator  B = ∑_j g_j  as a 2ⁿ × 2ⁿ Hermitian matrix.
"""
function bell_operator(n::Int, k::Int)
    dim = 2^n
    B   = zeros(ComplexF64, dim, dim)
    for j in 1:n
        B .+= stabilizer_matrix(j, n, k)
    end
    return Hermitian(B)
end

"""
    hypergraph_state(n, k) → Vector{ComplexF64}

Return the state vector of the k-uniform complete hypergraph state:

    |G_{n,k}⟩ = ∏_{e ∈ E} CZ_e |+⟩^⊗n .
"""
function hypergraph_state(n::Int, k::Int)
    edges = collect(combinations(1:n, k))
    dim   = 2^n
    ψ     = ones(ComplexF64, dim) / sqrt(dim)   # |+⟩^⊗n

    for edge in edges
        for state in 0:(dim - 1)
            # CZ_e applies (−1) when every qubit in e is |1⟩
            if all((state >> (n - i)) & 1 == 1 for i in edge)
                ψ[state + 1] *= -1
            end
        end
    end
    return ψ
end

# ─── Pauli expansion of the Bell functional ───────────────────────────────────

"""
    bell_functional_tensor(n, k) → Array{Float64, n}

Return the Bell functional in the correlation-tensor notation expected by
Ket.jl's `local_bound` and `tsirelson_bound`.

The array `FC` has shape `(3, 3, …, 3)` (n copies, one per party).  Index
mapping for party i:

  * 1  →  identity (I), used for marginal terms
  * 2  →  first dichotomic observable  = Z
  * 3  →  second dichotomic observable = X

So `FC[x₁,…,xₙ]` is the coefficient of the correlator
`⟨A₁^{x₁−1} ⊗ … ⊗ Aₙ^{xₙ−1}⟩` where `Aᵢ⁰ = I`, `Aᵢ¹ = Z`, `Aᵢ² = X`.

The functional encodes  β = ∑_j ⟨g_j⟩  using the Walsh–Pauli expansion

    g_j = ∑_{S ⊆ V\\{j}}  â_{j,S}  X_j ⊗ (⊗_{i∈S} Z_i) ⊗ I_rest ,

where the coefficients â_{j,S} are computed by discrete Fourier transform
over the (k−1)-fold CZ phases.
"""
function bell_functional_tensor(n::Int, k::Int)
    FC = zeros(fill(3, n)...)

    for j in 1:n
        others    = setdiff(1:n, [j])
        m         = n - 1                         # |others|
        sub_edges = collect(combinations(others, k - 1))

        for S in powerset(others)
            # Walsh–Pauli coefficient â_{j,S}
            c = 0.0
            for bitval in 0:(2^m - 1)
                # bit l encodes the qubit value of others[l]
                bits  = [(bitval >> (l - 1)) & 1 for l in 1:m]
                phase = 1
                for edge in sub_edges
                    if all(bits[findfirst(==(v), others)] == 1 for v in edge)
                        phase *= -1
                    end
                end
                walsh = prod((-1)^bits[findfirst(==(v), others)] for v in S; init = 1)
                c    += phase * walsh
            end
            c /= 2^m
            abs(c) > 1e-10 || continue

            # Build the index into FC:
            #   party j → X (index 3); parties in S → Z (index 2); rest → I (index 1)
            idx       = ones(Int, n)
            idx[j]    = 3
            for i in S
                idx[i] = 2
            end
            FC[idx...] += c
        end
    end
    return FC
end

# ─── Bound computation ────────────────────────────────────────────────────────

"""
    quantum_bound_exact(n, k) → Float64

Exact quantum bound: largest eigenvalue of the Bell operator B = ∑_j g_j.
"""
function quantum_bound_exact(n::Int, k::Int)
    B = bell_operator(n, k)
    return maximum(real(eigvals(B)))
end

"""
    local_bound_lhv(n, k) → Float64

Local (classical) bound of the Bell functional, computed exactly by Ket.jl's
branch-and-bound algorithm over all deterministic local strategies.
"""
function local_bound_lhv(n::Int, k::Int)
    FC = bell_functional_tensor(n, k)
    return local_bound(FC; correlation = true, marg = true)
end

"""
    quantum_bound_npa(n, k; level = 1) → Float64

Upper bound on the quantum value via the NPA (Navascués–Pironio–Acín)
semidefinite-programming hierarchy at the requested `level`, computed by
Ket.jl's `tsirelson_bound`.

Note: this function is provided as an independent cross-check for the
two-party case (n = 2).  For n ≥ 3 parties the NPA level needed to
capture all n-party correlators grows with n; use higher levels (e.g.
level = n) or the exact eigenvalue bound `quantum_bound_exact` instead.
"""
function quantum_bound_npa(n::Int, k::Int; level::Int = 1)
    FC    = bell_functional_tensor(n, k)
    qb, _ = tsirelson_bound(FC, level; verbose = false)
    return Float64(qb)
end

# ─── Main ─────────────────────────────────────────────────────────────────────

function main()
    println()
    println("k-Uniform Complete Hypergraph Bell Functional Bounds")
    println("=" ^ 49)
    println("  n    k   Quantum (exact)   Local (LHV)")
    println("-" ^ 49)

    cases = [(3, 2), (4, 2), (5, 2), (4, 3), (5, 3), (5, 4)]
    for (n, k) in cases
        n >= k >= 2 || continue

        qb_exact = quantum_bound_exact(n, k)
        lb       = local_bound_lhv(n, k)

        @printf("  %-3d  %-3d  %14.6f  %12.6f\n", n, k, qb_exact, lb)
    end
    println("=" ^ 49)
    println()
    println("Notes:")
    println("  * Quantum (exact) = largest eigenvalue of B = ∑_j g_j.")
    println("  * Local (LHV)     = classical maximum over deterministic")
    println("                      local strategies (Ket.jl branch-and-bound).")
end

main()
