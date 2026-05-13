module Simulation

using ..GeneticArchitectures: OneLocusDiploid, n_genotypes
using ..SelectionModels: SelectionModel, select!
using ..MigrationModels: MigrationModel, migrate!
using ..MatingModels: MatingModel, mate!

export simulate, secondary_contact, allele_frequencies

"""
    secondary_contact(arch::OneLocusDiploid; n_demes=101) -> Matrix{Float64}

Return an initial genotype frequency matrix for the secondary contact scenario.

The transect is split into three zones:

- **Left half** (demes 1 to `(n_demes - 1) ÷ 2`): fixed A1A1 (`freq[1, j] = 1`)
- **Contact deme** (`(n_demes + 1) ÷ 2`): all heterozygotes (`freq[2, j] = 1`)
- **Right half** (remaining demes): fixed A2A2 (`freq[3, j] = 1`)

With the default `n_demes = 101`, the contact zone is at deme 51, the geometric
centre of the transect. An odd number of demes is recommended so the centre is
an integer index and the two flanking regions are equal in length.

The contact deme is set to all heterozygotes rather than a 50/50 mix of the two
homozygotes, following the Pascal reference implementation convention. A
heterozygote-only contact is the natural starting point for studying introgression
between two differentiated populations meeting for the first time.

# Examples
```jldoctest
julia> using HybridZones

julia> arch = OneLocusDiploid();

julia> state = secondary_contact(arch);

julia> size(state)
(3, 101)

julia> state[:, 1]
3-element Vector{Float64}:
 1.0
 0.0
 0.0

julia> state[:, 51]
3-element Vector{Float64}:
 0.0
 1.0
 0.0

julia> state[:, 101]
3-element Vector{Float64}:
 0.0
 0.0
 1.0
```
"""
function secondary_contact(arch::OneLocusDiploid; n_demes::Int = 101)
    n_demes > 0 || throw(ArgumentError("n_demes must be positive, got $n_demes"))
    state = zeros(Float64, n_genotypes(arch), n_demes)
    contact = (n_demes + 1) ÷ 2
    left_end = (n_demes - 1) ÷ 2
    for j in 1:left_end
        state[1, j] = 1.0
    end
    state[2, contact] = 1.0
    for j in (contact + 1):n_demes
        state[3, j] = 1.0
    end
    return state
end

"""
    allele_frequencies(state::AbstractMatrix, arch::OneLocusDiploid) -> Vector{Float64}

Compute the frequency of the A2 allele in each deme from a genotype frequency matrix.

For `OneLocusDiploid`, the A2 frequency at deme `j` is:

    p_A2[j] = freq[3, j] + 0.5 × freq[2, j]

which counts A2A2 homozygotes fully and heterozygotes at half weight.

# Arguments
- `state`: 3×n_demes genotype frequency matrix (rows: A1A1, A1A2, A2A2; columns: demes)
- `arch`: genetic architecture; currently only `OneLocusDiploid` is supported

# Examples
```jldoctest
julia> using HybridZones

julia> arch = OneLocusDiploid();

julia> state = reshape([0.25, 0.50, 0.25], 3, 1);

julia> allele_frequencies(state, arch)
1-element Vector{Float64}:
 0.5
```
"""
function allele_frequencies(state::AbstractMatrix, ::OneLocusDiploid)
    n_demes = size(state, 2)
    freq_A2 = Vector{Float64}(undef, n_demes)
    @inbounds for j in 1:n_demes
        freq_A2[j] = state[3, j] + 0.5 * state[2, j]
    end
    return freq_A2
end

"""
    allele_frequencies(trajectory, arch::OneLocusDiploid) -> Matrix{Float64}

Compute A2 allele frequencies across all time points in a simulation trajectory.

Applies [`allele_frequencies(state, arch)`](@ref) slice-by-slice to a 3D trajectory
array of shape `(n_genotypes, n_demes, n_time)`.

# Returns
A matrix of shape `(n_demes, n_time)` where entry `[j, t]` is the frequency of A2
at deme `j` at time point `t`. For a trajectory returned by `simulate` with
`n_generations` generations, `t = 1` is the initial state and `t = n_generations + 1`
is the final state.

# Arguments
- `trajectory`: 3D array of shape `(n_genotypes, n_demes, n_time)` from [`simulate`](@ref)
- `arch`: genetic architecture
"""
function allele_frequencies(
        trajectory::AbstractArray{T, 3}, arch::OneLocusDiploid
) where {T <: Real}
    n_demes = size(trajectory, 2)
    n_time = size(trajectory, 3)
    result = Matrix{Float64}(undef, n_demes, n_time)
    for t in 1:n_time
        result[:, t] = allele_frequencies(@view(trajectory[:, :, t]), arch)
    end
    return result
end

"""
    simulate(arch, mat, sel, mig, initial_state; n_generations, lifecycle_order=:mate_migrate_select, save_trajectory=true)

Run a forward simulation of a hybrid zone for `n_generations` generations.

Each generation applies mating, selection, and migration in the order specified
by `lifecycle_order`. The default `:mate_migrate_select` ordering — mating, then
migration, then selection — matches Mallet's Pascal reference implementation and
is the recommended default for warning-color hybrid zone models.

## Lifecycle orderings

- `:mate_migrate_select` (default): mating restores Hardy-Weinberg within demes,
  then migration redistributes frequencies across demes, then selection acts on
  local frequencies. This matches Jim Mallet's standard Pascal ordering.
- `:mate_select_migrate`: mating, then selection, then migration. An alternative
  ordering for organisms where dispersal occurs after selection.
- `:select_migrate_mate`: selection first, then migration, then mating at the
  end of the generation. Biologically valid for organisms where dispersal occurs
  late in the life cycle after selection has acted, with mating completing the
  generation.

Other values raise an `ArgumentError`.

## Breaking change from previous API

The previous signature was `simulate(arch, sel, mig, initial_state; ...)`.
The new signature inserts `mat::MatingModel` as the second positional argument:
`simulate(arch, mat, sel, mig, initial_state; ...)`. No backward-compatible
default method is provided — the behavior change (HW restoration each generation)
is significant enough to require explicit user choice.

## Buffer reuse

The inner loop preallocates working buffers and updates them in place, avoiding
per-generation heap allocation. For `save_trajectory = true`, two intermediate
scratch buffers are allocated alongside the trajectory array. For
`save_trajectory = false`, three buffers are rotated per operation each generation.

## Migration convention

Migration is applied row-by-row (one genotype at a time) by calling
`migrate!` for each genotype row. This keeps `migrate!` operating on
`AbstractVector` as documented, matching the per-genotype stepping-stone logic
in the Pascal reference implementation.

# Arguments
- `arch::OneLocusDiploid`: genetic architecture
- `mat::MatingModel`: mating model (e.g. `RandomMating()`)
- `sel::SelectionModel`: selection model
- `mig::MigrationModel`: migration model
- `initial_state::AbstractMatrix`: `n_genotypes × n_demes` genotype frequency matrix;
  each column must sum to approximately 1.0

# Keyword arguments
- `n_generations::Int`: number of generations to simulate; must be positive
- `lifecycle_order::Symbol`: `:mate_migrate_select` (default), `:mate_select_migrate`,
  or `:select_migrate_mate`
- `save_trajectory::Bool`: if `true` (default), return all generations; if `false`,
  return only the final state (more memory-efficient for long runs)

# Returns
- `save_trajectory = true`: `Array{Float64, 3}` of shape
  `(n_genotypes, n_demes, n_generations + 1)`. Slice `[:, :, 1]` is `initial_state`;
  slice `[:, :, g + 1]` is the state after generation `g`.
- `save_trajectory = false`: `Matrix{Float64}` of shape `(n_genotypes, n_demes)`
  with the state after `n_generations` generations.

# Examples
```jldoctest
julia> using HybridZones

julia> arch = OneLocusDiploid(dominance = :codominant);

julia> mat = RandomMating();

julia> sel = FrequencyDependentSelection(s = 0.3);

julia> mig = BinomialStepping(30.0);

julia> initial = secondary_contact(arch);

julia> trajectory = simulate(arch, mat, sel, mig, initial; n_generations = 10);

julia> size(trajectory)
(3, 101, 11)

julia> trajectory[:, :, 1] == initial
true

julia> final_cline = allele_frequencies(trajectory[:, :, end], arch);

julia> length(final_cline)
101
```
"""
function simulate(
        arch::OneLocusDiploid,
        mat::MatingModel,
        sel::SelectionModel,
        mig::MigrationModel,
        initial_state::AbstractMatrix;
        n_generations::Int,
        lifecycle_order::Symbol = :mate_migrate_select,
        save_trajectory::Bool = true
)
    n_generations > 0 ||
        throw(ArgumentError("n_generations must be positive, got $n_generations"))
    lifecycle_order ∈ (:mate_migrate_select, :mate_select_migrate, :select_migrate_mate) ||
        throw(ArgumentError(
            "lifecycle_order must be :mate_migrate_select, :mate_select_migrate, " *
            "or :select_migrate_mate, got :$lifecycle_order"
        ))
    ng = n_genotypes(arch)
    size(initial_state, 1) == ng ||
        throw(DimensionMismatch(
            "initial_state has $(size(initial_state, 1)) genotype rows; " *
            "expected $ng for $(typeof(arch))"
        ))
    n_demes = size(initial_state, 2)
    for j in 1:n_demes
        col_sum = sum(@view(initial_state[:, j]))
        isapprox(col_sum, 1.0; atol = 1e-8) ||
            throw(ArgumentError(
                "column $j of initial_state sums to $col_sum; expected ~1.0"
            ))
    end

    if save_trajectory
        trajectory = zeros(Float64, ng, n_demes, n_generations + 1)
        trajectory[:, :, 1] .= initial_state
        buf_a = zeros(Float64, ng, n_demes)
        buf_b = zeros(Float64, ng, n_demes)
        for g in 1:n_generations
            current = @view(trajectory[:, :, g])
            next = @view(trajectory[:, :, g + 1])
            _apply_lifecycle!(next, current, buf_a, buf_b, mat, sel, mig, arch,
                lifecycle_order, ng)
        end
        return trajectory
    else
        # Three rotating buffers: buf_a holds the current state, buf_b and buf_c
        # are intermediates. _apply_lifecycle! writes the result back into buf_a
        # (out === current is safe because each lifecycle step fully consumes its
        # input before writing its output to a separate buffer).
        buf_a = similar(initial_state, Float64)
        buf_a .= initial_state
        buf_b = similar(buf_a)
        buf_c = similar(buf_a)
        for _ in 1:n_generations
            _apply_lifecycle!(buf_a, buf_a, buf_b, buf_c, mat, sel, mig, arch,
                lifecycle_order, ng)
        end
        return buf_a
    end
end

# Apply one full lifecycle, writing the result into `out`.
# `current` is the read-only input state; `tmp1` and `tmp2` are scratch buffers.
# Safe to call with out === current: each step finishes reading its input
# before the next step writes to a different buffer.
function _apply_lifecycle!(
        out::AbstractMatrix,
        current::AbstractMatrix,
        tmp1::AbstractMatrix,
        tmp2::AbstractMatrix,
        mat::MatingModel,
        sel::SelectionModel,
        mig::MigrationModel,
        arch::OneLocusDiploid,
        lifecycle_order::Symbol,
        ng::Int
)
    if lifecycle_order === :mate_migrate_select
        mate!(tmp1, current, mat, arch)
        for i in 1:ng
            migrate!(@view(tmp2[i, :]), @view(tmp1[i, :]), mig)
        end
        select!(out, tmp2, sel, arch)
    elseif lifecycle_order === :mate_select_migrate
        mate!(tmp1, current, mat, arch)
        select!(tmp2, tmp1, sel, arch)
        for i in 1:ng
            migrate!(@view(out[i, :]), @view(tmp2[i, :]), mig)
        end
    else
        # :select_migrate_mate
        select!(tmp1, current, sel, arch)
        for i in 1:ng
            migrate!(@view(tmp2[i, :]), @view(tmp1[i, :]), mig)
        end
        mate!(out, tmp2, mat, arch)
    end
    return out
end

end # module Simulation
