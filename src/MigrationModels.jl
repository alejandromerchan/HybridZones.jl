module MigrationModels

using Distributions: Binomial, pdf

export MigrationModel, BinomialStepping
export migrate!, migration_variance, max_distance

"""
    MigrationModel

Abstract supertype for migration model specifications. Concrete subtypes carry
the migration kernel and parameters needed to redistribute frequencies across
demes via `migrate!`.
"""
abstract type MigrationModel end

"""
    BinomialStepping <: MigrationModel

Binomial stepping-stone migration model as used by Mallet & Barton (1989).

The migration kernel is a Binomial(2n, α) PMF evaluated at displacements
-n:n from each source deme, where `n` is the maximum migration distance
derived from the migration variance `σ²` and symmetry parameter `α`:

    n = ceil(Int, σ² / (2α(1-α)))

The kernel is stored in full: `kernel[n+1+i]` holds the probability of
displacement `i` for `i ∈ -n:n`, so the kernel has length `2n+1` and sums
to 1.0.

Boundary demes receive the frequency of the nearest in-bounds deme to fill
the kernel's reach beyond the transect edge (reflective boundary, matching
the Pascal source's `Reflection` procedure).

!!! note "Asymmetric case not yet implemented"
    `α = 0.5` gives a symmetric kernel. The constructor accepts other `α`
    values (which shift the kernel mean), but the `migrate!` implementation
    uses symmetric indexing (`kernel[n+1+i]`) that is correct only when the
    kernel is symmetric. Asymmetric support is deferred to a future PR.

# Fields
- `σ²::Float64`: migration variance per generation
- `α::Float64`: symmetry parameter (0.5 = symmetric)
- `n::Int`: maximum migration distance (half-width of kernel)
- `kernel::Vector{Float64}`: migration weights of length `2n+1`

# Examples
```jldoctest
julia> using HybridZones

julia> m = BinomialStepping(10.0)
BinomialStepping(σ²=10.0, α=0.5, n=20)

julia> m isa MigrationModel
true

julia> migration_variance(m)
10.0

julia> max_distance(m)
20

julia> isapprox(sum(HybridZones.MigrationModels.kernel(m)), 1.0; atol=1e-12)
true
```
"""
struct BinomialStepping <: MigrationModel
    σ²::Float64
    α::Float64
    n::Int
    kernel::Vector{Float64}
end

function BinomialStepping(σ²::Real; α::Real = 0.5)
    σ² > 0 || throw(ArgumentError("σ² must be positive, got $σ²"))
    0 < α < 1 || throw(ArgumentError("α must be in (0, 1), got $α"))
    n = ceil(Int, σ² / (2 * α * (1 - α)))
    dist = Binomial(2n, α)
    kernel = [pdf(dist, n + i) for i in (-n):n]
    sum(kernel) > 0 || throw(
        ArgumentError("migration kernel is degenerate (all-zero) for σ²=$σ², α=$α, n=$n"),
    )
    return BinomialStepping(Float64(σ²), Float64(α), n, kernel)
end

function Base.show(io::IO, m::BinomialStepping)
    return print(io, "BinomialStepping(σ²=$(m.σ²), α=$(m.α), n=$(m.n))")
end

"""
    migration_variance(m::BinomialStepping) -> Float64

Return the migration variance per generation `σ²`.
"""
migration_variance(m::BinomialStepping) = m.σ²

"""
    kernel(m::BinomialStepping) -> Vector{Float64}

Return the full migration weight vector of length `2n+1`. Element `n+1+i`
is the probability of displacement `i` for `i ∈ -n:n`.
"""
kernel(m::BinomialStepping) = m.kernel

"""
    max_distance(m::BinomialStepping) -> Int

Return `n`, the maximum displacement (in demes) that the kernel reaches.
"""
max_distance(m::BinomialStepping) = m.n

"""
    migrate!(freq_next, freq_current, m::BinomialStepping)

Apply one generation of binomial stepping-stone migration in-place.

Writes post-migration frequencies into `freq_next` from `freq_current`. Both
vectors must have the same length (number of demes). Demes beyond the transect
boundaries are treated as having the frequency of the nearest edge deme
(reflective boundary).

# Arguments
- `freq_next::AbstractVector`: output buffer; overwritten on return
- `freq_current::AbstractVector`: input frequencies, one value per deme
- `m::BinomialStepping`: migration model supplying kernel and range

# Examples
```jldoctest
julia> using HybridZones

julia> m = BinomialStepping(10.0);

julia> freq = fill(0.3, 100);

julia> freq_next = similar(freq);

julia> migrate!(freq_next, freq, m);

julia> all(isapprox.(freq_next, 0.3; atol=1e-10))
true
```
"""
function migrate!(
        freq_next::AbstractVector,
        freq_current::AbstractVector,
        m::BinomialStepping
)
    length(freq_next) == length(freq_current) ||
        throw(DimensionMismatch("freq_next and freq_current must have equal length"))
    n_demes = length(freq_current)
    n = m.n
    @inbounds for j in 1:n_demes
        s = 0.0
        for i in (-n):n
            source = j + i
            w = m.kernel[n + 1 + i]
            if source < 1
                s += w * freq_current[1]
            elseif source > n_demes
                s += w * freq_current[n_demes]
            else
                s += w * freq_current[source]
            end
        end
        freq_next[j] = s
    end
    return freq_next
end

end # module MigrationModels
