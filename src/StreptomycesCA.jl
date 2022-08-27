module StreptomycesCA

using Random, Distributions, StatsBase

export Parameters

Base.@kwdef struct Parameters
	L::Int     # Size of square lattice sides
	η::Int     # Neighbourhood size
	τ_s::Int   # Growth cycle duration
	ξ::Float64 # Fraction of spores to seed a growth cycle

	n_init_spores::Int       # Number of initial bacteria
	init_genome_length::Int  # Size of initial genome
	init_min_dist::Int       # Minimum Euclidean distance between initial bacteria
	max_genome_size::Int     # Maximum genome size

	# Replication parameters
	α_g::Float64 # Max replication probability per time unit
	h_g::Int     # Number of growth genes for half max growth rate
	β_r::Float64 # Antibiotic resistance factor
	μ_d::Float64 # Duplication/deletion probability per gene
	μ_f::Float64 # Fragile site-induced deletion probability (per fragile site)
	μ_n::Float64 # New fragile site formation probability (per genome)
	μ_a::Float64 # Antibiotic type mutation probability (per antibiotic gene)

	# Antibiotics parameters
	α_a::Float64 # Max antibiotic production probability per unit time
	h_a::Int     # Number of antibiotic genes for half max production rate
	r_a::Int     # Max distance of antibiotic placement
	β_g::Float64 # Antibiotic production decrease due to trade-off

	# Other
	p_mov::Float64 # Probability of migration into empty adjacent site
end

include("gene.jl")
include("bacterium.jl")
include("world.jl")

end # module
