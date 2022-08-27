### A Pluto.jl notebook ###
# v0.19.11

using Markdown
using InteractiveUtils

# ╔═╡ a97832bf-5d0b-40c5-860c-c319c22634d4
begin
	using Plots, LaTeXStrings
	
	plot_font = "Computer Modern"
	default(
		fontfamily=plot_font,
		linewidth=1.5, 
		framestyle=:box, 
		label=nothing,
		grid=true,
		dpi=300
	)
	
	Plots.scalefontsizes()
	Plots.scalefontsizes(1.3)
end

# ╔═╡ 336b315a-14cd-11ed-33ba-633ac3d588a8
# begin
# 	import Pkg
# 	# activate the shared project environment
# 	Pkg.activate(Base.current_project())
# 	# instantiate, i.e. make sure that all packages are downloaded
# 	Pkg.instantiate()
# end

# ╔═╡ a359d795-b97f-4472-905f-38f0489b691c
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

# ╔═╡ b604c715-7942-426e-8622-4b6562a79509
begin
	using Distributions, StatsBase
	
	@enum Gene A G F
	ν = 16 # ν fixed to 16 for now to avoid bit masking hassle
	Antibiotic = UInt16 

	struct Genome
		genes::Vector{Gene}
		antibiotics::Vector{Antibiotic}
		genes_locs::Dict{Gene, Set{Int}}
	end

	function mutate!(genome::Genome, params::Parameters)
		mutate_new_fragile_gene!(genome, params)

		if length(genome.genes) > 0
			mutate_duplicate_delete!(genome, params)
			mutate_fragile_site_deletion!(genome, params)
			mutate_antibiotic_genes!(genome, params)
		end
	end

	function mutate_duplicate_delete!(genome::Genome, params::Parameters)
		n_mutations = rand(Binomial(length(genome.genes), params.μ_d))
		if iszero(n_mutations)
			return
		end
		
		mutate_genes = sample(1:length(genome.genes), n_mutations, replace=false)
		mutate_idxs = zeros(Int, params.max_genome_size)
		for gene_idx in mutate_genes
			mutate_idxs[gene_idx] = gene_idx
		end

		deletes = rand(Float64, length(mutate_genes)) .<= 0.5

		for (i, delete) in enumerate(deletes)
			gene_idx = mutate_idxs[mutate_genes[i]]
			if delete
				# Deletion
				delete!(genome.genes_locs[genome.genes[gene_idx]], gene_idx)
				deleteat!(genome.genes, gene_idx)
				mutate_idxs[gene_idx:end] .-= 1
			else
				# Duplication
				if length(genome.genes) >= params.max_genome_size
					continue
				end
				duplicated_gene = genome.genes[gene_idx]
				target_locus = rand(1:length(genome.genes))
				
				insert!(genome.genes, target_locus, duplicated_gene)
				push!(genome.genes_locs[duplicated_gene], target_locus)
				mutate_idxs[target_locus:end] .+= 1
				
				if duplicated_gene === A
					# is this right? copy antibiotic type if an antibiotic
					#   gene is duplicated
					genome.antibiotics[target_locus] = genome.antibiotics[gene_idx]
				end
			end
		end
	end

	function mutate_new_fragile_gene!(genome::Genome, params::Parameters)
		if rand(Float64) <= params.μ_n && length(genome.genes) < params.max_genome_size
			target_locus =
				if length(genome.genes) == 0
					1
				else
					rand(1:length(genome.genes))
				end
			
			insert!(genome.genes, target_locus, F)
			push!(genome.genes_locs[F], target_locus)
		end
	end

	function mutate_fragile_site_deletion!(genome::Genome, params::Parameters)
		delete_sites = rand(Float64, length(genome.genes_locs[F])) .<= params.μ_f
		for (fragile_gene_idx, delete_site) in zip(genome.genes_locs[F], delete_sites)
			if !delete_site
				continue
			end

			# delete region offstream of the genomic location of the fragile site
			end_offstream_idx = fragile_gene_idx + 1
			while true
				at_end = end_offstream_idx >= length(genome.genes)
				if at_end
					end_offstream_idx = length(genome.genes)
					break
				elseif genome.genes[end_offstream_idx + 1] === F
					break
				end
				
				end_offstream_idx += 1
			end

			remove_idxs = collect(fragile_gene_idx:end_offstream_idx)
			remove_genes = genome.genes[remove_idxs]
			deleteat!(genome.genes, remove_idxs)
			
			for (removed_gene, gene_idx) in zip(remove_genes, remove_idxs)
				delete!(genome.genes_locs[removed_gene], gene_idx)
			end
		end
	end

	function mutate_antibiotic_genes!(genome::Genome, params::Parameters)
		n_bits = ν::Int64
		trigger_mutation = rand(Float64, length(genome.genes_locs[A])) .<= params.μ_a
		
		for (i, ab_gene_idx) in enumerate(genome.genes_locs[A])
			if trigger_mutation[i]
				genome.antibiotics[ab_gene_idx] ⊻= 1 << rand(0:(n_bits - 1))
			end
		end
	end
end

# ╔═╡ 3bac1164-975a-411d-a12e-bbe6958ebed5
function generate_random_genome(size::Int, max_genome_size::Int)::Genome
	genes_locs = Dict{Gene, Set{Int}}(A => Set([]), G => Set([]), F => Set([]))
	antibiotics = zeros(Antibiotic, max_genome_size)
	
	possible_genes = (A, G) # no fragile sites on initialization
	rands = rand(possible_genes, size)
	genes = Vector{Gene}(undef, size)

	for (i, r) in enumerate(rands)
		genes[i] = r
		push!(genes_locs[genes[i]], i)
		if genes[i] === A
			# generate random antibiotic genotype
			antibiotics[i] = rand(typemin(Antibiotic):typemax(Antibiotic))
		end
	end

	Genome(genes, antibiotics, genes_locs)
end

# ╔═╡ a213a07f-65c1-4636-b646-8b6503019513
function genome_from_sequence(genes::Vector{Gene}, max_genome_size::Int)::Genome
	genes_locs = Dict{Gene, Set{Int}}(A => Set([]), G => Set([]), F => Set([]))
	antibiotics = zeros(Antibiotic, params.max_genome_size)

	for (i, gene) in enumerate(genes)
		push!(genes_locs[genes[i]], i)
		if genes[i] === A
			# generate random antibiotic genotype
			antibiotics[i] = rand(typemin(Antibiotic):typemax(Antibiotic))
		end
	end
	Genome(genes, antibiotics, genes_locs)
end

# ╔═╡ 4739d59b-a538-4736-9716-7ba295b770ee
begin
	intrinsic_replication_rate(g::Int, α_g::Float64, h_g::Int)::Float64 = α_g*g/(g + h_g)
	ab_production_rate(a::Int, g::Int, α_a::Float64, h_a::Int, β_g::Float64)::Float64 = α_a*a/(a + h_a) * exp(-β_g*g)

	Position = Tuple{Int, Int}
	
	mutable struct Bacterium
		id::Int
		colony_id::Int
		
		genome::Genome
		pos::Position
	
		intrinsic_replication_rate::Float64
		ab_production_rate::Float64
	
		function Bacterium(id::Int, colony_id::Int, genome::Genome, pos::Position, params::Parameters)
			n_growth_genes = length(genome.genes_locs[G])
			n_ab_genes = length(genome.genes_locs[A])
			new(
				id, colony_id, genome, pos,
				intrinsic_replication_rate(n_growth_genes, params.α_g, params.h_g),
				ab_production_rate(n_ab_genes, n_growth_genes, params.α_a, params.h_a, params.β_g)
			)
		end
	end

	function mutate!(bacterium::Bacterium, params::Parameters)
		mutate!(bacterium.genome, params)
		n_growth_genes = length(bacterium.genome.genes_locs[G])
		n_ab_genes = length(bacterium.genome.genes_locs[A])
		
		bacterium.intrinsic_replication_rate = intrinsic_replication_rate(n_growth_genes, params.α_g, params.h_g)
		bacterium.ab_production_rate = ab_production_rate(n_ab_genes, n_growth_genes, params.α_a, params.h_a, params.β_g)
	end

	@inline function replication_rate(bacterium::Bacterium, R::AbstractFloat)::Float64
		bacterium.intrinsic_replication_rate*R
	end

	@inline function move!(bacterium::Bacterium, new_pos::Position)
		bacterium.pos = new_pos
	end
end

# ╔═╡ 55a6715c-1a24-45b0-ab38-6cd52c1edfa8
begin
	Grid2D = Matrix

	@inline function neighbours(grid::Grid2D, (i, j)::Position)::Set{Position}
		L = size(grid, 1)
		# only considering Moore for now
		Set([
			(mod1(i - 1, L), mod1(j + 1, L)),
			(i, mod1(j + 1, L)),
			(mod1(i + 1, L), mod1(j + 1, L)),
			(mod1(i + 1, L), j),
			(mod1(i + 1, L), mod1(j - 1, L)),
			(i, mod1(j - 1, L)),
			(mod1(i - 1, L), mod1(j - 1, L)),
			(mod1(i - 1, L), j)
		])
	end

	@inline function empty_neighbours(grid::Grid2D{Int}, pos::Position)::Set{Position}
		empty_neighbours(grid, pos, neighbours(grid, pos))
	end
	
	@inline function empty_neighbours(grid::Grid2D{Int}, pos::Position, neighs::Set{Position})::Set{Position}
		Set([n for n in neighs if grid[n...] === 0])
	end
end

# ╔═╡ 0192ed14-dd0f-40de-b830-8987313857c8
begin
	using Random
	
	mutable struct World
		bacterium_grid::Grid2D{Int}
		ab_grid::Grid2D{Set{Antibiotic}}
	
		bacteria::Dict{Int, Bacterium}
		population_size::Int
		n::Int # total individual counter

		ab_circle_sites::Set{Position}
	
		params::Parameters

		function euclid_distance(pos_a::Position, pos_b::Vector{Int})::Float64
			sqrt((pos_a[1] - pos_b[1])^2 + (pos_a[2] - pos_b[2])^2)
		end

		function generate_ab_circle_sites(r_a::Int)::Set{Position}
			center = (0, 0)
			Set([(i, j) for i in -r_a:r_a, j in -r_a:r_a
			     if euclid_distance(center, [i, j]) <= r_a])
		end

		function World(params::Parameters)
			L = params.L
			
			init_bacterium_grid = zeros(Int, L, L)
			init_ab_grid = [Set{Antibiotic}() for i in 1:L, j in 1:L]
			init_bacteria = Dict{Int, Bacterium}()
	
			for i in 1:params.n_init_spores
				genome = generate_random_genome(params.init_genome_length, params.max_genome_size)

				# Initialize bacterium on a site that is not within `init_min_dist`
				#   of any other bacteria
				pos = Vector{Int}()
				# brute force this for now
				while true
					pos = rand(1:L, 2)
					if init_bacterium_grid[pos...] !== 0
						continue
					end
					
					too_close = false
					for (_, bacteria) in init_bacteria
						if euclid_distance(bacteria.pos, pos) <= params.init_min_dist
							too_close = true
							break
						end
					end

					if !too_close
						break
					end
				end
				init_bacterium_grid[pos...] = i
				init_bacteria[i] = Bacterium(i, i, genome, (pos...,), params)
			end
			
			init_population_size = init_n = length(init_bacteria)
			ab_circle_sites = generate_ab_circle_sites(params.r_a)
	
			new(
				init_bacterium_grid, init_ab_grid, init_bacteria, init_population_size, init_n, ab_circle_sites, params
			)
		end

		function World(params::Parameters, prev_world::World)
			# Initialize world from previous growth cycle
			L = params.L
			
			init_bacterium_grid = zeros(Int, L, L)
			init_ab_grid = [Set{Antibiotic}() for i in 1:L, j in 1:L]

			# randomly choose ξ survivors
			n_survivors = params.ξ*length(prev_world.bacteria) |> round |> Int
			if n_survivors == 0
				error("extinction")
			end
			survivors = rand(deepcopy(prev_world.bacteria), n_survivors)

			init_bacteria = survivors
			for bacterium in init_bacteria
				init_bacterium_grid[bacterium.pos...] = bacterium.pos
			end

			init_n = prev_world.n # take over individual count from previous cycle
			init_population_size = length(init_bacteria)
			ab_circle_sites = prev_world.ab_circle_sites
	
			new(
				init_bacterium_grid, init_ab_grid, init_bacteria, init_population_size, init_n, ab_circle_sites, params
			)
		end
	end

	function produce_antibiotics!(world::World, bacterium::Bacterium)
		(b_i, b_j) = bacterium.pos
		deposit_at_sites = rand(Float64, length(world.ab_circle_sites)) .<= bacterium.ab_production_rate
		
		for ((i, j), deposit) in zip(world.ab_circle_sites, deposit_at_sites)
			if deposit
				# deposit uniformly random antibiotic genotype
				ab_gene_loc = rand(bacterium.genome.genes_locs[A])
				ab_gene = bacterium.genome.antibiotics[ab_gene_loc]
				push!(
					world.ab_grid[mod1(b_i + i, world.params.L) , mod1(b_j + j, world.params.L)],
					ab_gene
				)
			end
		end
	end

	function ab_susceptibility_score(bacterium::Bacterium, antibiotics::Set{Antibiotic})::Int
		genome = bacterium.genome
		if isempty(genome.genes_locs[A])
			return 999
		end
		bacterium_antibiotics = getindex.(Ref(genome.antibiotics), collect(genome.genes_locs[A]))

		dists = [
			count_ones.(site_ab .⊻ bacterium_antibiotics) for site_ab in antibiotics
		]

		sum(minimum.(dists))
	end

	function replicate!(world::World, site::Position, candidates::Set{Position}, params::Parameters)::Union{Bacterium, Nothing}
		# not sure if correct implementation
		candids_ordered = collect(candidates)
		replication_probs = zeros(Float64, length(candids_ordered))

		for (i, pos) in enumerate(candids_ordered)
			bacterium = world.bacteria[world.bacterium_grid[pos...]]
			site_antibiotics = world.ab_grid[bacterium.pos...]
			susceptibility = ab_susceptibility_score(bacterium, site_antibiotics)
			R = exp(-world.params.β_r*susceptibility^2)
			replication_probs[i] = replication_rate(bacterium, R)
		end

		total_growth = sum(replication_probs)
		rn = rand(Float64)*params.η
		if rn >= total_growth
			return
		end

		c = 1
		cum_prob = replication_probs[1]
		while cum_prob < rn
			c += 1
			cum_prob += replication_probs[c]
		end
		
		winner_pos = candids_ordered[c]
		winner_bacterium = world.bacteria[world.bacterium_grid[winner_pos...]]

		world.n += 1
	
		daughter = deepcopy(winner_bacterium)
		daughter.pos = site
		daughter.id = world.n
		
		mutate!(daughter, params)
		add_bacterium!(world, daughter)

		return daughter
	end

	function add_bacterium!(world::World, bacterium::Bacterium)
		world.bacterium_grid[bacterium.pos...] = bacterium.id
		world.bacteria[bacterium.id] = bacterium
		world.population_size += 1
	end

	function delete_bacterium!(world::World, bacterium::Bacterium)
		world.bacterium_grid[bacterium.pos...] = 0
		delete!(world.bacteria, bacterium.id)
		world.population_size -= 1
	end

	function move_bacterium!(world::World, bacterium::Bacterium, new_pos::Position)
		world.bacterium_grid[bacterium.pos...] = 0
		move!(bacterium, new_pos)
		world.bacterium_grid[new_pos...] = bacterium.id
	end

	function update!(world::World; n::Int=1, data_collector::Function=nothing)
		collect_data = !isnothing(data_collector)
		if collect_data
			data_collector(world, 0) # capture initial state (t = 0)
		end
		
		for t in 1:n
			order = world.bacteria |> keys |> collect |> shuffle
			rands = rand(Float64, 2, length(order))

			for (j, move_r, death_r) in zip(order, rands[1,:], rands[2,:])
				bacterium = world.bacteria[j]

				# antibiotics production
				produce_antibiotics!(world, bacterium)

				# replication
				free_adjacent = empty_neighbours(world.bacterium_grid, bacterium.pos)
				
				if !isempty(free_adjacent)
					replication_site = rand(free_adjacent)
					site_neighbours = neighbours(world.bacterium_grid, replication_site)
					site_neighbours_free = empty_neighbours(world.bacterium_grid, replication_site, site_neighbours)
					occupied_sites = setdiff(site_neighbours, site_neighbours_free)
						
					daughter = replicate!(world, replication_site, occupied_sites, world.params)
					if !isnothing(daughter)
						delete!(free_adjacent, daughter.pos)
					end
				end

				# movement
				if !isempty(free_adjacent)
					if move_r <= world.params.p_mov
						move_bacterium!(world, bacterium, rand(free_adjacent))
					end
				end

				# death
				site_antibiotics = world.ab_grid[bacterium.pos...]
				susceptibility = ab_susceptibility_score(bacterium, site_antibiotics)
				R = exp(-world.params.β_r*susceptibility^2)
				if death_r <= 1 - R
					delete_bacterium!(world, bacterium)
				end
			end
			
			if collect_data
				data_collector(world, t)
			end
		end
	end
end

# ╔═╡ bcbba418-c3b8-4e27-b647-5c7df4c8c1ec
begin
	using FileIO, JLD2
	
	mutable struct ResultsData
		bacterium_grids::Vector{Grid2D{Int}}
		ab_grids::Vector{Grid2D{Set{Antibiotic}}}
	
		bacteria_states::Vector{Dict{Int, Bacterium}}
		ts::Vector{Int}

		i::Int
	
		function ResultsData(n_snapshots::Int)
			new(
				Vector{Grid2D{Int}}(undef, n_snapshots),
				Vector{Grid2D{Set{Antibiotic}}}(undef, n_snapshots),
	
				Vector{Dict{Int, Bacterium}}(undef, n_snapshots),
				Vector{Int}(undef, n_snapshots),

				1
			)
		end
	end
	
	function collect_all!(rd::ResultsData, world::World, t::Int, t_snapshots::Set{Int})
		if t in t_snapshots
			rd.bacterium_grids[rd.i] = deepcopy(world.bacterium_grid)
			rd.ab_grids[rd.i] = deepcopy(world.ab_grid)

			rd.bacteria_states[rd.i] = deepcopy(world.bacteria)
			rd.ts[rd.i] = t

			rd.i += 1
		end
	end

	function persist(rd::ResultsData, params::Parameters, filename::String)
		jldsave(filename, data=rd, parameters=params)
	end
end

# ╔═╡ 2229542e-15fb-4f19-8149-23534cdbff1f
parameters = Parameters(
	L=150,
	η=8,
	τ_s=2500,
	ξ=0.001,
	
	n_init_spores=60,
	init_genome_length=10,
	init_min_dist=15,
	max_genome_size=1024,

	α_g=0.1,
	h_g=10,
	β_r=0.3,
	μ_d=0.001,
	μ_f=0.01,
	μ_n=0.01,
	μ_a=0.005,

	α_a=1.,
	h_a=3,
	r_a=3,
	β_g=1.,

	p_mov=0.01
)

# ╔═╡ 0809f6a7-129f-4ed3-adac-4ef0b7a6d71d
# using BenchmarkTools

# ╔═╡ 5f2309ed-e4bb-4524-b2fa-f3c637958b5c
# begin
# 	@benchmark update!(world, n=5) setup=(world=World(parameters))
# end

# ╔═╡ 60244234-02be-44be-944b-262080511aff
# cycle = 1
# n_cycles = 20

# t_snapshots = Set(0:100:parameters.τ_s)
# rd = ResultsData(length(t_snapshots))
# data_collector = (world::World, t::Int) -> collect_all!(rd, world, t, t_snapshots)

# Random.seed!(19911)
# prev_world = World(parameters)
# update!(prev_world, n=parameters.τ_s, data_collector=data_collector)
# persist(rd, parameters, "./cycles/$(cycle).jld2")

# for n in 2:n_cycles
# 	cycle += 1
# 	rd = ResultsData(length(t_snapshots))
# 	world = World(parameters, prev_world)
# 	GC.gc()

# 	update!(prev_world, n=parameters.τ_s, data_collector=data_collector)
# 	persist(rd, parameters, "./cycles/$(cycle).jld2")

# 	prev_world = world
# end

# ╔═╡ 5aa1261b-a782-42ac-96e9-a1459681c22f
begin
	data_file = "./sample_run.jld2"
	if !isfile(data_file)
		t_snapshots = Set(0:100:parameters.τ_s)
		rd = ResultsData(length(t_snapshots))
		data_collector = (world::World, t::Int) -> collect_all!(rd, world, t, t_snapshots)
	
		Random.seed!(83257325)
		world = World(parameters)
		@time update!(world, n=parameters.τ_s, data_collector=data_collector)
		persist(rd, parameters, "sample_run.jld2")
	else
		rd = load("sample_run.jld2")["data"]
	end
end

# ╔═╡ da9b2151-8836-47c1-a4b7-e51f599c6fdf
begin
	using Statistics

	gene_types = [A, G, F]
	avg_counts = Dict(
		gene_type => Vector{Float64}(undef, length(rd.ts))
		for gene_type in gene_types
	)
	stdev_counts = Dict(
		gene_type => Vector{Float64}(undef, length(rd.ts))
		for gene_type in gene_types
	)
	
	for (i, t) in enumerate(rd.ts)
		bacteria = rd.bacteria_states[i]

		counts = Dict(
			gene_type => Vector{Int}(undef, length(bacteria))
			for gene_type in gene_types
		)
	
		for (j, bacterium) in enumerate(values(bacteria))
			for gene_type in gene_types
				gene_count = length(bacterium.genome.genes_locs[gene_type])
				counts[gene_type][j] = gene_count
			end
		end
	
		for gene_type in gene_types
			avg_counts[gene_type][i] = mean(counts[gene_type])
			stdev_counts[gene_type][i] = std(counts[gene_type])
		end
		
	end

	freqs = Dict(
		gene_type => Dict{Int, Int}() for gene_type in gene_types
	)
	
	for bacterium in values(rd.bacteria_states[end])
		for gene_type in gene_types
			gene_count = length(bacterium.genome.genes_locs[gene_type])
			freqs[gene_type][gene_count] = get(freqs[gene_type], gene_count, 0) + 1
			end
	end
	
end

# ╔═╡ 3a92e564-d9a3-4c0d-9404-4a338f9f0c74
begin
	snapshots = 2:2:length(rd.bacterium_grids)
	bacteria_heatmaps = [
		heatmap(
			map(b -> begin
				if b !== 0
					rd.bacteria_states[i][b].colony_id
				else
					-2
				end
			end, rd.bacterium_grids[i]),
			color=[:white, palette(:rainbow, parameters.n_init_spores)...],
			legend=false, axis=nothing
		)
		for i in snapshots
	]
	plot(bacteria_heatmaps..., framestyle=:box, size=(1100, 750), layout=(3, 5),
	     title=[latexstring("\$t = $(t)\$") for t in rd.ts[snapshots]'], dpi=300,
	     plot_title="Bacterium colonies")
	# savefig("bacteria.pdf")
end

# ╔═╡ 97b48d2f-3dc0-45a0-bd88-8177b50ef5ec
begin
	ab_heatmaps = [
		heatmap(map(ab -> begin
			if length(ab) !== 0
				length(ab)
			else
				0
			end
		end, rd.ab_grids[i]), legend=false, axis=nothing)
		for i in snapshots
	]

	plot(ab_heatmaps..., framestyle=:box, size=(1100, 750), layout=(3, 5),
	     title=[latexstring("\$t = $(t)\$") for t in rd.ts[snapshots]'], dpi=300,
	     plot_title="Antibiotic frequency")
	# savefig("antibiotics.pdf")
end

# ╔═╡ e647befe-b586-4934-8df1-3c5a193bb5a1
begin
	bacteria_growth_rate_heatmaps = [
		heatmap(map(b -> begin
			if b !== 0
				rd.bacteria_states[i][b].intrinsic_replication_rate
			else
				0
			end
		end, rd.bacterium_grids[i]),
		color=[:black, palette(:gist_heat)...], legend=false, axis=nothing)
		for i in snapshots
	]
	plot(bacteria_growth_rate_heatmaps..., framestyle=:box, size=(1100, 750),
		  layout=(3, 5), title=[latexstring("\$t = $(t)\$") for t in rd.ts[snapshots]'], dpi=300, alpha=0.1, plot_title="Intrinsic replication rate")
	# savefig("replication_rates.pdf")
end

# ╔═╡ 1b3b126b-e1c2-4e02-a193-4f7855ab5898
begin
	bacteria_ab_rate_heatmaps = [
		heatmap(map(b -> begin
			if b !== 0
				rd.bacteria_states[i][b].ab_production_rate
			else
				0
			end
		end, rd.bacterium_grids[i]),
		color=[:black, palette(:gist_heat)...], legend=false, axis=nothing)
		for i in snapshots
	]
	plot(bacteria_ab_rate_heatmaps..., framestyle=:box, size=(1100, 750),
		  layout=(3, 5), title=[latexstring("\$t = $(t)\$") for t in rd.ts[snapshots]'], dpi=300, alpha=0.1, plot_title="AB production rate")
end

# ╔═╡ 48dabc12-6ef6-4057-a07f-1970041cbb40
begin
	plot(rd.ts, avg_counts[A], ribbon=stdev_counts[A], label="A")
	plot!(rd.ts, avg_counts[G], ribbon=stdev_counts[G], label="G")
	plot!(rd.ts, avg_counts[F], ribbon=stdev_counts[F], label="F",
		  xaxis=L"$t$", yaxis="Avg. number of genes")
	plot!(rd.ts, sum(values(avg_counts)), ribbon=sum(values(stdev_counts)), label="Total")
end

# ╔═╡ e1561751-ccd5-446d-9f41-2213d7e95db4
plot(rd.ts, length.(rd.bacteria_states), xaxis=L"$t$", yaxis="Population size")

# ╔═╡ 36ffaf44-1520-49b1-b0de-38fddfe4e180
begin
	scatter(keys(freqs[G]) |> collect, values(freqs[G]) |> collect, alpha=0.8, label="G")
	scatter!(keys(freqs[A]) |> collect, values(freqs[A]) |> collect, alpha=0.8, label="A")
	scatter!(keys(freqs[F]) |> collect, values(freqs[F]) |> collect, alpha=0.8, label="F", xaxis="Gene count", yaxis="Frequency")
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
FileIO = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
JLD2 = "033835bb-8acc-5ee8-8aae-3f567f8a3819"
LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
StatsBase = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"

[compat]
Distributions = "~0.25.66"
FileIO = "~1.15.0"
JLD2 = "~0.4.22"
LaTeXStrings = "~1.3.0"
Plots = "~1.31.5"
StatsBase = "~0.33.21"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.3"
manifest_format = "2.0"

[[deps.Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "195c5505521008abea5aee4f96930717958eac6f"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.4.0"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

[[deps.Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f641eb0a4f00c343bbc32346e1217b86f3ce9dad"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.1"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "80ca332f6dcb2508adba68f22f551adb2d00a624"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.15.3"

[[deps.ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "38f7a08f19d8810338d4f5085211c7dfa5d5bdd8"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.4"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "ded953804d019afa9a3f98981d99b33e3db7b6da"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.0"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "Random"]
git-tree-sha1 = "1fd869cc3875b57347f7027521f561cf46d1fcd8"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.19.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "SpecialFunctions", "Statistics", "TensorCore"]
git-tree-sha1 = "d08c20eef1f2cbc6e60fd3612ac4340b89fea322"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.9.9"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

[[deps.Compat]]
deps = ["Dates", "LinearAlgebra", "UUIDs"]
git-tree-sha1 = "924cdca592bc16f14d2f7006754a621735280b74"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.1.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[deps.Contour]]
git-tree-sha1 = "d05d9e7b7aedff4e5b51a029dced05cfb6125781"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.2"

[[deps.DataAPI]]
git-tree-sha1 = "fb5f5316dd3fd4c5e7c30a24d50643b73e37cd40"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.10.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "d1fff3a548102f48987a52a2e0d114fa97d730f0"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.13"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[deps.DensityInterface]]
deps = ["InverseFunctions", "Test"]
git-tree-sha1 = "80c3e8639e3353e5d2912fb3a1916b8455e2494b"
uuid = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
version = "0.4.0"

[[deps.Distributions]]
deps = ["ChainRulesCore", "DensityInterface", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SparseArrays", "SpecialFunctions", "Statistics", "StatsBase", "StatsFuns", "Test"]
git-tree-sha1 = "aafa0665e3db0d3d0890cdc8191ea03dc279b042"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.66"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "5158c2b41018c5f7eb1470d558127ac274eca0c9"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.1"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[deps.DualNumbers]]
deps = ["Calculus", "NaNMath", "SpecialFunctions"]
git-tree-sha1 = "5837a837389fccf076445fce071c8ddaea35a566"
uuid = "fa6b7ba4-c1ee-5f82-b5fc-ecf0adba8f74"
version = "0.6.8"

[[deps.EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3f3a2501fa7236e9b911e0f7a588c657e822bb6d"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.3+0"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bad72f730e9e91c08d9427d5e8db95478a3c323d"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.4.8+0"

[[deps.Extents]]
git-tree-sha1 = "5e1e4c53fa39afe63a7d356e30452249365fba99"
uuid = "411431e0-e8b7-467b-b5e0-f676ba4f2910"
version = "0.1.1"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "Pkg", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "ccd479984c7838684b3ac204b716c89955c76623"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.2+0"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "94f5101b96d2d968ace56f7f2db19d0a5f592e28"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.15.0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FillArrays]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "Statistics"]
git-tree-sha1 = "246621d23d1f43e3b9c368bf3b72b2331a27c286"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "0.13.2"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "87eb71354d8ec1a96d4a7636bd57a7347dde3ef9"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.10.4+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pkg", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "d972031d28c8c8d9d7b41a536ad7bb0c2579caca"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.8+0"

[[deps.GR]]
deps = ["Base64", "DelimitedFiles", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Printf", "Random", "RelocatableFolders", "Serialization", "Sockets", "Test", "UUIDs"]
git-tree-sha1 = "037a1ca47e8a5989cc07d19729567bb71bfabd0c"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.66.0"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Pkg", "Qt5Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "c8ab731c9127cd931c93221f65d6a1008dad7256"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.66.0+0"

[[deps.GeoInterface]]
deps = ["Extents"]
git-tree-sha1 = "fb28b5dc239d0174d7297310ef7b84a11804dfab"
uuid = "cf35fbd7-0cd7-5166-be24-54bfbe79505f"
version = "1.0.1"

[[deps.GeometryBasics]]
deps = ["EarCut_jll", "GeoInterface", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "a7a97895780dab1085a97769316aa348830dc991"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.3"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "a32d672ac2c967f3deb8a81d828afc739c838a06"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.68.3+2"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "Dates", "IniFile", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "ed47af35905b7cc8f1a522ca684b35a212269bd8"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.2.0"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.HypergeometricFunctions]]
deps = ["DualNumbers", "LinearAlgebra", "OpenLibm_jll", "SpecialFunctions", "Test"]
git-tree-sha1 = "709d864e3ed6e3545230601f94e11ebc65994641"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.11"

[[deps.IniFile]]
git-tree-sha1 = "f550e6e32074c939295eb5ea6de31849ac2c9625"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.1"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "b3364212fb5d870f724876ffcd34dd8ec6d98918"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.7"

[[deps.IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[deps.IterTools]]
git-tree-sha1 = "fa6287a4469f5e048d763df38279ee729fbd44e5"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.4.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLD2]]
deps = ["FileIO", "MacroTools", "Mmap", "OrderedCollections", "Pkg", "Printf", "Reexport", "TranscodingStreams", "UUIDs"]
git-tree-sha1 = "81b9477b49402b47fbe7f7ae0b252077f53e4a08"
uuid = "033835bb-8acc-5ee8-8aae-3f567f8a3819"
version = "0.4.22"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b53380851c6e6664204efb2e62cd24fa5c47e4ba"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.2+0"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "Printf", "Requires"]
git-tree-sha1 = "1a43be956d433b5d0321197150c2f94e16c0aaa0"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.16"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "7739f837d6447403596a75d19ed01fd08d6f56bf"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.3.0+3"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "42b62845d70a619f063a7da093d995ec8e15e778"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+1"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "3eb79b0ca5764d4799c06699573fd8f533259713"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.4.0+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "361c2b088575b07946508f135ac556751240091c"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.17"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "5d4d2d9904227b8bd66386c1138cf4d5ffa826bf"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "0.4.9"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "3d3e902b31198a27340d0bf00d6ac452866021cf"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.9"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "Random", "Sockets"]
git-tree-sha1 = "d9ab10da9de748859a7780338e1d6566993d1f25"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.3"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[deps.Measures]]
git-tree-sha1 = "e498ddeee6f9fdb4551ce855a46f54dbd900245f"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.1"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "bf210ce90b6c9eed32d25dbcae1ebc565df2687f"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.2"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "a7c3d1da1189a1c2fe843a3bfa04d18d20eb3211"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.1"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e60321e3f2616584ff98f0a4f18d98ae6f89bbb3"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.17+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.PCRE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b2a7af664e098055a7529ad1a900ded962bca488"
uuid = "2f80f16e-611a-54ab-bc61-aa92de5b98fc"
version = "8.44.0+0"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "cf494dca75a69712a72b80bc48f59dcf3dea63ec"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.16"

[[deps.Parsers]]
deps = ["Dates"]
git-tree-sha1 = "0044b23da09b5608b4ecacb4e5e6c6332f833a7e"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.3.2"

[[deps.Pixman_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b4f5d02549a10e20780a24fce72bea96b6329e29"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.40.1+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "8162b2f8547bc23876edd0c5181b27702ae58dce"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.0.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "9888e59493658e476d3073f1ce24348bdc086660"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.3.0"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "GeometryBasics", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun", "Unzip"]
git-tree-sha1 = "05873db92e703f134649d88b8a164f3b7acb4d73"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.31.5"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "47e5f437cc0e7ef2ce8406ce1e7e24d44915f88d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.3.0"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Qt5Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "xkbcommon_jll"]
git-tree-sha1 = "c6c0f690d0cc7caddb74cef7aa847b824a16b256"
uuid = "ea2cea3b-5b76-57ae-a6ef-0a8af62496e1"
version = "5.15.3+1"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "78aadffb3efd2155af139781b8a8df1ef279ea39"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.4.2"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RecipesBase]]
git-tree-sha1 = "6bf3f380ff52ce0832ddd3a2a7b9538ed1bcca7d"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.2.1"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "RecipesBase"]
git-tree-sha1 = "e7eac76a958f8664f2718508435d058168c7953d"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.6.3"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "22c5201127d7b243b9ee1de3b43c408879dff60f"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "0.3.0"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "bf3188feca147ce108c76ad82c2792c57abe7b1f"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.0"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "68db32dff12bb6127bac73c209881191bf0efbb7"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.3.0+0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "f94f779c94e58bf9ea243e77a37e16d9de9126bd"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.1.1"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "874e8867b33a00e784c8a7e4b60afe9e037b74e1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.1.0"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "b3363d7460f7d098ca0912c69b082f75625d7508"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.0.1"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "d75bda01f8c31ebb72df80a46c88b25d1c79c56d"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.1.7"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "StaticArraysCore", "Statistics"]
git-tree-sha1 = "23368a3313d12a2326ad0035f0db0c0966f438ef"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.5.2"

[[deps.StaticArraysCore]]
git-tree-sha1 = "66fe9eb253f910fe8cf161953880cfdaef01cdf0"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.0.1"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f9af7f195fb13589dd2e2d57fdb401717d2eb1f6"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.5.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "d1bf48bfcc554a3761a133fe3a9bb01488e06916"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.21"

[[deps.StatsFuns]]
deps = ["ChainRulesCore", "HypergeometricFunctions", "InverseFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "5783b877201a82fc0014cbf381e7e6eb130473a4"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.0.1"

[[deps.StructArrays]]
deps = ["Adapt", "DataAPI", "StaticArrays", "Tables"]
git-tree-sha1 = "ec47fb6069c57f1cee2f67541bf8f23415146de7"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.11"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits", "Test"]
git-tree-sha1 = "5ce79ce186cc678bbb5c5681ca3379d1ddae11a1"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.7.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "216b95ea110b5972db65aa90f88d8d89dcb8851c"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.6"

[[deps.URIs]]
git-tree-sha1 = "e59ecc5a41b000fa94423a578d29290c7266fc10"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.4.0"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unzip]]
git-tree-sha1 = "34db80951901073501137bdbc3d5a8e7bbd06670"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.1.2"

[[deps.Wayland_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "3e61f0b86f90dacb0bc0e73a0c5a83f6a8636e23"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.19.0+0"

[[deps.Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4528479aa01ee1b3b4cd0e6faef0e04cf16466da"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.25.0+0"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "58443b63fb7e465a8a7210828c91c08b92132dff"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.9.14+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "5be649d550f3f4b95308bf0183b82e2582876527"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.6.9+4"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4e490d5c960c314f33885790ed410ff3a94ce67e"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.9+4"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fe47bd2247248125c428978740e18a681372dd4"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.3+4"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6783737e45d3c59a4a4c4091f5f88cdcf0908cbb"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.0+3"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "daf17f441228e7a3833846cd048892861cff16d6"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.13.0+3"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "926af861744212db0eb001d9e40b5d16292080b2"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.0+4"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "4bcbf660f6c2e714f87e960a171b119d06ee163b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.2+4"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "5c8424f8a67c3f2209646d4425f3d415fee5931d"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.27.0+4"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "79c31e7844f6ecf779705fbc12146eb190b7d845"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.4.0+3"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e45044cd873ded54b6a5bac0eb5c971392cf1927"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.2+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3a2ea60308f0996d26f1e5354e10c24e9ef905d4"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.4.0+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "ece2350174195bb31de1a63bea3a41ae1aa593b6"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "0.9.1+5"
"""

# ╔═╡ Cell order:
# ╠═336b315a-14cd-11ed-33ba-633ac3d588a8
# ╠═a359d795-b97f-4472-905f-38f0489b691c
# ╠═b604c715-7942-426e-8622-4b6562a79509
# ╠═3bac1164-975a-411d-a12e-bbe6958ebed5
# ╠═a213a07f-65c1-4636-b646-8b6503019513
# ╠═4739d59b-a538-4736-9716-7ba295b770ee
# ╠═55a6715c-1a24-45b0-ab38-6cd52c1edfa8
# ╠═0192ed14-dd0f-40de-b830-8987313857c8
# ╠═bcbba418-c3b8-4e27-b647-5c7df4c8c1ec
# ╠═2229542e-15fb-4f19-8149-23534cdbff1f
# ╠═0809f6a7-129f-4ed3-adac-4ef0b7a6d71d
# ╠═5f2309ed-e4bb-4524-b2fa-f3c637958b5c
# ╠═60244234-02be-44be-944b-262080511aff
# ╠═5aa1261b-a782-42ac-96e9-a1459681c22f
# ╠═a97832bf-5d0b-40c5-860c-c319c22634d4
# ╠═3a92e564-d9a3-4c0d-9404-4a338f9f0c74
# ╠═97b48d2f-3dc0-45a0-bd88-8177b50ef5ec
# ╠═e647befe-b586-4934-8df1-3c5a193bb5a1
# ╠═1b3b126b-e1c2-4e02-a193-4f7855ab5898
# ╠═da9b2151-8836-47c1-a4b7-e51f599c6fdf
# ╠═48dabc12-6ef6-4057-a07f-1970041cbb40
# ╠═e1561751-ccd5-446d-9f41-2213d7e95db4
# ╠═36ffaf44-1520-49b1-b0de-38fddfe4e180
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
