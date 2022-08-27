export Grid2D, World, update!

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
            init_bacterium_grid, init_ab_grid, init_bacteria, init_population_size,
            init_n, ab_circle_sites, params
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
        survivors = deepcopy(rand(collect(values(prev_world.bacteria)), n_survivors))

        init_bacteria = Dict(bacterium.id => bacterium for bacterium in survivors)
        for bacterium in values(init_bacteria)
            init_bacterium_grid[bacterium.pos...] = bacterium.id
        end

        init_n = prev_world.n # take over individual count from previous cycle
        init_population_size = length(init_bacteria)
        ab_circle_sites = prev_world.ab_circle_sites

        new(
            init_bacterium_grid, init_ab_grid, init_bacteria, init_population_size,
            init_n, ab_circle_sites, params
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
