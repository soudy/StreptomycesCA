using Random, StreptomycesCA, FileIO, JLD2

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

parameters = Parameters(
    L=150,
    η=8,
    τ_s=2500,
    ξ=0.001,

    n_init_spores=30,
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

    α_a=1.0,
    h_a=3,
    r_a=3,
    β_g=1.0,

    p_mov=0.01
)

cycle = 1
n_cycles = 20

t_snapshots = Set([parameters.τ_s])
rd = ResultsData(length(t_snapshots))
data_collector(world::World, t::Int) = collect_all!(rd, world, t, t_snapshots)

Random.seed!(19911)
prev_world = World(parameters)

println("Cycle $(cycle)/$(n_cycles)")
tt = @timed update!(prev_world, n=parameters.τ_s, data_collector=data_collector)
println("  Took $(tt[:time]), persisting...")
persist(rd, parameters, "./data/run1/$(cycle).jld2")

for n in 2:n_cycles
    global cycle += 1
    global rd = ResultsData(length(t_snapshots))
    world = World(parameters, prev_world)
    data_collector(world::World, t::Int) = collect_all!(rd, world, t, t_snapshots)
    GC.gc()

    println("Cycle $(cycle)/$(n_cycles)")

    global tt = @timed update!(world, n=parameters.τ_s, data_collector=data_collector)
    println("  Took $(tt[:time]), persisting...")
    persist(rd, parameters, "./data/run1/$(cycle).jld2")

    global prev_world = world
end
