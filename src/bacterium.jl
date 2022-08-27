export Position, Bacterium

function intrinsic_replication_rate(g::Int, α_g::Float64, h_g::Int)::Float64
    return α_g*g/(g + h_g)
end

function ab_production_rate(a::Int, g::Int, α_a::Float64, h_a::Int, β_g::Float64)::Float64
    return α_a*a/(a + h_a) * exp(-β_g*g)
end

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
    return bacterium.intrinsic_replication_rate*R
end

@inline function move!(bacterium::Bacterium, new_pos::Position)
    bacterium.pos = new_pos
end
