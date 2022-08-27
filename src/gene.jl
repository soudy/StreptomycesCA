export Gene, Antibiotic, Genome, generate_random_genome, genome_from_sequence

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
