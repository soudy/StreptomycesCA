### A Pluto.jl notebook ###
# v0.19.11

using Markdown
using InteractiveUtils

# ╔═╡ f6c5ef58-2636-11ed-39bf-4f30092caea2
begin
	import Pkg
	# activate the shared project environment
	Pkg.activate(Base.current_project())
	# instantiate, i.e. make sure that all packages are downloaded
	Pkg.instantiate()
end

# ╔═╡ db4a793c-b0e0-47d1-9043-1bc039a005fe
using StreptomycesCA, JLD2

# ╔═╡ c8dd482c-6a96-49f1-b18c-e8333f53ab3f
using Plots, LaTeXStrings

# ╔═╡ 65b80d2c-a5de-4ac0-a6d1-9b37828de7b1
rd = load("../data/run1/20.jld2")["data"]

# ╔═╡ 39235129-2384-4748-8f88-414b44d7a13a
begin
	bacteria_heatmaps = [
		heatmap(
			map(b -> begin
				if b !== 0
					rd.bacteria_states[i][b].colony_id
				else
					-2
				end
			end, rd.bacterium_grids[i]),
			color=[:white, palette(:rainbow, 30)...],
			legend=false, axis=nothing
		)
		for i in length(rd.bacterium_grids)
	]
	plot(bacteria_heatmaps..., framestyle=:box, size=(500, 500), layout=(3, 5),
	     dpi=300)
	# savefig("bacteria.pdf")
end

# ╔═╡ d05589d0-e36b-44c0-a182-f9f61bece1f1
begin
	ab_heatmaps = [
		heatmap(map(ab -> begin
			if length(ab) !== 0
				length(ab)
			else
				0
			end
		end, i), legend=false, axis=nothing)
		for i in [rd.ab_grids[1]]
	]

	plot(ab_heatmaps..., framestyle=:box, size=(500, 500), layout=(3, 5), dpi=300)
	# savefig("antibiotics.pdf")
end

# ╔═╡ Cell order:
# ╠═f6c5ef58-2636-11ed-39bf-4f30092caea2
# ╠═db4a793c-b0e0-47d1-9043-1bc039a005fe
# ╠═65b80d2c-a5de-4ac0-a6d1-9b37828de7b1
# ╠═c8dd482c-6a96-49f1-b18c-e8333f53ab3f
# ╠═39235129-2384-4748-8f88-414b44d7a13a
# ╠═d05589d0-e36b-44c0-a182-f9f61bece1f1
