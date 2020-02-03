include("force_layout_fire.jl")
import(JSON)

# test graph
chains = [["A", "B", "C", "D", "E"],
          ["R", "E", "T", "ER", "erz", "ZE"],
          ["C", "re", "eze"]]

# ====================
# Load graph from Json
files = filter(x -> endswith(x, "_chains.json"), readdir("./data/"))

file_name = Juno.selector(files)
path = string("./data/", file_name)
data = read(path, String)
chains = JSON.parse(data, dicttype=Array{Array{String}});
sort!(chains, by=x->length(x), rev=true)
println("number of chains: ", length(chains))

#chains = chains[1:3]
graph = create_graph(chains);
place_on_circle!(graph)
reset_forces!(graph)
println("number of chains: ", length(chains))


fire = Fire()
for k in 1:100
    step!(graph, fire)
    end;

graphplot(graph, chains)
