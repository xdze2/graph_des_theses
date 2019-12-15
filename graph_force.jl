
using(Optim)


struct Graph
    chains
    nodes
end



## Force Layout 
# Spring
function attractive_gradient(x1, y1, x2, y2)
    u = x2 - x1
    v = y2 - y1
    return (-u, -v, u, v)
    end;
function attractive_energy(x1, y1, x2, y2)
    u = x2 - x1
    v = y2 - y1
    return (u^2 + v^2)/2
    end;

# Elec
const EPSILON = 1e-4
function repulsive_energy(x1, y1, x2, y2)
    u = x2 - x1
    v = y2 - y1
    return 1/sqrt(u^2 + v^2 + EPSILON)
    end;

function repulsive_gradient(x1, y1, x2, y2)
    u = x2 - x1
    v = y2 - y1
    d = sqrt(u^2 + v^2 + EPSILON)^3
    return (u/d, v/d, -u/d, -v/d)
    end;

# d/dx(-1/sqrt(x^2 + y^2)) = x/(x^2 + y^2)^(3/2)

function graph_energy(coords, graph)
    E = 0
    # Attractive
    for chain in graph.chains
        for (n1, n2) in zip(chain[1:end-1], chain[2:end])
            i1, i2 = graph.nodes[n1], graph.nodes[n2]
            x1, y1 = coords[1, i1], coords[2, i1]
            x2, y2 = coords[1, i2], coords[2, i2]
            E += attractive_energy(x1, y1, x2, y2)
            end;
        end;
    
    # Repulsive
    # Todo: insert Barnes Hut here...
    for i in 1:length(graph.nodes)
        x1, y1 = coords[1, i], coords[2, i]
        for j in i+1:length(graph.nodes)
            x2, y2 = coords[1, j], coords[2, j]
            E += repulsive_energy(x1, y1, x2, y2)
            end;
        end;
    return E
    end;

function graph_gradient!(G, coords, graph)
    fill!(G, 0.f0)
    # Attractive
    for chain in graph.chains
        for (n1, n2) in zip(chain[1:end-1], chain[2:end])
            i1, i2 = graph.nodes[n1], graph.nodes[n2]
            x1, y1 = coords[1, i1], coords[2, i1]
            x2, y2 = coords[1, i2], coords[2, i2]
            grad = attractive_gradient(x1, y1, x2, y2)
            G[1, i1] += grad[1]
            G[2, i1] += grad[2]
            G[1, i2] += grad[3]
            G[2, i2] += grad[4]
            end;
        end;
    
    # Repulsive
    # Todo: insert Barnes Hut here...
    for i in 1:length(graph.nodes)
        x1, y1 = coords[1, i], coords[2, i]
        for j in i+1:length(graph.nodes)
            x2, y2 = coords[1, j], coords[2, j]
            grad = repulsive_gradient(x1, y1, x2, y2)
            G[1, i] += grad[1]
            G[2, i] += grad[2]
            G[1, j] += grad[3]
            G[2, j] += grad[4]
            end;
        end;
    end;



function run_optimiz(graph)

    # Convert dict of coords to an array + dict of indexes
    nodes_idx = Dict(n => idx  for (idx, n) in enumerate(keys(graph.nodes)) );
    coords = zeros(Float32, 2, length(nodes_idx))
    for (n, idx) in pairs(nodes_idx)
        coords[:, idx] = graph.nodes[n]
        end;

    graph_solve = Graph(graph.chains, nodes_idx);

    println("graph_energy ", graph_energy(coords, graph_solve))

    G = zeros(Float32, size(coords)...)
    result = Optim.optimize(x->graph_energy(x, graph_solve),
                            (G, x)->graph_gradient!(G, x, graph_solve),
                            coords)

    # convert back to a "normal" graph
    coords = result.minimizer;
    nodes =  Dict(n => coords[:, idx]
                    for (n, idx) in pairs(graph_solve.nodes))

    graph = Graph(graph.chains, nodes);
    return graph
    end;