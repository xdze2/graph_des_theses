#using(BenchmarkTools)
import (Plots)

# https://github.com/mauro3/Parameters.jl
# import Pkg; Pkg.add("Parameters")
using Parameters


"""
Graph the graph
"""
function graphplot(graph, chains)
    fig = Plots.plot(legend = false, aspect_ratio = :equal)
    #scatter!(coords[:, 1], coords[:, 2], linewidth=2)
    for chain in chains
        chain_xy = map(x -> graph[x].position, chain)
        filter!(xy -> ~isnothing(xy), chain_xy)
        if length(chain_xy) > 0
            x = getindex.(chain_xy, 1)
            y = getindex.(chain_xy, 2)
            Plots.plot!(x, y, linewidth = 3)
        end
    end
    fig
end;


"""
Node struct
"""
mutable struct Node
    neighbors
    position::Union{Nothing,Tuple{Float64,Float64}}
    velocity::Union{Nothing,Tuple{Float64,Float64}}
    force::Union{Nothing,Tuple{Float64,Float64}}
    energy::Union{Nothing,Float64}
    previous_force::Union{Nothing,Tuple{Float64,Float64}}
end;

Node() = Node([], nothing, nothing, nothing, nothing, nothing)


"""
Init graph data structure
"""
function create_graph(chains)
    graph = Dict()
    for chain in chains
        node = get!(graph, chain[1], Node())
        push!(node.neighbors, (nothing, chain[2]))

        for k = 2:length(chain)-1
            node = get!(graph, chain[k], Node())
            push!(node.neighbors, (chain[k-1], chain[k+1]))
        end

        node = get!(graph, chain[end], Node())
        push!(node.neighbors, (chain[end-1], nothing))
    end

    return graph
end;

"""
Init nodes coordinates on a circle
"""
function place_on_circle!(graph, node_ids = [])
    if length(node_ids) == 0
        node_ids = keys(graph)
    end
    N = length(node_ids)
    radius = N / (2 * pi)
    delta_theta = 2 * pi / N
    for (k, node_id) in enumerate(node_ids)
        node = graph[node_id]
        theta = k * delta_theta
        node.position = (radius * cos(theta), radius * sin(theta))
        node.velocity = (0.0, 0.0)
    end

end;

"""
    Velocity Verlet
http://students.iitk.ac.in/projects/wiki/lib/exe/fetch.php?media=2014as:verlet.pdf
"""
function velocity_verlet!(graph, dt)
    dt2_2 = dt^2 / 2

    nodes = filter(n -> ~isnothing(n.position), collect(values(graph)))
    #nodes = values(graph)
    # Positions
    for node in nodes
        dX = (dt2_2) .* node.force
        if norm(dX) > 20.0
            dX = dX ./ norm(dX)
            end;
        node.position = node.position .+ dX
        #node.previous_force = node.force
    end

    # Forces
    # reset_forces!(graph)
    # all_attractive_forces!(graph)
    # all_repulsive_forces!(graph)
    energy = compute_forces!(graph) # and reset

    # Velocities
    for node in nodes
        node.velocity = node.velocity .+
                        (node.force .+ node.previous_force) .* (0.5dt)
    end

    return energy
end;

dot(a, b) = sum(a .* b);
norm(a) = sqrt(dot(a, a));
unit_vector(a) = a ./ norm(a);

@with_kw mutable struct Fire
    N_min = 5
    f_inc = 1.1
    f_dec = 0.5
    alpha_start = 0.1
    f_alpha = 0.99
    dt_max = 0.5
    # variable
    dt = dt_max / 5.0
    N = 0
    alpha = alpha_start
end

#Fire()

# Fire minimization
# Bitzek, Erik, et al. "Structural relaxation made simple."
# Physical review letters 97.17 (2006): 170201.

function step!(graph, fire)

    # loop
    nrj = velocity_verlet!(graph, fire.dt)

    nodes = filter(n -> ~isnothing(n.position), collect(values(graph)))
    #nodes = values(graph)
    P = sum(dot(node.force, node.velocity) for node in nodes)

    for node in nodes
        hat_F = unit_vector(node.force)
        v = node.velocity
        node.velocity = (1 - fire.alpha).*v .+ (fire.alpha*norm(v)).*hat_F
    end

    if P > 0 && fire.N > fire.N_min
        fire.dt = min(fire.f_inc*fire.dt, fire.dt_max)
        fire.alpha = fire.alpha * fire.f_alpha
        fire.N = 0
    elseif P > 0
        fire.N += 1
    end

    if P <= 0
        fire.dt = fire.f_dec * fire.dt
        fire.alpha = fire.alpha_start
        foreach(node -> node.velocity = (0.0, 0.0), nodes)
        fire.N = 0
    end
    return fire.dt, nrj
end;



"""
    reset_forces!(graph)
"""
function reset_forces!(graph)
    for node in values(graph)
        if ~isnothing(node.position)
            node.previous_force = node.force
            node.force = (0.0, 0.0)
        end
    end
end;


const coeff = 5.0
function three_body_force(u, v)
    #  points a, b, c
    #  u = b - a
    #  v = c - b
    #  F->a = -dE/da = dE/du
    #  F->c = -dE/dc = -dE/dv
    #  F->b = -F->a - F->c = -dE/du + dE/dv
    #
    #  E = -cos( (a, b, c) )

    x1, y1 = u
    x2, y2 = v
    norm2_u = x1^2 + y1^2
    norm2_v = x2^2 + y2^2
    cross = x1 * y2 - x2 * y1

    cross_nunv = coeff * cross / sqrt(norm2_u * norm2_v)
    s3 = cross_nunv / norm2_u
    s4 = cross_nunv / norm2_v

    dEdu = (s3 * y1, -s3 * x1)
    dEdv = (-s4 * y2, s4 * x2)
    return dEdu, .-dEdu .+ dEdv, .-dEdv
end;

# ========  Forces  ========
# Spring (+2 power law):
#   E = d2 = x2 + y2
#   F = (-dE/dx, -dE/dy) = (-2*x, -2*y)
#   u = b - a
#   F->b
spring_force(u) = -2 .* u;
spring_energy(u) = dot(u, u);

# Repulsion (-2 power law):
#   E = 1/d2 = 1/(x2 + y2)
#   F = (-dE/dx, -dE/dy) = (-2x/d4, -2y/d4)
const EPSILON = 1e-4
function repulsive_force(u)
    d4 = (sum(u .^ 2) + EPSILON)^2
    return 2 / d4 .* u
end;
function repulsive_energy(u)
    d2 = (sum(u .^ 2) + EPSILON)
    return 1 / d2
end;


# test
# repulsive_force_Coulomb((1., 1.))

function compute_forces!(graph)
    energy = 0.0
    # Compute all forces
    reset_forces!(graph)
    three_body = false
    # Link forces
    for node in values(graph)
        if isnothing(node.position)
            continue
        end

        for (left_name, right_name) in node.neighbors
            if !isnothing(left_name) && !isnothing(graph[left_name].position)
                left = graph[left_name]
                # Compute pair force
                u = node.position .- left.position
                F = spring_force(u)
                energy += spring_energy(u)
                left.force = left.force .- F
                node.force = node.force .+ F

                if three_body &&
                   !isnothing(right_name) &&
                   !isnothing(graph[right_name].position)
                    right = graph[right_name]
                    v = right.position .- node.position
                    # Compute three body force
                    # energy += ...
                    Fa, Fb, Fc = three_body_force(u, v)
                    left.force = left.force .+ Fa
                    node.force = node.force .+ Fb
                    right.force = right.force .+ Fc
                end

            end

        end
    end

    # Repulsive forces
    for (k, node_a) in enumerate(values(graph))
        if isnothing(node_a.position)
            continue
        end
        for node_b in Iterators.take(values(graph), k - 1)
            if ~isnothing(node_b.position)
                u = node_b.position .- node_a.position
                energy += repulsive_energy(u)
                Fb = repulsive_force(u)
                node_a.force = node_a.force .- Fb
                node_b.force = node_b.force .+ Fb
            end
        end
    end

    return energy
end;

# function all_repulsive_forces!(graph)
#     # Repulsive forces
#     for (k, node_a) in enumerate(values(graph))
#         if isnothing(node_a.position)
#             continue
#             end;
#         for node_b in Iterators.take(values(graph), k-1)
#              if ~isnothing(node_b.position)
#                 u = node_b.position .- node_a.position
#                 Fb = repulsive_force_Coulomb(u)
#                 node_a.force = node_a.force .- Fb
#                 node_b.force = node_b.force .+ Fb
#                 end;
#             end;
#         end;
#     end;
