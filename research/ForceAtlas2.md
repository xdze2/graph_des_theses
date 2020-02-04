# ForceAtlas2 Algorithm

ForceAtlas2, a Continuous Graph Layout Algorithm for Handy Network Visualization Designed for the Gephi Software, M. Jacomy, T. Venturini, S. Heymann, and M. Bastian, PLoS ONE, vol. 9, no. 6, p. e98679, Jun. 2014, doi: 10.1371/journal.pone.0098679.


## Performance Optimization, the issue of speed:
- speed = time step or step length
- oscillation = temperature
- tradeoff speed vs. precision

Linked ideas:
- local temperature GEM [ref 5]
- adaptive cooling, Yifan Hu [ref 7]
- simulated annealing [ref 6]

ForceAtlas2: local temperatures + adaptive cooling  + continuous layout (i.e. live interactive rendering)

> Our strategy is to determine the **speed
of each node by observing its oscillation**, like in GEM [5]. But our
implementation is actually quite different.

"in ForceAtlas2 the speed is different for every node" i.e. local temperature



- The swing of the node `i`:

        swg(i) = | Fi(t) - Fi(t-dt) |

--- curvature of the potential well

- node speed `s`:

                        k sG
        s(i) =  ---------------------
                1 + sG sqrt[ swg(i) ]

i.e. swing ++ --> speed --


- the global speed (adaptative cooling):

                   traction(G)
        sG = tau ---------------
                      swg(G)

`f(G)` means weighted by degrees sum over all nodes,   
and traction is `|Fi(t) + Fi(t+dt)|`...
