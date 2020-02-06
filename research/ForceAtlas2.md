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


- Displacement of a node `ì`:

        D(i) = s(i)*F(i)
        i.e.  xi = xi + si*Fi


- The swing of the node `i`:

        swg(i) = | Fi(t) - Fi(t-dt) |

--- curvature of the potential well

- node speed `s`:

                        k sG
        s(i) =  ---------------------
                1 + sG sqrt[ swg(i) ]

`sG` is the global speed.

(i.e. swing ++ --> speed --)

With the constant parameter `k = 0.1`.
In addition, there is a limit on the local speed: `s(i) < k_max /|F_i|` with `k_max=10`.



- the global speed (adaptative cooling):

                   traction(G)
        sG = tau ---------------
                      swg(G)

`f(G)` means weighted by (degrees + 1) sum over all nodes,   
and traction is `|Fi(t) + Fi(t+dt)|/2` which is the force average...

> the ratio `tau` represents the tolerance to swinging and is set by the
user.

> This amount is set by the user as ‘‘Tolerance (speed)’’. In Gephi’s
implementation, we set three default values: 0.1 under 5000 nodes,
1 up to 50000 nodes and 10 above 50000. We now describe how
this feature works.

> NB: During our tests we observed that an excessive rise of the
global speed could have a negative impact. That is why we limited
the increase of global speed s (t) (G) to 50% of the previous step
s (t{1) (G).

## Notes
- named "damped dynamics method" from the [LAMMPS min_style doc ](https://lammps.sandia.gov/doc/min_style.html)
- see Adagrad, RMSProp, Adadelta and Adam methods (book Algorithms for Optimization, M. J. Kochenderfer & T. A. Wheeler)

- What about simulated annealing?  pure monte-carlo, i.e. without solving the dynamics, local random sampling is the phase space
