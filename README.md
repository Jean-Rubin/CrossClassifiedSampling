# Cross Classified Sampling variance estimation

This repository is associated to the work done in the following paper : [Rubin, J. and Chauvet, G. (2024) - Asymptotic properties of cross-classified sampling designs]().

Suppose that we would like to sample a finite population $\mathcal{U} = \prod_{d=1}^D \mathcal{U}_d$ which can be seen as the cartesian product of $D$ finite populations.
One approach would be to use a cross-classified sampling (CCS) design, which produce $D$ independent samples $S_d$ in the population $\mathcal{U}_d$, $d = 1, \dots, D$ and use as a sample of $\mathcal{U}$ the cartesian product of them: 
$$S = \prod_{d = 1}^D S_d$$

This repository is a Julia project allowing to assess the performance of multiple variance estimations when producing indicators of a simulated population sampled with CCS.
We can distinguish two families of methods :
- Linearization methods. There are three different methods corresponding to different levels of simplifications of the general variance formula of an Horvitz-Thompson estimator.
- Weighted bootstrap methods. There are three methods, one from Skinner (2015), one based on Gross (1980) and one based on Rao, J. Wu, C., and Yue, K. (1992).


## Run the code

Clone the project, then execute in the Julia REPL
```shell
julia> ]activate
```
to activate the environment. You can then execute
```shell
julia> include("main.jl")
```
to run the main script that will produce the simulations.
