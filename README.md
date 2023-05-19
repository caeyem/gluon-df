# gluon-df
Gluon transverse momenta distribution generator

This a repo containing Monte-Carlo Simulator for Gluon transverse momenta distributions in the saturation region.
Generators for additional related distributions spanning the theory of spin-glass structures, random breaking etc. are also provided.

## Requirements:
* <code>python</code> (min 3.7m with numpy, scipy, matplotlib) and
* <code>matlab/ octave</code> (tested on 2021b, with statistics and machine learning toolbox)

## Included files:

* <code>mft_gluon_spinglass_common.m</code> : Sets up common global variables and parameters for running mft_gluon_spinglass_\<model\> models
* <code>mft_gluon_spinglass_beta.m</code> : Computes probability distribution functions of weights of maximum gluon transverse momenta cluster valleys (using an analogous spin glass phase from polymer theory and mean field theory -mft), along with second maximum distributions and the distribution of gluon/ sping glass overlap functions. This is achieved using monte-carlo sampling based on samples derived from a beta distribution.
* <code>mft_gluon_spinglass_piecewisedf.m</code> : Similar to the mft_gluon_spinglass_beta.m, but tests the setup of sampling distributions for an alternate approach based on monte-carlo sampling using a piecewise linear approximate distribution.
* <code>mft_gluon_spinglass_piecewise_params.m</code> :  Similar to the mft_gluon_spinglass_beta.m, implementation using an alternate approach based on monte-carlo sampling from a piecewise linear approximate distribution and inverse cdf.
* <code>mft_gluon_spinglass_intermediatedf.m</code> : Supporting code to validate intermediate distributions that taken together yield the overlap distributions from mft_gluon_spinglass_beta.m
* <code>generate_overlap_distribution.py</code> : Implementation of mft_gluon_spinglass_beta.m in python
* <code>gluon_spinglass_theo_overlapdf.m</code> : Generates theoretical overlap distributions for Spin glass structure using alternate model.
* <code>gluon_theo_overlapdf.m</code> : Generates theoretical overlap distributions for Random Map model, Spin glass structures and Breaking Intervals.
* <code>spinglass_finite_moments.m</code> : Symbolic evaluation of spin glass overlap distributions using method of moments (finite)
* <code>gluon_scratchpad.m</code> : Scratch pad for distributions related to gluon saturation

* <code>random_breaking.m</code> : Theoretical and monte-carlo simulated distributions for randomly broken intervals
* <code>random_energy_model.m</code> : Theoretical and monte-carlo simulated distributions for random energy model (spin glass)
* <code>random_map.m</code> : Theoretical and monte-carlo simulated distributions for random map disorder system with deterministic dynamics using inverse transform function
* <code>attractor_random_map.m</code> : Attractor random map model support code
* <code>kauffman_multivalley.m</code> : Support code exploring multivalley structure in the Kauffman model
* <code>magnetism.m</code> : Support code for equivalent models of magnetism

To test, run the params setup file <code>mft_gluon_spinglass_common.m</code>, followed by any of other programs for various gluon, spinglass models.

