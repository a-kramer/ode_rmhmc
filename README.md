RMHMC for steady state ODE models
=================================

This collection of matlab and octave scripts was used to perform
simulations for the publication «Hamiltonian Monte Carlo Methods for
Efficient Parameter Estimation in Steady State Dynamical Systems».
You should read the publication first. A link will be provided as soon
as the manuscriupt becomes available.

We will improve the documentation and hope that this project will
mature to an implementation that can be used more easily by everyone.

How to replicate the results
============================

The reference implementations of RMHMC, HMC and SMMALA require the
matlab toolbox [SBPOP](http://www.sbtoolbox2.org/main.php) to be
installed.

The steady state adapted algorithms do not integrate the model and do
not require this package.

The model files for the steady state adapted algorithms (NR prefix)
have been built using octave forge's [symbolic package](http://octave.sourceforge.net/symbolic/). 
However, the symbolic calculations required can be done using any other
software capable of symbolic calculations.

The models are named _MAPK_, _Mma_ and _Mifa_. You will find setup files for
each of them (called ```MAPKsetup.m```, ```MmaSetup.m``` and
```MifaSetup.m```). These files set variables that are not algorithm
specific. But some of the values may impact an algorithm's performance
(such as the step size).

If your SBPOP package is installed correctly you should be able to run
the procedures called ```Run_«model»_using_«algorithm».m```.

The only files you need to investigate are the setup files. They
contain prior setups, step sizes and sample sizes.

Inspect the simulation results
===============================

You have acces to the values of all defined options as used by us
during sampling by loading the ```.mat``` files in the ```Results```
folders and checking the options variable.

How to set up your own example
==============================

This is not very convenient yet, since a lot of the model setup is not
automated. If you wish to try it out anyway, then inspect the
```symmodel``` octave files, e.g. ```Models/Insulin/Mma_symmodel.m```
and adjust them to your needs. If you wish you can also inspect the
model files themselves, e.g. ```make_Mma_model.m``` and produce a
similar file for your model using Maple, Maxima, Mathematica or any
other symbolic software package.

If you really want to compare performance, you will have to build an
SBToolBox2 model file as well. But, keep in mind, that the proposed
method (using Newton Raphson and steady state sensitivity analysis)
should be strictly superior to the reference implementation in all
cases (the drawback beeing the restriction to steady states). So,
there is no need to compare RMHMC to the steady state adapted version
described in the publication. You will probably want to compare this
method to the one you have been using in the past.

Files
=====

All files not mentioned here have been used to process the results and
need not necessarily be investigated but are provided nevertheless for
the curious.