# FSAR
This is the R simulation code for the FSAR method.

FSAR_functions.R contains all the necessary functions for data generation and estimation.

Main.R outputs the simulation results. To run the stochastic block model, set the argument "net.type" in "get.func.network" to be "SBM"; To run the power-law distributed model, set the same argument to be "power". 

Plot.R outputs the plot for the estimation of $\rho(t)$ and $\beta(t)$.

Result_SBM100-50.rda is the simulation result for the stochastic block model when the sample size is 100 and maximum number of points observed for a single trajectory is 50.

