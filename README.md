# 2d-euler-gas-simulation
Computational model of breaking wave for NYU computation physics class in 2009.

## 2nd Order Euler Equations in 2D with Gravity

### Instructions

1. make
1. python run.py

You can change the initial conditions in 
param.cfg, just change the values for velocity, density and pressure
for Mr and Ml. In each initial condition, specified by InitialType there 
are two gases, one corresponding to Mr, and the other to Ml.

In run.py you can specify the resolution, the BC type, the Initial Condition
type and , the time step, and the number of runs. Also you can specify where
you want to save the figures it creates to.

In runAgain.py, you can run the data from where it was just left off, you need to specify where to save the figures to, the time step, and the number of graphs wanted.

It is important to have a file in the results folder that corresponds to the destination specified in run.py and runAgain.py
