# Monte-Carlo Simulation for Multiple Assets
This is a MC simulation program which incorporate general type of stochastic-differential equations (SDE). Drifts and volatility terms are expressed in terms of the power of state variables. For that reason, it is possible to simulate various types of stochastic models if one know how to construct a SDE.

# Usage
- Compilation
  - go to "build" directory and run "cmake .." and "make"
- Run
  - The input files are located in "input" directory. 
  - To run the code using "BSM.inp", in the "build" directory, type "./MultiAsset ../input/BSM.inp"
  - The output is saved in "output" directory. The values of state variables throughout the simulations will be saved. For example, "x0.csv" contains the result of 0-th state variable at each time step. Each row represents a path, and each column shows values at a time step.
- Input file
  - The input directoy is "input"
  - The string after '#' on each line is commented.
- Example (correlated Black-Scholes Model)
  - An example of two correlated Black-Scholes model with a constraint is given
  - Declaration of drift term consists of (index of LHS state variable, coefficient for a drift, index of RHS state variable, the power of RHS)
  - drift 0,-0.105,0,0 # Equivalent to dx[0] += -0.105 * x[0]^0 * dt = -0.105 * dt
  - drift 1,-0.20125,1,0 # Equivalent to dx[1] += -0.20125 * x[1]^0 * dt = -0.20125 * dt
  - Volatility terms in a SDE are defined in the same way. 
  - volatility 0,0.1,0,0 # Equivalent to dx[0] += 0.1 * x[0]^0 * dW[0] = 0.1 * dW[0] 
  - volatility 1,0.05,1,0 # Equivalent to dx[1] += 0.05 * x[1]^0 * dW[1] = 0.05 * dW[1] 
  - x0 0.0 0.01 # This specifies the two initial values for the state variables. x[0] = 0.0; x[1] = 0.01
  - T 1 # The maturity is 1-year
  - Nt 365 # The number of annual time steps
  - num_paths 10000 # The number of paths for a MC simulation
  - correlation 0,1,0.75 # Constant correlation between state variables and means corr(x[0],x[1]) = 0.75
  - constraint min,1,0.0 # x[1] > 0.0; Applying a constraint is useful for non-negative variables like stochastic volatilities

# License
Copyright © 

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

