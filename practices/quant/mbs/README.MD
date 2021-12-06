# Mortgage-Backed Securities
The model for pricing a MBS with prepayment risk is implemented. The intensity process for the prepayment risk takes the form of the CIR-process. Because the price of a MBS is expressed as the function of the expectations of the discounted intensity process, this complicates the solution. To solve this problem, the finite-difference scheme with the Crank-Nicolson method is chosen as this scheme is efficient for 1-D problems.

# Julia Codes
- fd.jl: The implementation of the Crank-Nicolson scheme
- fd_functions.jl: Functions for constructing the finite-difference schemes

# Reference
- INTENSITY-BASED VALUATION OF RESIDENTIAL MORTGAGES: AN ANALYTICALLY TRACTABLE MODEL, VYACHESLAV GOROVOY, Mathematical Finance, Vol. 17, No. 4 (October 2007), 541–573
