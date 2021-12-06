# Heston Stochastic-Local Volatility

Heston SLV is a model which combines the parameteric local and stochastic local volatilities. The conventional Heston model is know to not be able to fit a volatility surface with short maturities. In order to circumvent this problem, the combination of local and stochastic volatilities has been created. The SLV model of the Heston type is chosen since the well-know Heston model can be easily extended by applying the leverage function in terms of an underlying asset proce. As the dynamics of a contract under the risk-neutral measure could be derived, the price of a contract can be interpreted as the solution to this PDE. Because there is a cross-term between the stock price and volatility processes, the altering-direction implicit (ADI) is deployed. For the sake of simplicity, the method of douglas is used.

# Python codes
- BS.py: formulas to compute the BS model and the implied volatility.
- fd_adi.py: It contains the functions for the finite-difference scheme and the class of tridiagonal matrix. Tridiagonal matrix can be useful since it requires much less memory than dense matrices when one should deal with high dimensional problems.

# Reference
- Mathematical Modeling and Computation in Finance: With Exercises and Python and MATLAB Computer Codes,
by C.W. Oosterlee and L.A. Grzelak, World Scientific Publishing Europe Ltd, 2019.
- https://en.wikipedia.org/wiki/Alternating-direction_implicit_method
