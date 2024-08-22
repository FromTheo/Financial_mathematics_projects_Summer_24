*These financial mathematics projects were completed during the summer of 2024. They reflect an in-depth exploration of the methods and tools used in quantitative finance, as well as practical analyses for a better understanding of financial models.*

## Table of contents : 

- [Project 1 : Integral Approximation and Bond Pricing](#project-1--integral-approximation-and-bond-pricing)
- [Project 2 : Vanilla Options for Black-Scholes Model](#project-2--vanilla-options-for-black-scholes-model)
- [Project 3 : Nonlinear Solvers for Bond Yields and Implied Volatility](#project-3--nonlinear-solvers-for-bond-yields-and-implied-volatility)
- [Project 4 : Greeks and Hedging](#project-4--greeks-and-hedging)


## Project 1 : Integral Approximation and Bond Pricing

1. **Methods for Approximating an Integral**
   - 1.1 Midpoint rule 
   - 1.2 Trapezoidal rule 
   - 1.3 Simpson's rule 
   - 1.4 Approximation with tolerance threshold

2. **Bond pricing**
   - 2.1 The Interest Rate Market
     - 2.1.1 The Continuous Interest Rate
     - 2.1.2 The Compound Interest Rate
     - 2.1.3 Relationships Between Rates
     - 2.1.4 Yield of a bond, par yield
   - 2.2 Duration
   - 2.3 Convexity 
   - 2.4 Approximation with tolerance threshold

## Project 2 : Vanilla Options for Black-Scholes Model

1. **Approximation of the Cumulative Distribution Function (cdf) of a Standard Gaussian**
   - 1.1 Zelen & Severo's formula (1964)
   - 1.2 Byrc's formula (2002B)
   - 1.3 Page's formula (1977)
   - 1.4 Bagby's formula (1995)
   - 1.5 Error Calculation and Execution Time

2. **Simulation of Random Variables with Distribution N(0,1)**
   - 2.1 Box Muller method
   - 2.2 Marsaglia method 
   - 2.3 Time measurement

3. **Black-Scholes Model**

4. **Black-Scholes formulas: Europeans Put & Call**

5. **Options strategies**
   - 5.1 Spreads
     - 5.1.1 Bull spread
     - 5.1.2 Bear spread
     - 5.1.3 Butterfly spread
     - 5.1.4 Box spread
     - 5.1.5 Calendar spread
     - 5.1.6 Diagonal spread
   - 5.2 Combinations
     - 5.2.1 Straddles
     - 5.2.2 Strips and straps
     - 5.2.3 Strangles
  
## Project 3 : Nonlinear Solvers for Bond Yields and Implied Volatility

1. **Numerical Methods for Finding a Zero in 1D**
   - 1.1 Bisection method 
   - 1.2 Newton's method   
   - 1.3 Secant method

2. **Numerical Methods for Finding a Zero in N-Dimension"**
   - 2.1 N-dimensional Newton's method 
   - 2.2 Approximate Newton's method 

3. **Computing bond yields**

4. **Computing implied volatility**
   - 4.1 Case of the call option
   - 4.2 Case of the put option 

5. **Bootstrapping for finding zero rate curves**

## Project 4 : Greeks and Hedging 

1. **Delta (Δ)**
   - 1.1 Calculation for a call / put 
   - 1.2 Implementation 
   - 1.3 Finite Difference Approximations
   - 1.4 Δ-hedging

2. **Gamma (Γ)** 
   - 2.1 Calculation for a call / put 
   - 2.2 Implementation
   - 2.3 Maximality for ATM Options
   - 2.4 Finite Difference Approximations
   - 2.5 Γ-hedging

3. **Thêta (θ)**
   - 3.1 Calculation for a call / put 
   - 3.2 Implementation
   - 3.3 Finite Difference Approximations

4. **Rhô (ρ)**
   - 4.1 Calculation for a call / put 
   - 4.2 Implementation
   - 4.3 Finite Difference Approximations

5. **Vega** 
   - 5.1 Calculation for a call / put 
   - 5.2 Implementation and Increase with Maturity 
   - 5.3 Finite Difference Approximations
   - 5.4 Vanna 
   - 5.5 Volga 

6. **Black-Scholes PDE** 
   - 6.1 Proof
   - 6.2 Financial Interpretation and Link with the Greeks

**Appendix : Proof of the formula(☆)**
