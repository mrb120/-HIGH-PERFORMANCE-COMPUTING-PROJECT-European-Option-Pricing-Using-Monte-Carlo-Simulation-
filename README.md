# -HIGH-PERFORMANCE-COMPUTING-PROJECT-European-Option-Pricing-Using-Monte-Carlo-Simulation-
This report consist of comparison of using shared memory, distributed memory and hybrid approach 
applied to the problem of European call option using the Monte Carlo simulation method. Monte Carlo 
simulation is used for modelling uncertainty and predicting future outcomes. A European option is a 
financial contract that gives the holder the right to buy a stock at a fixed price, but only on a specific 
future date. To simulate the option’s value, the Black Scholes model is used, which accounts for factors 
like current stock price, volatility, interest rate, and time to maturity.  
 
However, Monte Carlo simulation involves generating millions of possible future stock prices using 
random numbers and calculating the average potential profit (payoff) across all simulations.  
This report analyses the accuracy, and execution time comparison for each parallelization method with 
the serial code.  
For the comparison following parallelization techniques were used. 
  Shared Memory – OpenMP 
  Distribute Memory -MPI 
  Hybrid Approach - MPI + OpenMP 
