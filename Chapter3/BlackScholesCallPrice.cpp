//
//  
//  Black-Scholes Call Price
//
/*
    In this program we will implement the Black-Scholes call option pricing formula.
 
    The Black-Scholes model is a mathematical model that is used to find the theoretical price
    of European-style options where this type of option is a financial derivative contract
    that can be exercised only on the expiration date of the option.
    
    A person who holds the option has the right to buy or sell the asset at the specified
    strike price only at the predetermined expiration date.
 
    From [C++ for Financial Mathematics - John Armstrong]
 
    "The Black-Scholes Formulae - if one uses a q measure geometric Brownian
    motion model to price stock options then the price of a European call option
    at time 0 with maturity T and strike K is:
    
    C = N(d1)S0 - N(d2)Kexp(-rT)
 
    where N is the cumulative distribution function of the standard normal
    distribution." [Definition provided in C++ for Financial Mathematics - John Armstrong]
 
    Note that:
 
    d1 = 1/(sigma * sqrt(T)) * (log(S/K) + (r + sigma**2 / 2) * sqrt(T))
 
    d2 = 1/(sigma * sqrt(T)) * (log(S/K) + (r - sigma**2 / 2) * sqrt(T))
 
    Further note that in exercise 3.9.4 we are shown that the cumulative normal
    distribution can be approximated by the following:
 
    If x >= 0, set k = 1 / (1 + 0.2316419*x)
 
    Then a good approximation can be found by:
 
    1 - 1/(sqrt(2*PI)) * exp(-x**2 / 2) * k * (0.319381530 +
    k * (-0.356563682 + k * (1.781477937 +
    k * (-1.821255978 + 1.330274429k))))
 
    If x < 0 then we can use the same formula to find 1 - N(-x)
 
    In this program we will approximate the cumulative normal distribution using
    the above information.
 
    Breaking down the formula we have the following:
    
    C(S, K, r, T, sigma) = S * N(d1) - K * exp(-rT) * N(d2)
 
    Where C is the call option price, S is the current price of the asset, K
    is the strike price of the option, r is the risk free interest rate, T is the
    time to expiration (in years), sigma is the volatility of the underlying
    asset's returns, and N(d1) and N(d2) are the cumulative distribution functions
    of the standard normal distribution where we have d1, d2 defined above.
 */
//

#include <iostream>
#include <string>

class BSCP {
    
public:
    // Member function to approximate the cumulative normal distribution
    double normcdf(double);
    // Member function for computing d1
    double computeD1(double, double, double, double, double);
    // Member function for computing d2
    double computeD2(double, double, double, double, double);
    // Member function for computing the Black-Scholes Call Price
    double callPrice(double, double, double, double, double, double);
    
private:
    // Private member variables
    double x;                    // Value to pass to the normcdf function
    double volatility;           // Volatility - degree of fluctuation in the underlying asset's price
    double time;                 // Time - Time to expiration measured in years
    double stockPrice;           // Current Stock Price - market price of the underlying stock
    double strikePrice;          // Strike Price - price at which the option holder can buy the underlying asset
    double riskFreeRate;         // Risk-free interest rate
};

/**
 *
 * @brief Calculate the cumulative distribution function (CDF) of the standard normal distribution.
 *
 * This function computes the CDF of the standard normal distribution for the given input 'x'.
 * It uses the following approximation formula for efficiency:
 * CDF(x) ≈ 1 - (1 / sqrt(2 * π)) * exp(-0.5 * x * x) * approx,
 * where 'approx' is an approximation polynomial derived from the error function.
 *
 * @param x The value for which to calculate the CDF.
 * @return The CDF value for the given 'x'.
 *
 * This function handles both positive and negative 'x' values appropriately.
 *
 * Example usage:
 * @code
 * double cdf = BSCP::normcdf(1.5);
 * @endcode
 */
double BSCP::normcdf(double x){
    const double g = 0.2316419;
    double k = (1.0 / (1.0 + g * x));
    double approx = k*(0.319381530 + k*(-0.356563782 + k*(1.781477937 + k*(-1.821255978 + 1.330274429*k))));
    
    if (x >= 0.0) {
        return (1.0 - (1.0/(pow(2*M_PI,0.5)))*exp(-0.5*x*x) * approx);
    } else {
        return 1.0 - normcdf(-x);
    }
}

/**
 * @brief Calculate the 'd1' parameter for the Black-Scholes option pricing model.
 *
 * This function computes the 'd1' parameter used in the Black-Scholes option pricing model.
 * The 'd1' parameter is essential for calculating the option's price.
 *
 * @param sigma The volatility of the underlying asset's returns.
 * @param T The time to expiration (in years) of the option.
 * @param S The current price of the underlying asset.
 * @param K The strike price of the option.
 * @param r The risk-free interest rate.
 * @return The calculated 'd1' parameter.
 *
 * Example usage:
 * @code
 * double d1 = BSCP::computeD1(0.2, 1.0, 100.0, 95.0, 0.05);
 * @endcode
 */
double BSCP::computeD1(double sigma, double T, double S, double K, double r) {
    double d1;
    d1 = (log(S / K) + (r + (sigma*sigma / 2)) * sqrt(T)) / (sigma * sqrt(T));
    return d1;
}

/**
 * @brief Calculate the 'd2' parameter for the Black-Scholes option pricing model.
 *
 * This function computes the 'd2' parameter used in the Black-Scholes option pricing model.
 * The 'd2' parameter is essential for calculating the option's price.
 *
 * @param sigma The volatility of the underlying asset's returns.
 * @param T The time to expiration (in years) of the option.
 * @param S The current price of the underlying asset.
 * @param K The strike price of the option.
 * @param r The risk-free interest rate.
 * @return The calculated 'd2' parameter.
 *
 * Example usage:
 * @code
 * double d2 = BSCP::computeD2(0.2, 1.0, 100.0, 95.0, 0.05);
 * @endcode
 */
double BSCP::computeD2(double sigma, double T, double S, double K, double r) {
    double d2;
    d2 = (log(S / K) + (r - (sigma*sigma / 2)) * sqrt(T)) / (sigma * sqrt(T));
    return d2;
}

/**
 * @brief Calculate the price of a European call option using the Black-Scholes model.
 *
 * This function computes the price of a European call option using the Black-Scholes option pricing model.
 *
 * @param d1 The 'd1' parameter calculated for the option.
 * @param S The current price of the underlying asset.
 * @param d2 The 'd2' parameter calculated for the option.
 * @param K The strike price of the option.
 * @param r The risk-free interest rate.
 * @param T The time to expiration (in years) of the option.
 * @return The calculated price of the call option.
 *
 * Example usage:
 * @code
 * double d1 = BSCP::computeD1(0.2, 1.0, 100.0, 95.0, 0.05);
 * double d2 = BSCP::computeD2(0.2, 1.0, 100.0, 95.0, 0.05);
 * double callPrice = BSCP::callPrice(d1, 100.0, d2, 95.0, 0.05, 1.0);
 * @endcode
 */
double BSCP::callPrice(double d1, double S, double d2, double K, double r, double T){
    double C;
    C = (normcdf(d1) * S) - (normcdf(d2) * K * std::exp(-r * T));
    return C;
}

int main(int argc, const char * argv[]) {
    // Create a Black-Scholes call price object
    BSCP option;
    
    // Set values for Sigma (Volatility), T (Time to Expiration), S (Current Stock Price)
    // K (Strike Price), and r (Risk-free interest rate)
    double sigma = 0.15;
    double T = 1.0;
    double S = 300.00;
    double K = 200.00;
    double r = 0.03;
    
    // Compute d1 and d2
    double d1 = option.computeD1(sigma, T, S, K, r);
    double d2 = option.computeD2(sigma, T, S, K, r);
    
    // Compute the call price
    double callPrice = option.callPrice(d1, S, d2, K, r, T);
    
    // Display volatility, expiration, stock price, strike price, interest, and the computer call price
    std::cout << "Volatility: " << sigma << std::endl;
    std::cout << "Time to Expiration (Years): " << T << std::endl;
    std::cout << "Current Stock Price: " << S << std::endl;
    std::cout << "Strike Price: " << K << std::endl;
    std::cout << "Risk-Free Interest Rate: " << r << std::endl;
    std::cout << "Black-Scholes Call Price: " << callPrice << std::endl;
    
    return 0;
}
