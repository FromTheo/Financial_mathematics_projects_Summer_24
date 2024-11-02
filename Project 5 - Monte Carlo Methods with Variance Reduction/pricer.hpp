#include <ctime> 
#include <cmath>
#include <random> 
#include <vector>
#include <fstream>
#include <utility> 
#include <iostream>
#include <Eigen/Core> 
#include <Eigen/Cholesky> 

#include "variance_calculator.hpp"

template<int d> 
class PriceIndex; 

template<int d> 
std::ostream& operator<<(std::ostream& out, const PriceIndex<d>& P); 

template<int d> 
Eigen::Matrix<double, d, d> generateCorrelationMatrix(double rho) 
{
    Eigen::Matrix<double, d, d> correlationMatrix = Eigen::Matrix<double, d, d>::Identity();
    for (int i = 0; i < d; i++) 
    {
        for (int j = i + 1; j < d; j++) 
        {
                correlationMatrix(i, j) = rho;
                correlationMatrix(j, i) = rho;
        }
    }
    return correlationMatrix;
}

double norm_cdf(double x) 
{
    return 0.5 * (1.0 + std::erf(x / std::sqrt(2.0)));
}

template<int d> 
class PriceIndex
{
    private :
    double T;                                               // Maturité
    double r;                                               // Taux d'intérêt
    Eigen::Matrix<double, d, 1>  s0;                        // Prix initiaux 
    Eigen::Matrix<double, d, 1>  sigma;                     // Volatilités
    Eigen::Matrix<double, d, d> correlationMatrix;          // Matrice de corrélation 
    Eigen::Matrix<double, d, 1>  weights;                   // Poids de chaque actif sous-jacent dans l'indice 

    public :    
    PriceIndex(double T, double r, Eigen::Matrix<double, d, 1> s0, Eigen::Matrix<double, d, 1>  sigma, Eigen::Matrix<double, d, d> correlationMatrix, Eigen::Matrix<double, d, 1>  weights) : T(T), r(r), s0(s0), sigma(sigma), correlationMatrix(correlationMatrix), weights(weights) {}

    Eigen::Matrix<double, d, 1> simulateAssetPrices(Eigen::Matrix<double, d, 1> Z) const; 

    friend std::ostream& operator<< <>(std::ostream& out, const PriceIndex<d>& P); 

    template<int d_, typename RNG>
    friend std::pair<double, double> MonteCarlo_CallBasket(PriceIndex<d_>& P, double K, RNG& G, int N);

    template<int d_, typename RNG>
    friend std::pair<double, double> MonteCarlo_CallBasket_Antithetic(PriceIndex<d_>& P, double K, RNG& G, int N);

    template<int d_, typename RNG>
    friend std::pair<double, double> MonteCarlo_CallBasket_ControlVar(PriceIndex<d_>& P, double K, RNG& G, int N);

    template<int d_, typename RNG>
    friend std::pair<double, double> MonteCarlo_CallBasket_ControlVar_Antithetic(PriceIndex<d_>& P, double K, RNG& G, int N);

    template<int d_>
    friend double closed_formula(PriceIndex<d_>& P, double K);
};


template<int d> 
Eigen::Matrix<double, d, 1> PriceIndex<d>::simulateAssetPrices(Eigen::Matrix<double, d, 1> Z) const
{   
    Eigen::Matrix<double, d, 1> S_T;
    Eigen::Matrix<double, d, d> L = correlationMatrix.llt().matrixL();   // Décomposition de Cholesky 

    Eigen::Matrix<double, d, 1> W_T = L * Z * sqrt(T);
    for (int i = 0; i < d; i++) 
    {
        S_T(i) = s0(i) * std::exp((r - 0.5 * sigma(i) * sigma(i)) * T + sigma(i) * W_T(i));
    }
    return S_T; 
}

template<int d> 
std::ostream& operator<<(std::ostream& out, const PriceIndex<d>& P)
{   
    out << P.T << std::endl; 
    out << P.r << std::endl; 
    out << P.s0 << std::endl; 
    out << P.sigma << std::endl; 
    out << P.correlationMatrix << std::endl; 
    out << P.weights << std::endl; 
    return out; 
}

template<int d, typename RNG>
std::pair<double, double> MonteCarlo_CallBasket(PriceIndex<d>& P, double K, RNG& G, int N)
{
    double discount_factor = std::exp(-P.r*P.T); 
    VarianceCalculator V;
    double sum = 0.0; 
    std::normal_distribution<double> norm(0.0,1.0);
    for(int i = 0; i < N; i++)
    {
        Eigen::Matrix<double, d, 1> Z;
        for (int j = 0; j < d; j++) Z(j) = norm(G);
        Eigen::Matrix<double, d, 1> S_T = P.simulateAssetPrices(Z); 
        double I_T = P.weights.dot(S_T); 
        double x = std::max(I_T - K, 0.0); 
        sum += x;
        V.update(x); 
    }
    return std::make_pair((sum / N)*discount_factor, V.variance());
}

template<int d, typename RNG>
std::pair<double, double> MonteCarlo_CallBasket_Antithetic(PriceIndex<d>& P, double K, RNG& G, int N)
{
    double discount_factor = std::exp(-P.r*P.T); 
    VarianceCalculator V; 
    double sum = 0.0;
    std::normal_distribution<double> norm(0.0, 1.0);
    for (int i = 0; i < N; ++i)
    {
        Eigen::Matrix<double, d, 1> Z;
        for (int j = 0; j < d; j++) Z(j) = norm(G);

        Eigen::Matrix<double, d, 1> S_T_pos = P.simulateAssetPrices(Z);
        Eigen::Matrix<double, d, 1> S_T_neg = P.simulateAssetPrices(-Z);

        double I_T_pos = P.weights.dot(S_T_pos);
        double I_T_neg = P.weights.dot(S_T_neg);

        double payoff_pos = std::max(I_T_pos - K, 0.0);
        double payoff_neg = std::max(I_T_neg - K, 0.0);
        sum += (payoff_pos + payoff_neg) / 2.0;
        V.update((payoff_pos + payoff_neg) / 2.0); 
    }

    return std::make_pair((sum / N) * discount_factor, V.variance()); 
}

template<int d>
double closed_formula(PriceIndex<d>& P, double K)
{
    double discount_factor = std::exp(-P.r*P.T); 
    double mu = 0.0; 
    for(int i = 0; i < d; ++i)
    {
        mu += std::log(P.s0(i)) + (P.r-P.sigma(i)*P.sigma(i)/2)*P.T;
    }
    mu /= d; 

    double gamma2 = 0.0; 
    for(int i = 0; i <d ; ++i)
    {
        for(int j = 0; j < d; ++j)
        {
            gamma2 += P.sigma(i)*P.sigma(j)*P.correlationMatrix(i,j); 
        }
    }
    gamma2 /= d*d;
    gamma2 *= P.T;  

    double D = (gamma2 + mu - std::log(K))/std::sqrt(gamma2); 
    return (std::exp(gamma2/2+mu)*norm_cdf(D) - K*norm_cdf(D-std::sqrt(gamma2))) * discount_factor;  
}

template<int d, typename RNG>
std::pair<double, double>  MonteCarlo_CallBasket_ControlVar(PriceIndex<d>& P, double K, RNG& G, int N)
{
    VarianceCalculator V; 
    double sum = 0; 
    double discount_factor = std::exp(-P.r*P.T); 
    
    std::normal_distribution<double> norm(0.0, 1.0);

    for (int i = 0; i < N; ++i)
    {
        Eigen::Matrix<double, d, 1> Z;
        for (int j = 0; j < d; j++) Z(j) = norm(G);

        Eigen::Matrix<double, d, 1> S_T = P.simulateAssetPrices(Z);
        double I_T = P.weights.dot(S_T);
        double payoff_X = std::max(I_T - K, 0.0);
        sum += payoff_X; 

        double prod = 1.0; 
        for(int j = 0; j < d; ++j)
        {
            prod *= S_T(j); 
        }
        
        double payoff_Y = std::max(std::pow(prod, 1.0/d) - K, 0.0);
        sum -= payoff_Y; 
        V.update(payoff_X - payoff_Y); 
    }

    return std::make_pair((sum / N) * discount_factor + closed_formula(P,K), V.variance());
}

template<int d, typename RNG>
std::pair<double, double> MonteCarlo_CallBasket_ControlVar_Antithetic(PriceIndex<d>& P, double K, RNG& G, int N)
{
    VarianceCalculator V; 
    double discount_factor = std::exp(-P.r*P.T); 
    double sum = 0.0; 
    double closed_form_val = closed_formula(P, K); 

    std::normal_distribution<double> norm(0.0, 1.0);

    for (int i = 0; i < N; ++i)
    {
        Eigen::Matrix<double, d, 1> Z;
        for (int j = 0; j < d; j++) Z(j) = norm(G);

        Eigen::Matrix<double, d, 1> S_T_pos = P.simulateAssetPrices(Z);
        double I_T_pos = P.weights.dot(S_T_pos);
        double payoff_X_pos = std::max(I_T_pos - K, 0.0);

        Eigen::Matrix<double, d, 1> S_T_neg = P.simulateAssetPrices(-Z);
        double I_T_neg = P.weights.dot(S_T_neg);
        double payoff_X_neg = std::max(I_T_neg - K, 0.0);

        double payoff_X = (payoff_X_pos + payoff_X_neg) / 2.0;

        double prod_pos = 1.0; 
        for(int j = 0; j < d; ++j)
        {
            prod_pos *= S_T_pos(j); 
        }
        double payoff_Y_pos = std::max(std::pow(prod_pos, 1.0/d) - K, 0.0);

        double prod_neg = 1.0; 
        for(int j = 0; j < d; ++j)
        {
            prod_neg *= S_T_neg(j); 
        }
        double payoff_Y_neg = std::max(std::pow(prod_neg, 1.0/d) - K, 0.0);
        double payoff_Y = (payoff_Y_pos + payoff_Y_neg) / 2.0;

        double adjusted_payoff = payoff_X - payoff_Y;
        sum += adjusted_payoff;

        V.update(adjusted_payoff);
    }

    return std::make_pair((sum / N) * discount_factor + closed_form_val, V.variance());
}