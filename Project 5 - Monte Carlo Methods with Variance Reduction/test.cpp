#include "pricer.hpp"
#include <chrono> 
using namespace std; 
using namespace Eigen; 

int main()
{
    // --- DÉCLARATION DES PARAMÈTRES --- 
    const int d = 10; 
    const double rho = 0.5; 
    Matrix<double, d, d> correlationMatrix = generateCorrelationMatrix<d>(rho); 
    Vector<double, d> weights;
    weights.setConstant(1.0/d); 
    Matrix<double, d, 1> s0;
    Matrix<double, d, 1> sigma; 
    s0.setConstant(100.0); 
    sigma.setConstant(0.3);

    PriceIndex<d> I(1.0, 0.1, s0, sigma, correlationMatrix, weights); 
    
    // -- TEST DE LA METHODE simulateAssetPrices(G, Z) --- 
    Matrix<double, d, 1> Z; 
    normal_distribution<double> norm(0.0,1.0); 
    mt19937_64 G(time(nullptr)); 
    for(int i = 0; i < d; ++i) {Z(i) = norm(G);}
    Matrix<double, d, 1> S_T = I.simulateAssetPrices(Z); 
    double I_T = weights.dot(S_T);
    cout << "Valeur de l'indice : " << I_T << endl; 


    const vector<int> K_values{80, 90, 100, 110, 120}; 
    const int N = 100000; 

    cout << "Monte Carlo avec : " << N << " simulations : " << endl;


    // 3. Méthode de Monte Carlo Classique 
    auto t1 = std::chrono::system_clock::now();
    for(const int& K : K_values)
    {
        std::pair<double, double> val = MonteCarlo_CallBasket<d>(I, K, G, N);  
        cout << "K = " << K << " : " << val.first << endl; 
        cout << "De variance " << val.second << endl;
    }
    auto t2 = std::chrono::system_clock::now();
    std::chrono::duration<double> diff = t2-t1;
    std::cout << "Temps d'éxécution :" << diff.count() << "s." << std::endl;
    cout << endl; 
    cout << endl;


    // 4.1 Méthode avec variables antithétiques
    cout << "VARIABLE ANTITHETIQUE : " << endl;
    auto t3 = std::chrono::system_clock::now(); 
    for(const int& K : K_values)
    {
        std::pair<double, double> val = MonteCarlo_CallBasket_Antithetic<d>(I, K, G, N);  
        cout << "K = " << K << " : " << val.first << endl; 
        cout << "De variance " << val.second << endl;
    }
    auto t4 = std::chrono::system_clock::now();
    std::chrono::duration<double> diff2 = t4-t3;
    std::cout << "Temps d'éxécution :" << diff2.count() << "s." << std::endl;
    cout << endl; 
    cout << endl;

    // 4.2 Méthode avec variable de controle
    cout << "VARIABLE DE CONTROLE : " << endl;
    auto t5 = std::chrono::system_clock::now(); 
    cout << "closed_formula = " << closed_formula(I, 100) << endl; 
    for(const int& K : K_values)
    {
        std::pair<double, double> val = MonteCarlo_CallBasket_ControlVar<d>(I, K, G, N);  
        cout << "K = " << K << " : " << val.first << endl; 
        cout << "De variance " << val.second << endl;

    }
    auto t6 = std::chrono::system_clock::now();
    std::chrono::duration<double> diff3 = t6-t5;
    std::cout << "Temps d'éxécution :" << diff3.count() << "s." << std::endl;
    cout << endl; 
    cout << endl;

    // 4.3 Mélange des deux 
    cout << "MELANGE DES DEUX : " << endl;
    auto t7 = std::chrono::system_clock::now();
    for(const int& K : K_values)
    {
        std::pair<double, double>  val = MonteCarlo_CallBasket_ControlVar_Antithetic<d>(I, K, G, N);  
        cout << "K = " << K << " : " << val.first << endl; 
        cout << "De variance " << val.second << endl;
    }
    auto t8 = std::chrono::system_clock::now();
    std::chrono::duration<double> diff4 = t8-t7;
    std::cout << "Temps d'éxécution :" << diff4.count() << "s." << std::endl;
     
    return 0; 
}