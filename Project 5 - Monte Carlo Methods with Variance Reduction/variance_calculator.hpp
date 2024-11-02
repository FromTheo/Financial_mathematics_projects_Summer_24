#include <cmath>

class VarianceCalculator // Algorithme de Welford
{
    private : 
    int count;              // Nb. d'échantillons
    double mean;            // Moyenne courante 
    double partial_sum;     // Somme partienne pour calculer la variance 

    public:
    VarianceCalculator() : count(0), mean(0.0), partial_sum(0.0) {}

    void update(double x) 
    {
        count++;
        double delta = x - mean;
        mean += delta / count;
        double delta2 = x - mean;
        partial_sum += delta * delta2;
    }

    double variance() const 
    {
        return (count > 1) ? partial_sum / (count - 1) : 0.0;
    }

    double meanValue() const {return mean;}
};
