#include "NASGExactRiemann.hpp"
#include <vector>
#include <iostream>

int main(int argc, char *argv[]){
    // NASG Variables
    double gamma = 1.4;
    double Pinf = 0.0;
    double b = 0.0;

    // Primitive Vectors: density, velocity, pressure
    std::vector<double> PrimL(3);
    PrimL[0] = 5.99924; //Density
    PrimL[1] = 19.5975; //Velocity
    PrimL[2] = 460.894; //Pressure

    std::vector<double> PrimR(3);
    PrimR[0] = 5.99242; //Density
    PrimR[1] = -6.19633; //Velocity
    PrimR[2] = 46.0950; //Pressure

    double P, U;

    calcStar(P, U, PrimL, PrimR, gamma, Pinf, b);

    std::cout << "P*: " << P << ", U*: " << U << std::endl;

}