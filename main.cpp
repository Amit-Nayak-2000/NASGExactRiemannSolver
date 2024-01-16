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
    PrimL[0] = 1.0; //Density
    PrimL[1] = 0.0; //Velocity
    PrimL[2] = 1.0; //Pressure

    std::vector<double> PrimR(3);
    PrimR[0] = 0.125; //Density
    PrimR[1] = 0.0; //Velocity
    PrimR[2] = 0.1; //Pressure

    double P, U;

    calcStar(P, U, PrimL, PrimR, gamma, Pinf, b);

}