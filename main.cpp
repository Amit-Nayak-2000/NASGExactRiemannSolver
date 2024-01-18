#include "NASGExactRiemann.hpp"
#include <vector>
#include <iostream>
#include <fstream>

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

    std::cout << "P*: " << P << ", U*: " << U << std::endl;

    double t = 0.25;
    std::vector<double> x = linspace(-0.9, 1, 99);
    std::vector<double> Prim(3);

    std::ofstream outputFile("output" + std::to_string(t) + ".csv");
    outputFile << "\"x\",\"Density\",\"Velocity\",\"Pressure\"\n";

    for(int i=0; i < 99; i++){
        sample(Prim, P, U, x[i]/t, PrimL, PrimR, Pinf, b, gamma);
        outputFile << x[i] << "," << Prim[0] << "," << Prim[1] << "," << Prim[2] << "\n";
    }

    outputFile.close();

    return 0;

}
