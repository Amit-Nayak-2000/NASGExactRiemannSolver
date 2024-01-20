#include "NASGExactRiemann.hpp"
#include <vector>
#include <iostream>
#include <fstream>

int main(int argc, char *argv[]){
    // NASG Variables
    double gamma = 1.4;
    double Pinf = 0.0;
    double b = 0.0;
    double q = 0.0;

    // Primitive Vectors: density, velocity, pressure, internal energy
    std::vector<double> PrimL(4);
    PrimL[0] = 1.0; //Density
    PrimL[1] = 0.0; //Velocity
    PrimL[2] = 0.01; //Pressure
    //Calculated from above params
    PrimL[3] = (PrimL[2] + gamma*Pinf)*(1/PrimL[0] - b)/(gamma - 1) + q;

    std::vector<double> PrimR(4);
    PrimR[0] = 1.0; //Density
    PrimR[1] = 0.0; //Velocity
    PrimR[2] = 100; //Pressure
    PrimR[3] = (PrimR[2] + gamma*Pinf)*(1/PrimR[0] - b)/(gamma - 1) + q;

    double P, U;

    calcStar(P, U, PrimL, PrimR, gamma, Pinf, b);

    std::cout << "P*: " << P << ", U*: " << U << std::endl;

    double t = 0.035;
    int datasize = 200;
    std::vector<double> x = linspace(-1, 2, datasize);
    std::vector<double> Prim(3);

    std::ofstream outputFile("output" + std::to_string(t) + ".csv");
    outputFile << "\"x\",\"Density\",\"Velocity\",\"Pressure\",\"Internal Energy\"\n";

    for(int i=0; i < datasize; i++){
        sample(Prim, P, U, x[i]/t, PrimL, PrimR, Pinf, b, gamma, q);
        outputFile << x[i] << "," << Prim[0] << "," << Prim[1] << "," << Prim[2] << "," << Prim[3] << "\n";
    }

    outputFile.close();

    return 0;

}
