#include <iostream>
#include <cmath>
#include <vector>

void evalPressure(double &F, double &FD, double P2, double P1, double spv1, double gamma, double Pinf, double b){
    if (P2 <= P1){ //Rarefaction
        F = (2*std::sqrt(gamma*(P1+Pinf)*(spv1-b))/(gamma-1))*(std::pow((P2 + Pinf)/(P1 + Pinf), (gamma-1)/(2*gamma)) - 1);
        FD = std::sqrt(gamma*(P1+Pinf)*(spv1-b))*std::pow((P2 + Pinf)/(P1 + Pinf), (gamma-1)/(2*gamma))/(gamma*(P2+Pinf));
    }
    else{ //Shock
        F = (P2 - P1)*std::sqrt(2.0)/std::sqrt(((P1 + 2*Pinf + P2)*gamma - P1 + P2)/(spv1-b)); 
        FD = -(P2 - P1)*std::sqrt(2.0)*(gamma + 1)/(2*(spv1-b)*std::pow(((P1 + 2*Pinf + P2)*gamma - P1 + P2)/(spv1-b), 1.5)); 
    }
}

void calcStar(double &P, double &U, const std::vector<double> &PrimL, const std::vector<double> &PrimR, double gamma, double Pinf, double b){

    double tol = 1e-6;
    int Niter = 30;

    double Pstart = 0.55;
    double Pold = Pstart;
    double Udiff = PrimR[1] - PrimL[1];
    double spvL = 1/PrimL[0];
    double spvR = 1/PrimR[0];
    double FL, FLD, FR, FRD, change;

    for(int i = 0; i < Niter; i++){
        evalPressure(FL, FLD, Pold, PrimL[2], spvL, gamma, Pinf, b);
        evalPressure(FR, FRD, Pold, PrimR[2], spvR, gamma, Pinf, b);

        std::cout << FL << " " << FR << " " << FLD << " " << FRD << std::endl;

        P = Pold - (FL + FR + Udiff)/(FLD + FRD);

        if(P < 0.0){
            std::cout << "Negative Pressure Generated" << std::endl;
        }

        change = 2.0*std::abs((P-Pold)/(P+Pold));
        // std::cout << change << std::endl;

        if(change <= tol){
            break;
            std::cout << "Converged" << std::endl;
        }
        else{
            Pold = P;
        }

        
    }

    U = 0.5*(PrimL[1] + PrimR[1] + FR - FL);
    

    
}