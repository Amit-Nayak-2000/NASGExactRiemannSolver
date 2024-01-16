#include <iostream>
#include <cmath>
#include <vector>

void evalPressure(double &F, double &FD, double P2, double P1, double spv1, double gamma, double Pinf, double b){
    if (P2 <= P1){ //Rarefaction//OG funcs work, just derivatives needed
        F = (2*std::sqrt(gamma*(P1+Pinf)*(spv1-b))/(gamma-1))*(std::pow((P2 + Pinf)/(P1 + Pinf), (gamma-1)/(2*gamma)) - 1);
        // F = (2*std::sqrt(gamma*P1*spv1)/(gamma-1))*(std::pow((P2/P1), (gamma -1)/(2*gamma)) -1); //test OG
        // FD = spv1*std::sqrt(gamma*P1*spv1)*std::pow(P2/P1, -(gamma+1)/(2*gamma)); //test derivative
        FD = (std::sqrt(gamma*(P1+Pinf)*(spv1-b))/(gamma*P1))*std::pow((P2+Pinf)/(P1+Pinf), -(gamma+1)/(2*gamma));
    }
    else{ //Shock
        F = (P2 - P1)*std::sqrt(2.0)/std::sqrt(((P1 + 2*Pinf + P2)*gamma - P1 + P2)/(spv1-b)); 
        // double A = (2*spv1)/(gamma + 1);
        // double B = (gamma-1)*P1/(gamma+1);
        // F = (P2 - P1)*std::sqrt(A/(P2 + B)); //test OG
        double NUM = ((3*P1 + P2 + 4*Pinf)*gamma - P1 + P2)*std::sqrt(2);
        double DENOM = ((2*P1 + 2*P2 + 4*Pinf)*gamma - 2*P1 + 2*P2)*std::sqrt(( (-P1 -P2 -2*Pinf)*gamma + P1 - P2 )/(b-spv1));
        // FD = std::sqrt(A/(P2 + B))*(1- (P2-P1)/(2*(B+P2))); //test derivative
        FD = NUM/DENOM;
    }
}

void calcStar(double &P, double &U, const std::vector<double> &PrimL, const std::vector<double> &PrimR, double gamma, double Pinf, double b){

    double tol = 1e-6;
    int Niter = 30;

    double Pstart = 253.494;
    double Pold = Pstart;
    double Udiff = PrimR[1] - PrimL[1];
    double spvL = 1/PrimL[0];
    double spvR = 1/PrimR[0];
    double FL, FLD, FR, FRD, change;

    for(int i = 0; i < Niter; i++){
        evalPressure(FL, FLD, Pold, PrimL[2], spvL, gamma, Pinf, b);
        evalPressure(FR, FRD, Pold, PrimR[2], spvR, gamma, Pinf, b);



        P = Pold - (FL + FR + Udiff)/(FLD + FRD);

        if(P < 0.0){
            std::cout << "Negative Pressure Generated" << std::endl;
            return;
        }

        change = 2.0*std::abs((P-Pold)/(P+Pold));

        if(change <= tol){
            std::cout << "Converged with " << i+1 << " iterations." << std::endl;
            break;
        }
        else{
            Pold = P;
        }

        
    }

    U = 0.5*(PrimL[1] + PrimR[1] + FR - FL);
    

    
}