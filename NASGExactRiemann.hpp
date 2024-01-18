#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>

void evalPressure(double &F, double &FD, double P2, double P1, double spv1, double gamma, double Pinf, double b){
    if (P2 <= P1){ //Rarefaction
        F = (2*std::sqrt(gamma*(P1+Pinf)*(spv1-b))/(gamma-1))*(std::pow((P2 + Pinf)/(P1 + Pinf), (gamma-1)/(2*gamma)) - 1);
        FD = (std::sqrt(gamma*(P1+Pinf)*(spv1-b))/(gamma*P1))*std::pow((P2+Pinf)/(P1+Pinf), -(gamma+1)/(2*gamma));
    }
    else{ //Shock
        F = (P2 - P1)*std::sqrt(2.0)/std::sqrt(((P1 + 2*Pinf + P2)*gamma - P1 + P2)/(spv1-b)); 
        double NUM = ((3*P1 + P2 + 4*Pinf)*gamma - P1 + P2)*std::sqrt(2);
        double DENOM = ((2*P1 + 2*P2 + 4*Pinf)*gamma - 2*P1 + 2*P2)*std::sqrt(( (-P1 -P2 -2*Pinf)*gamma + P1 - P2 )/(b-spv1));
        FD = NUM/DENOM;
    }
}

void calcStar(double &P, double &U, const std::vector<double> &PrimL, const std::vector<double> &PrimR, double gamma, double Pinf, double b){

    double tol = 1e-6;
    int Niter = 30;

    double Udiff = PrimR[1] - PrimL[1];
    double spvL = 1/PrimL[0];
    double spvR = 1/PrimR[0];
    double Ppv = 0.5*(PrimL[2]+PrimR[2]) - 0.125*(Udiff)*(PrimL[0]+PrimR[0])*(std::sqrt(gamma*(PrimL[2] + Pinf)*(spvL - b)) + std::sqrt(gamma*(PrimR[2] + Pinf)*(spvR - b)));
    double Pold = std::max(tol, Ppv);
    double FL, FLD, FR, FRD, change;

    for(int i = 0; i < Niter; i++){
        evalPressure(FL, FLD, Pold, PrimL[2], spvL, gamma, Pinf, b);
        evalPressure(FR, FRD, Pold, PrimR[2], spvR, gamma, Pinf, b);



        P = Pold - (FL + FR + Udiff)/(FLD + FRD);

        if(P < 0.0){
            std::cout << "Negative Pressure Generated at " << i+1 << " iterations." << std::endl;
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
//0 rho, 1 U, 2 P
void sample(std::vector<double> &Prim, double Pstar, double Ustar, double S, const std::vector<double> &PrimL, const std::vector<double> &PrimR, double Pinf, double b, double gamma){
    double spvL = 1/PrimL[0];
    double spvR = 1/PrimR[0];

    double CL = std::sqrt(gamma*(PrimL[2] + Pinf)*(spvL - b));
    double CR = std::sqrt(gamma*(PrimR[2] + Pinf)*(spvR - b));

    if(S <= Ustar){ //Left side of contact surface
        if(Pstar > PrimL[2]){ //Left Shock
            double SL = PrimL[1] - CL*std::sqrt( ((gamma+1)/(2*gamma))*((Pstar + Pinf)/(PrimL[2] + Pinf)) + (gamma-1)/(2*gamma) );
            if(S <= SL){ //Primitive Vector is left state
                Prim[0] = PrimL[0];
                Prim[1] = PrimL[1];
                Prim[2] = PrimL[2];
            }
            else{ //Primitive Vector is * state
                Prim[0] = PrimL[0] / (((1 - pow(((gamma - 1)/(gamma + 1)),2)) / ((Pstar + Pinf)/(PrimL[2] + Pinf) + (gamma - 1)/(gamma + 1))) + (gamma - 1)/(gamma + 1));
                Prim[1] = Ustar;
                Prim[2] = Pstar;
            }
        }
        else{ //Left Fan
            if(S <= (PrimL[1] - CL)){ //head of fan (before)
                Prim[0] = PrimL[0];
                Prim[1] = PrimL[1];
                Prim[2] = PrimL[2];
            }
            else{
                double CLstar = CL*std::pow(Pstar/PrimL[2], (gamma-1)/(2*gamma));
                if(S >= (Ustar - CLstar)){ //tail of fan (after)
                    Prim[0] = 1/( ((std::pow((CL + (PrimL[1] - Ustar)*(gamma - 1)/2), 2)) / (gamma*(Pstar + Pinf)) ) + b );
                    Prim[1] = Ustar;
                    Prim[2] = Pstar;
                }
                else{ //Inside Fan
                std::cout << "H" << std::endl;
                    double Cratio = (2/(gamma - 1) - (PrimL[1] + S)/CL)/(2/(gamma-1) - 1);
                    Prim[0] = PrimL[0]*std::pow(Cratio, 2/(gamma-1));
                    Prim[1] = (2*(CL + S)/(gamma - 1) - PrimL[1]) / (2/(gamma - 1) - 1);
                    Prim[2] = (PrimL[2] + Pinf)*std::pow(Cratio, (2*gamma)/(gamma-1)) - Pinf;
                }
            }
        }
    }
    else{ //Right side of contact surface
        if(Pstar > PrimR[2]){ //Right Shock
            double SR = PrimR[1] + CR*std::sqrt( ((gamma+1)/(2*gamma))*((Pstar + Pinf)/(PrimR[2] + Pinf)) + (gamma-1)/(2*gamma) );
            if(S > SR){ //State is Right State
            // std::cout << "here1" << std::endl;
                Prim[0] = PrimR[0];
                Prim[1] = PrimR[1];
                Prim[2] = PrimR[2];
            }
            else{ //State is * State
            // std::cout << "here2" << std::endl;
                Prim[0] = PrimR[0] / (((1 - std::pow((gamma - 1)/(gamma + 1),2)) / ((Pstar + Pinf)/(PrimR[2] + Pinf) + (gamma - 1)/(gamma + 1))) + (gamma - 1)/(gamma + 1));
                Prim[1] = Ustar;
                Prim[2] = Pstar;
            }        
        }
        else{ //Right Fan
            // std::cout << "here3" << std::endl;
            double Cratio = (2/(gamma - 1) + (PrimR[1] - S)/CR)/(2/(gamma-1) - 1);
            Prim[0] = PrimR[0]*std::pow(Cratio, 2/(gamma-1));
            Prim[1] = (2*(CR - S)/(gamma - 1) + PrimR[1]) / (-2/(gamma - 1) + 1);
            Prim[2] = (PrimR[2] + Pinf)*std::pow(Cratio, (2*gamma)/(gamma-1)) - Pinf;
        }
    }
}

std::vector<double> linspace(double a, double b, int N){
    double dx = (b-a)/(N-1);
    std::vector<double> result(N);

    for(int i = 0; i < N; i++){
        result[i] = a + i*dx;
    }

    return result;
}