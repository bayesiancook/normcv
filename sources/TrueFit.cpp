#include "MultiNormalModel.hpp"

double logBF(double theta, double tau0, double tau, double n)  {
    return 0.5*(log(tau0/(tau0 + n*tau)) + n*tau/(tau0 + n*tau)*(1 + n*tau*theta*theta));
}

int main (int argc, char* argv[])   {

    int p = atoi(argv[1]);
    int n = atoi(argv[2]);
    double f = atof(argv[3]);
    double thetamax = atof(argv[4]);
    double thetastep = atof(argv[5]);
    double tau = atof(argv[6]);
    double tau0 = atof(argv[7]);

    for (double theta=0; theta<thetamax; theta+=thetastep) {

        int ntrain = int(n*(1-f));
        int nvalid = int(n*f);

        double r0 = p*theta*theta;
        double r1 = p*(n*tau + tau0*tau0*theta*theta)/(tau0 + n*tau)/(tau0 + n*tau);

        double logbf = logBF(theta, tau0, tau, n) / n;
        double logcv = (logBF(theta, tau0, tau, n) - logBF(theta, tau0, tau, ntrain)) / nvalid;
        double loocv = logBF(theta, tau0, tau, n) - logBF(theta, tau0, tau, n-1);

        cout << theta << '\t' << r1 << '\t' << r0 << '\t' << r0-r1 << '\t' << logbf << '\t' << logcv << '\t' << loocv << '\n';
    }
}
