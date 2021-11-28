#include "MultiNormalModel.hpp"

double logBF(double theta, double tau0, double tau, double n)  {
    return 0.5*(log(tau0/(tau0 + n*tau)) + n*tau/(tau0 + n*tau)*(1 + n*tau*theta*theta));
}

int main (int argc, char* argv[])   {

    int p = atoi(argv[1]);
    int nmin = atoi(argv[2]);
    int nmax = atoi(argv[3]);
    double nstep = atof(argv[4]);
    double f = atof(argv[5]);
    double theta = atof(argv[6]);
    double tau = atof(argv[7]);
    double tau0 = atof(argv[8]);

    cout << "n\tR1\tR0\tDR01\tlogBF\tlogCV\tlogSiteCV\tlogLOO\n";
    for (int n = nmin; n<=nmax; n*=nstep)   {

        int ntrain = int(n*(1-f));
        int nvalid = int(n*f);

        double r0 = p*theta*theta;
        double r1 = p*(n*tau + tau0*tau0*theta*theta)/(tau0 + n*tau)/(tau0 + n*tau);

        double logbf = logBF(theta, tau0, tau, n) / n;
        double logcv = (logBF(theta, tau0, tau, n) - logBF(theta, tau0, tau, ntrain)) / nvalid;
        double loocv = logBF(theta, tau0, tau, n) - logBF(theta, tau0, tau, n-1);
        double logsitecv = logBF(theta, tau0, tau, ntrain+1) - logBF(theta, tau0, tau, ntrain);

        cout << n << '\t' << r1 << '\t' << r0 << '\t' << r0-r1 << '\t' << logbf << '\t' << logcv << '\t' << logsitecv << '\t' << loocv << '\n';
    }
}
