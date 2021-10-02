#include "MultiNormalModel.hpp"

int main (int argc, char* argv[])   {

    int p = atoi(argv[1]);
    int n = atoi(argv[2]);
    int m = atoi(argv[3]);
    double theta = atof(argv[4]);
    double tau = atof(argv[5]);
    double tau0 = atof(argv[6]);
    int nsample = atoi(argv[7]);
    int nrep = atoi(argv[8]);
    int ref = atoi(argv[9]);

    double meantruel = 0;
    double meanssl = 0;
    double meanstdev = 0;
    double meaness = 0;
    double meanerr2 = 0;
    double meanestbias = 0;

    for (int rep=0; rep<nrep; rep++)    {

        NormalModel model(p, n+m, theta, tau, tau0);

        double refl = model.GetLogCV0(n,m);

        double truel = model.GetLogCV(n,m);
        if (ref)    {
            truel -= refl;
        }

        meantruel += truel;

        vector<double> sitecv1(m,0);
        vector<double> sitecv2(m,0);
        double var1, var2, ess1, ess2;
        double ssl1 = model.GetSteppingLogCV(n,m,nsample, sitecv1, var1, ess1);
        double ssl2 = model.GetSteppingLogCV(n,m,nsample, sitecv2, var2, ess2);

        double var = 0;
        for (int i=0; i<m; i++) {
            var += 0.5 * (sitecv1[i] - sitecv2[i]) * (sitecv1[i] - sitecv2[i]);
        }
        double estbias = -0.5 * var;
        double stdev = sqrt(var);
        meanestbias += estbias;
        meanstdev += stdev;

        // double dssl = 0.5 * (ssl1 + ssl2) + 0.5 * var;

        meanerr2 += (ssl1-truel)*(ssl1-truel);
        meaness += ess1;
        meanssl += ssl1;
    }

    meantruel /= nrep;
    meanssl /= nrep;
    meanstdev /= nrep;
    meanestbias /= nrep;
    meaness /= nrep;
    meanerr2 /= nrep;
    double meanerr = sqrt(meanerr2);
    double bias = meanssl - meantruel;

    cout << meantruel << '\t' << meanssl << '\t' << bias << '\t' << meanestbias << '\t' << meanerr << '\t' << meanstdev << '\t' << meaness << '\n';
}
