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
    double meandssl = 0;
    double meaness = 0;
    double meanerr2 = 0;
    double meanderr2 = 0;
    double meanesterr2 = 0;
    double meandesterr2 = 0;
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
        cerr << var << '\t' << m*var1 << '\t' << m*var2 << '\t' << ess1 << '\t' << ess2 << '\n';
        double estbias = -0.5 * var;
        double esterr2 = var + estbias*estbias;

        double dssl = 0.5 * (ssl1 + ssl2) + 0.5 * var;
        double desterr2 = var;

        meanssl += ssl1;
        meandssl += dssl;
        meanestbias += estbias;
        meanesterr2 += esterr2;
        meandesterr2 += desterr2;
        meanerr2 += (ssl1-truel)*(ssl1-truel);
        meanderr2 += (dssl-truel)*(dssl-truel);
        meaness += ess1;

    }

    meantruel /= nrep;
    meanssl /= nrep;
    meandssl /= nrep;
    meanestbias /= nrep;
    meaness /= nrep;
    meanerr2 /= nrep;
    meanderr2 /= nrep;
    meandesterr2 /= nrep;
    double meanerr = sqrt(meanerr2);
    double meanesterr = sqrt(meanesterr2);
    double meanderr = sqrt(meanderr2);
    double meandesterr = sqrt(meandesterr2);
    double bias = meanssl - meantruel;
    double dbias = meandssl - meantruel;

    cout << meantruel << '\t' << meanssl << '\t' << bias << '\t' << meanestbias << '\t' << meanerr << '\t' << meanesterr << '\t' << meaness << '\n';
    cout << meantruel << '\t' << meandssl << '\t' << dbias << '\t' << 0 << '\t' << meanderr << '\t' << meandesterr << '\t' << meaness << '\n';
}
