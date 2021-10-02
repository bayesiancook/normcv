#include "MultiNormalModel.hpp"
#include "Chrono.hpp"

int main (int argc, char* argv[])   {

    int p = atoi(argv[1]);
    int n = atoi(argv[2]);
    double theta = atof(argv[3]);
    double tau = atof(argv[4]);
    double tau0 = atof(argv[5]);
    int nsample = atoi(argv[6]);
    int nrep = atoi(argv[7]);
    int persite = atoi(argv[8]);

    double meantruel = 0;
    double meanssl = 0;
    double meandssl = 0;
    double meaness = 0;
    double meanerr2 = 0;
    double meanderr2 = 0;
    double meanesterr2 = 0;
    double meandesterr2 = 0;
    double meanestbias = 0;

    // double c = rnd::GetRandom().get_studentC(nsample-1);

    Chrono ch;
    ch.Start();
    for (int rep=0; rep<nrep; rep++)    {

        NormalModel model(p, n, theta, tau, tau0);

        vector<double> truesitecv(n,0);
        double truel = model.GetLooCV(truesitecv);
        if (persite)    {
            truel /= n;
        }
        meantruel += truel;

        vector<double> sitecv1(n,0);
        vector<double> sitecv2(n,0);
        double var1, var2, ess1, ess2;
        double ssl1 = model.GetLogCPO(nsample, sitecv1, var1, ess1);
        double ssl2 = model.GetLogCPO(nsample, sitecv2, var2, ess2);
        if (persite)    {
            ssl1 /= n;
            ssl2 /= n;
            var1 /= n;
            var2 /= n;
        }
        double ssl = 0.5 * (ssl1 + ssl2);

        double var = 0;
        for (int i=0; i<n; i++) {
            var += 0.5 * (sitecv1[i] - sitecv2[i]) * (sitecv1[i] - sitecv2[i]);
        }
        if (persite)    {
            var /= n;
        }

        double estbias = 0.5 * var;
        double esterr2 = var + estbias*estbias;

        double dssl = 0.5 * (ssl1 + ssl2) - 0.5 * var;
        double desterr2 = var;

        meanssl += ssl;
        meandssl += dssl;
        meanestbias += estbias;
        meanesterr2 += esterr2;
        meandesterr2 += desterr2;
        meanerr2 += (ssl-truel)*(ssl-truel);
        meanderr2 += (dssl-truel)*(dssl-truel);

        double ess = 0.5*(ess1 + ess2);
        meaness += ess;

    }
    ch.Stop();
    double time = ch.GetTime()/nrep;

    meantruel /= nrep;
    meanssl /= nrep;
    meandssl /= nrep;
    meanestbias /= nrep;
    meaness /= nrep;
    meanerr2 /= nrep;
    meanderr2 /= nrep;
    meanesterr2 /= nrep;
    meandesterr2 /= nrep;
    double meanerr = sqrt(meanerr2);
    double meanesterr = sqrt(meanesterr2);
    double meanderr = sqrt(meanderr2);
    double meandesterr = sqrt(meandesterr2);
    double bias = meanssl - meantruel;
    double dbias = meandssl - meantruel;

    cout << meantruel << '\t' << meanssl << '\t' << bias << '\t' << meanestbias << '\t' << meanerr << '\t' << meanesterr << '\t' << meaness << '\t' << time << '\n';
    cout << meantruel << '\t' << meandssl << '\t' << dbias << '\t' << 0 << '\t' << meanderr << '\t' << meandesterr << '\t' << meaness << '\t' << time << '\n';
}
