#include "MultiNormalModel.hpp"
#include "Chrono.hpp"

int main (int argc, char* argv[])   {

    int p = atoi(argv[1]);
    int nsite = atoi(argv[2]);
    double f = atof(argv[3]);
    int m = int(f*nsite);
    int n = nsite - m;

    double theta = atof(argv[4]);
    double tau = atof(argv[5]);
    double tau0 = atof(argv[6]);
    int nsample = atoi(argv[7]);
    int nrep = atoi(argv[8]);
    string name = argv[9];
    int persite = 0;

    double meantruel = 0;
    double meanssl = 0;
    double meandssl = 0;
    double meaness = 0;
    double meanerr2 = 0;
    double meanderr2 = 0;
    double meanesterr2 = 0;
    double meandesterr2 = 0;
    double meanestbias = 0;

    Chrono ch;
    ch.Start();
    for (int rep=0; rep<nrep; rep++)    {

        NormalModel model(p, n+m, theta, tau, tau0);
        if (name != "random")   {
            ostringstream s;
            s << name << rep << ".data";
            ifstream is(s.str().c_str());
            model.DataFromStream(is);
        }

        double truel = model.GetLogCV(n,m);
        if (persite)    {
            truel /= m;
        }
        meantruel += truel;

        vector<double> sitecv1(m,0);
        vector<double> sitecv2(m,0);
        double var1, var2, ess1, ess2;
        double ssl1 = model.GetSteppingLogCV(n,m,nsample, sitecv1, var1, ess1);
        double ssl2 = model.GetSteppingLogCV(n,m,nsample, sitecv2, var2, ess2);
        if (persite)    {
            ssl1 /= m;
            ssl2 /= m;
            var1 /= m;
            var2 /= m;
        }

        double ssl = 0.5*(ssl1 + ssl2);
        double var = 0;
        for (int i=0; i<m; i++) {
            var += 0.5 * (sitecv1[i] - sitecv2[i]) * (sitecv1[i] - sitecv2[i]);
        }
        if (persite)    {
            var /= m;
        }

        double estbias = -0.5 * var;
        double esterr2 = persite ? (var/m + estbias*estbias) : (var + estbias*estbias);

        double dssl = 0.5 * (ssl1 + ssl2) + 0.5 * var;
        double desterr2 = persite ? var/m : var;

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

    cout << "ss" << '\t' << p << '\t' << nsample << '\t' << meaness << '\t' << meantruel << '\t' << meanssl << '\t' << bias << '\t' << meanestbias << '\t' << meanerr << '\t' << meanesterr << '\t' << meanderr << '\t' << meandesterr << '\t' << time << '\n';
}

