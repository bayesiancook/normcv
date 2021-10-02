#include "MultiNormalModel.hpp"

int main (int argc, char* argv[])   {

    int p = atoi(argv[1]);
    int n = atoi(argv[2]);
    int m = atoi(argv[3]);
    double theta = atof(argv[4]);
    double tau = atof(argv[5]);
    double tau0 = atof(argv[6]);
    int nsample = atoi(argv[7]);
    int steppingnsample = atoi(argv[8]);
    int nrep = atoi(argv[9]);

    int ref = 0;

    double estmean = 0;
    double sestmean = 0;
    NormalModel model(p, n+m, theta, tau, tau0);
    double truel = model.GetLogCV(n,m);
    double refl = model.GetLogCV0(n,m);
    if (ref)    {
        truel -= refl;
    }

    double meanisl = 0;
    double varisl = 0;

    double meanssl = 0;
    double varssl = 0;

    double meandssl = 0;
    double vardssl = 0;

    double meanbias = 0;
    vector<double> sitecv1(m,0);
    vector<double> sitecv2(m,0);

    for (int rep=0; rep<nrep; rep++)    {
        double isl = model.GetISLogCV(n,m,nsample);
        double ssl1 = model.GetSteppingLogCV(n,m,steppingnsample, sitecv1);
        double ssl2 = model.GetSteppingLogCV(n,m,steppingnsample, sitecv2);
        double var = 0;
        for (int i=0; i<m; i++) {
            var += 0.5 * (sitecv1[i] - sitecv2[i]) * (sitecv1[i] - sitecv2[i]);
        }
        meanbias -= 0.5*var;

        double dssl = 0.5 * (ssl1 + ssl2) + 0.5 * var;

        if (ref)    {
            isl -= refl;
            ssl1 -= refl;
            dssl -= refl;
        }

        meanisl += isl;
        varisl += isl*isl;

        meanssl += ssl1;
        varssl += ssl1*ssl1;

        meandssl += dssl;
        vardssl += dssl*dssl;
    }

    meanisl /= nrep;
    varisl /= nrep;
    varisl -= meanisl*meanisl;
    double stdisl = sqrt(varisl / (nrep-1));
    double minisl = meanisl - 2*stdisl;
    double maxisl = meanisl + 2*stdisl;

    meanssl /= nrep;
    varssl /= nrep;
    varssl -= meanssl*meanssl;
    double stdssl = sqrt(varssl / (nrep-1));
    double minssl = meanssl - 2*stdssl;
    double maxssl = meanssl + 2*stdssl;

    meandssl /= nrep;
    vardssl /= nrep;
    vardssl -= meandssl*meandssl;
    double stddssl = sqrt(vardssl / (nrep-1));
    double mindssl = meandssl - 2*stddssl;
    double maxdssl = meandssl + 2*stddssl;

    meanbias /= nrep;

    cout << "true value   : " << truel << '\n';
    cout << "is  estimate : " << meanisl << '\t' << meanisl-truel << '\t' << stdisl << '\t' << varisl << '\n';
    cout << "ss  estimate : " << meanssl << '\t' << meanssl-truel << '\t' << meanbias << '\t' << stdssl << '\n';
    cout << "unbiased ss  : " << meandssl << '\t' << meandssl-truel << '\t' << meanbias << '\t' << stddssl << '\n';

    /*
    for (int rep=0; rep<nrep; rep++)    {
        NormalModel model(p, n+m, theta, tau, tau0);
        // double truel = model.GetDLogCV(n,m);
        // double isl = model.GetISDLogCV(n,m,nsample);
        // double sisl = model.GetSteppingDLogCV(n,m,steppingnsample);
        double truel = model.GetLogCV(n,m);
        double isl = model.GetISLogCV(n,m,nsample);
        double sisl = model.GetSteppingLogCV(n,m,steppingnsample);
        // cerr << truel << '\t' << isl << '\t' << isl-truel << '\t' << sisl << '\t' << sisl-truel << '\n';
        double d = isl - truel;
        m1 += d;
        m2 += d*d;
        double e = sisl - truel;
        n1 += e;
        n2 += e*e;
        truemean += truel;
        estmean += isl;
        sestmean += sisl;
    }
    cerr << '\n';
    truemean /= nrep;
    estmean /= nrep;
    sestmean /= nrep;
    m1 /= nrep;
    m2 /= nrep;
    double mean = m1;
    double var = (m2 - m1*m1)/(nrep-1)*nrep;
    double stdev = sqrt(var/(nrep-1));
    double min = mean - 2*stdev;
    double max = mean + 2*stdev;

    n1 /= nrep;
    n2 /= nrep;
    double smean = n1;
    double svar = (n2 - n1*n1)/(nrep-1)*nrep;
    double sstdev = sqrt(var/(nrep-1));
    double smin = smean - 2*sstdev;
    double smax = smean + 2*sstdev;
    cout << "mean true cv score : " << truemean << '\n';
    cout << '\n';
    cout << "IS mean estimate   : " << estmean << '\n';
    cout << "estimation bias    : " << mean << '\n';
    cout << "apparent stdev     : " << stdev << '\n';
    cout << "advertised CI      : " << min << '\t' << max << '\n';
    cout << '\n';
    cout << "SS mean estimate   : " << sestmean << '\n';
    cout << "estimation bias    : " << smean << '\n';
    cout << "apparent stdev     : " << sstdev << '\n';
    cout << "advertised CI      : " << smin << '\t' << smax << '\n';
    */
}
