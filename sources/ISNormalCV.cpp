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
    double meanisl = 0;
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

        double var = 0;
        double ess = 0;
        double isl = model.GetISLogCV(n,m,nsample,var,ess);
        if (ref)    {
            isl -= refl;
        }
        double stdev = sqrt(var);
        double estbias = -0.5*var;

        /*
        double mean = 0;
        double var = 0;
        for (int subrep=0; subrep<subnrep; subrep++)    {

            double isl = model.GetISLogCV(n,m,nsample);
            if (ref)    {
                isl -= refl;
            }

            mean += isl;
            var += isl*isl;

            meanerr2 += (isl-truel)*(isl-truel);
        }

        mean /= subnrep;
        var /= subnrep;
        var -= mean*mean;
        double stdev = sqrt(var/(subnrep-1)*subnrep);
        */

        meanerr2 += (isl-truel)*(isl-truel);
        meanstdev += stdev;
        meanestbias += estbias;
        meaness += ess;
        meanisl += isl;
    }

    meantruel /= nrep;
    meanisl /= nrep;
    meanstdev /= nrep;
    meanestbias /= nrep;
    meaness /= nrep;
    meanerr2 /= nrep;
    double meanerr = sqrt(meanerr2);
    double bias = meanisl - meantruel;

    cout << meantruel << '\t' << meanisl << '\t' << bias << '\t' << meanestbias << '\t' << meanerr << '\t' << meanstdev << '\t' << meaness << '\n';
}
