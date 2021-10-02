#include "MultiNormalModel.hpp"
#include "Chrono.hpp"

int main (int argc, char* argv[])   {

    int p = atoi(argv[1]);
    int n = atoi(argv[2]);
    int m = atoi(argv[3]);
    double theta = atof(argv[4]);
    double tau = atof(argv[5]);
    double tau0 = atof(argv[6]);
    int nsample = atoi(argv[7]);
    int nrep = atoi(argv[8]);
    int persite = atoi(argv[9]);

    double meantruel = 0;
    double meanisl = 0;
    double meandisl = 0;
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

        double truel = model.GetLogCV(n,m);
        if (persite)    {
            truel /= m;
        }
        meantruel += truel;

        double var = 0;
        double ess = 0;
        double isl = model.GetISLogCV(n,m,nsample,var,ess);
        if (persite)    {
            isl /= m;
            var /= m;
        }

        double estbias = -0.5*var;
        double esterr2 = persite ? (var/m + estbias*estbias) : (var + estbias*estbias);

        double disl = isl + 0.5 * var;
        double desterr2 = persite ? var/m : var;

        meanisl += isl;
        meandisl += disl;
        meanerr2 += (isl-truel)*(isl-truel);
        meanderr2 += (disl-truel)*(disl-truel);
        meanesterr2 += esterr2;
        meandesterr2 += desterr2;
        meanestbias += estbias;
        meaness += ess;
    }
    ch.Stop();
    double time = ch.GetTime()/nrep;

    meantruel /= nrep;
    meanisl /= nrep;
    meandisl /= nrep;
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
    double bias = meanisl - meantruel;
    double dbias = meandisl - meantruel;

    cout << meantruel << '\t' << meanisl << '\t' << bias << '\t' << meanestbias << '\t' << meanerr << '\t' << meanesterr << '\t' << meaness << '\t' << time << '\n';
    cout << meantruel << '\t' << meandisl << '\t' << dbias << '\t' << 0 << '\t' << meanderr << '\t' << meandesterr << '\t' << meaness << '\t' << time << '\n';
}
