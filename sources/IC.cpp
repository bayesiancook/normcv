#include "MultiNormalModel.hpp"

int main (int argc, char* argv[])   {

    int p = atoi(argv[1]);
    int n = atoi(argv[2]);
    double theta = atof(argv[3]);
    double tau = atof(argv[4]);
    double tau0 = atof(argv[5]);
    double overalpha = atof(argv[6]);
    int nrep = atoi(argv[7]);
    string name = argv[8];
    int persite = atoi(argv[9]);
    int relative = atoi(argv[10]);

    double meanloocv = 0;
    double meanaic = 0;
    double meanbic = 0;
    double meantic = 0;
    double meanric = 0;
    double meanrmsd0 = 0;
    double meanrmsd = 0;
    double meanmplrmsd = 0;
    double meanmplloocv = 0;
    double meanpen = 0;

    for (int rep=0; rep<nrep; rep++)    {

        NormalModel model(p, n, theta, tau, tau0, overalpha);
        if (name != "random")   {
            ostringstream s;
            s << name << rep << ".data";
            ifstream is(s.str().c_str());
            model.DataFromStream(is);
        }

        double rmsd0 = model.GetRMSD0();
        double rmsd = model.GetML_RMSD();
        meanrmsd0 += rmsd0;
        meanrmsd += rmsd;

        vector<double> sitecv(n,0);
        double loocv = model.GetMLLooCV(sitecv);
        double loocv0 = model.GetLooCV0(sitecv);
        if (relative)   {
            loocv -= loocv0;
        }
        if (persite)    {
            loocv /= n;
        }
        meanloocv += loocv;

        double aic = model.GetAIC();
        double aic0 = model.GetAIC0();
        double bic = model.GetBIC();
        double bic0 = model.GetBIC0();
        double tic = model.GetTIC();
        double tic0 = model.GetTIC0();

        if (relative)   {
            aic -= aic0;
            bic -= bic0;
            tic -= tic0;
        }
        if (persite)    {
            aic /= n;
            bic /= n;
            tic /= n;
        }

        meanaic += aic;
        meanbic += bic;
        meantic += tic;

        // penalized
        double pen = model.TunePenalization();
        meanpen += pen;
        double mplrmsd = model.GetMPL_RMSD();
        double mplloocv = model.GetMPLLooCV(sitecv);
        double ric = model.GetRIC();
        double ric0 = model.GetRIC0();

        if (relative)   {
            ric -= ric0;
            mplloocv -= loocv0; 
        }

        if (persite)    {
            ric /= n;
            mplloocv /= n;
        }

        meanmplrmsd += mplrmsd;
        meanric += ric;
        meanmplloocv += mplloocv;

    }

    meanrmsd0 /= nrep;
    meanrmsd /= nrep;

    meanloocv /= nrep;
    meanaic /= nrep;
    meanbic /= nrep;
    meantic /= nrep;

    meanpen /= nrep;
    meanmplrmsd /= nrep;
    meanmplloocv /= nrep;
    meanric /= nrep;

    cout << n << '\t' << meanrmsd0 << '\t' << meanrmsd << '\t' << meanmplrmsd << '\t' << meanloocv << '\t' << meanaic << '\t' << meantic << '\t' << meanmplloocv << '\t' << meanric << '\t' << meanbic << '\n';
}
