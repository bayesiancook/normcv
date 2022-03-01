
#include <iostream>
#include <iomanip>
#include "RRModel.hpp"

int main(int argc, char* argv[])    {

    string datafile = argv[1];
    int nmin = atoi(argv[2]);
    int nmax = atoi(argv[3]);
    int nstep = atoi(argv[4]);
    int nrep = atoi(argv[5]);

    RRModel model(datafile);

    cout << "nsite\trmsd0\trmsd\tpenrmsd\tlambda\tloocv\taic\ttic\tpenloocv\tric\tbic\n";
    exit(1);

    for (int nsite = nmin; nsite<=nmax; nsite+=nstep)    {

        double mean_lambda = 0;
        double mean_loocv_uni = 0;
        double mean_loocv_pen = 0;

        double mean_rmsd_ref = 0;
        double mean_rmsd_uni = 0;
        double mean_rmsd_pen = 0;

        double mean_aic = 0;
        double mean_bic = 0;
        double mean_tic = 0;
        double mean_ric = 0;

        for (int rep=0; rep<nrep; rep++)    {
            model.MakeRandomReplicate(nsite);

            double lambda = model.get_lambda();
            double vintra = model.get_vintra();
            double vtot = model.get_vtot();
            double vinter = vtot - vintra;
            double alpha = (vinter - vintra) / vtot;
            mean_lambda += alpha;

            mean_rmsd_ref += model.rmsd0();
            mean_rmsd_uni += model.rmsd(1.0, true);
            mean_rmsd_pen += model.rmsd(lambda, false);

            mean_loocv_uni += model.loocv(1.0, true) - model.loocv0();
            mean_loocv_pen += model.loocv(lambda, false) - model.loocv0();

            mean_aic += model.aic() - model.aic0();
            mean_bic += model.bic() - model.bic0();
            mean_tic += model.ric(1.0, true) - model.ric0();
            mean_ric += model.ric(lambda, false) - model.ric0();
        }

        mean_lambda /= nrep;

        mean_rmsd_ref /= nrep;
        mean_rmsd_uni /= nrep;
        mean_rmsd_pen /= nrep;

        mean_loocv_uni /= nrep;
        mean_loocv_pen /= nrep;

        mean_aic /= nrep;
        mean_bic /= nrep;
        mean_tic /= nrep;
        mean_ric /= nrep;

        cout << nsite << '\t' << mean_rmsd_ref << '\t' << mean_rmsd_uni << '\t' << mean_rmsd_pen << '\t' << mean_lambda << '\t' << mean_loocv_uni << '\t' << mean_aic << '\t' << mean_tic << '\t' << mean_loocv_pen << '\t' << mean_ric << '\t' << mean_bic << '\n';
        cout.flush();

    }
}

