#include "Random.hpp"
using namespace std;
#include <vector>
#include "rr.h"

class RRModel  {

    int full_nsite;
    vector<vector<double>> full_site_count;
    vector<vector<double>> full_site_beta;

    int nsite;
    vector<vector<double>> site_count;
    vector<vector<double>> site_beta;

    vector<double> tot_count;
    vector<double> tot_beta;

    vector<double> uni_rr;
    vector<double> ref_rr;
    vector<double> true_rr;

    public:

    RRModel(string datafile)   {
        ifstream is(datafile.c_str());
        is >> full_nsite;
        full_site_count.assign(full_nsite,vector<double>(nrr,0));
        full_site_beta.assign(full_nsite,vector<double>(nrr,0));
        tot_count.assign(nrr,0);
        tot_beta.assign(nrr,0);
        for (int i=0; i<full_nsite; i++) {
            for (int j=0; j<nrr; j++)   {
                is >> full_site_count[i][j];
                tot_count[j] += full_site_count[i][j];
            }
            for (int j=0; j<nrr; j++)   {
                is >> full_site_beta[i][j];
                tot_beta[j] += full_site_beta[i][j];
            }
        }

        uni_rr.assign(nrr,1.0);

        ref_rr.assign(nrr,0);
        true_rr.assign(nrr,0);
        double ref_mean = 0;
        double true_mean = 0;
        for (int j=0; j<nrr; j++)   {
            true_rr[j] = tot_count[j] / tot_beta[j];
            ref_rr[j] = LG_RR[j];
            true_mean += true_rr[j];
            ref_mean += ref_rr[j];
        }
        ref_mean /= nrr;
        true_mean /= nrr;
        for (int j=0; j<nrr; j++)   {
            ref_rr[j] /= ref_mean;
            // ref_rr[j] = 0.5*ref_rr[j] + 0.5*true_rr[j];
            // cout << true_rr[j] << '\t' << ref_rr[j] << '\n';
        }
        // exit(1);
    }

    void MakeRandomReplicate(int innsite)   {
        nsite = innsite;
        site_count.assign(nsite,vector<double>(nrr,0));
        site_beta.assign(nsite,vector<double>(nrr,0));
        int* subset = new int[nsite];
        rnd::GetRandom().DrawFromUrn(subset, nsite, full_nsite);
        for (int j=0; j<nrr; j++)   {
            tot_count[j] = 0;
            tot_beta[j] = 0;
        }
        for (int i=0; i<nsite; i++) {
            for (int j=0; j<nrr; j++)   {
                site_count[i][j] = full_site_count[subset[i]][j];
                site_beta[i][j] = full_site_beta[subset[i]][j];
                tot_count[j] += site_count[i][j];
                tot_beta[j] += site_beta[i][j];
            }
        }
    }

    double log_likelihood(const vector<double>& count, const vector<double>& beta, const vector<double>& rr)    {
        double tot = 0;
        for (int j=0; j<nrr; j++)   {
            if (beta[j])  {
                tot += -beta[j]*rr[j] + count[j]*log(beta[j]*rr[j]);
            }
            else    {
                if (count[j])   {
                    cerr << "error in log lik: inf\n";
                    exit(1);
                }
            }
        }
        return tot;
    }

    void estimate_rr(const vector<double>& count, const vector<double>& beta, double lambda, const vector<double>& rr0, vector<double>& rr)  {
        for (int j=0; j<nrr; j++)   {
            rr[j] = (count[j] + lambda*rr0[j]) / (beta[j] + lambda);
        }
    }

    double get_vintra() {
        double vintra = 0;
        for (int j=0; j<nrr; j++)   {
            if (tot_beta[j])    {
                vintra += tot_count[j] / tot_beta[j];
            }
        }
        return vintra;
    }

    double get_vtot() {
        double vtot = 0;
        for (int j=0; j<nrr; j++)   {
            if (tot_beta[j])    {
                double rr = tot_count[j] / tot_beta[j];
                vtot += tot_beta[j] * (rr - ref_rr[j])*(rr - ref_rr[j]);
            }
        }
        return vtot;
    }

    double get_lambda() {
        double vintra = 0;
        double vtot = 0;
        for (int j=0; j<nrr; j++)   {
            if (tot_beta[j])    {
                double rr = tot_count[j] / tot_beta[j];
                vtot += tot_beta[j] * (rr - ref_rr[j])*(rr - ref_rr[j]);
                vintra += tot_count[j] / tot_beta[j];
            }
        }
        double lambda = 10;
        if ((vtot > 0) && (vintra > 0) && (vintra < vtot))   {
            lambda = vintra / (vtot-vintra);
        }
        return lambda;
    }

    double loocv0()  {
        double mean_loocv = 0;
        vector<double> tmpcount(nrr,0);
        vector<double> tmpbeta(nrr,0);
        for (int i=0; i<nsite; i++) {
            for (int j=0; j<nrr; j++)   {
                tmpcount[j] = tot_count[j] - site_count[i][j];
                tmpbeta[j] = tot_beta[j] - site_beta[i][j];
            }
            double lnl = log_likelihood(site_count[i], site_beta[i], ref_rr);
            mean_loocv += lnl;
        }
        mean_loocv /= nsite;
        return mean_loocv;
    }

    double loocv(double lambda, const vector<double>& rr0)  {
        double mean_loocv = 0;
        vector<double> tmpcount(nrr,0);
        vector<double> tmpbeta(nrr,0);
        for (int i=0; i<nsite; i++) {
            for (int j=0; j<nrr; j++)   {
                tmpcount[j] = tot_count[j] - site_count[i][j];
                tmpbeta[j] = tot_beta[j] - site_beta[i][j];
            }
            vector<double> rr(nrr,0);
            estimate_rr(tmpcount, tmpbeta, lambda, rr0, rr);
            double lnl = log_likelihood(site_count[i], site_beta[i], rr);
            mean_loocv += lnl;
        }
        mean_loocv /= nsite;
        return mean_loocv;
    }

    double loocv(double lambda, bool uniform)   {
        if (uniform)    {
            return loocv(lambda, uni_rr);
        }
        return loocv(lambda, ref_rr);
    }

    double aic0()   {
        return log_likelihood(tot_count, tot_beta, ref_rr) / nsite;
    }

    double aic()    {
        vector<double> rr(nrr, 1.0);
        estimate_rr(tot_count, tot_beta, 1.0, uni_rr, rr);
        return (log_likelihood(tot_count, tot_beta, rr) - nrr) / nsite;
    }

    double bic0()   {
        return log_likelihood(tot_count, tot_beta, ref_rr) / nsite;
    }

    double bic()    {
        vector<double> rr(nrr, 1.0);
        estimate_rr(tot_count, tot_beta, 1.0, uni_rr, rr);
        return (log_likelihood(tot_count, tot_beta, rr) - 0.5 * nrr * log(nsite)) / nsite;
    }

    double traceIJ(double lambda, const vector<double>& rr0, const vector<double>& rr)    {

        vector<double> I(nrr,0);
        vector<double> J(nrr,0);

        for (int i=0; i<nsite; i++) {
            for (int j=0; j<nrr; j++)   {
                I[j] += (site_count[i][j] + lambda*rr0[j]/nsite) / rr[j] / rr[j];
                double d1 = -(site_beta[i][j] + lambda/nsite) + (site_count[i][j] + lambda*rr0[j]/nsite) / rr[j];
                J[j] += d1*d1;
            }
        }

        double T = 0;
        for (int j=0; j<nrr; j++)   {
            T += J[j] / I[j];
        }
        return T;
    }

    double ric0()   {
        return log_likelihood(tot_count, tot_beta, ref_rr) / nsite;
    }

    double ric(double lambda, bool uniform) {
        vector<double> rr(nrr, 1.0);
        if (uniform)    {
            estimate_rr(tot_count, tot_beta, lambda, uni_rr, rr);
            return (log_likelihood(tot_count, tot_beta, rr) - traceIJ(lambda, uni_rr, rr)) / nsite;
        }
        estimate_rr(tot_count, tot_beta, lambda, ref_rr, rr);
        return (log_likelihood(tot_count, tot_beta, rr) - traceIJ(lambda, ref_rr, rr)) / nsite;
    }

    double rmsd0()  {
        double tot = 0;
        for (int j=0; j<nrr; j++)   {
            tot += (ref_rr[j] - true_rr[j]) * (ref_rr[j] - true_rr[j]);
        }
        return tot;
    }

    double rmsd(double lambda, bool uniform)   {
        vector<double> rr(nrr, 1.0);
        if (uniform)    {
            estimate_rr(tot_count, tot_beta, lambda, uni_rr, rr);
        }
        else    {
            estimate_rr(tot_count, tot_beta, lambda, ref_rr, rr);
        }
        double tot = 0;
        for (int j=0; j<nrr; j++)   {
            tot += (rr[j] - true_rr[j]) * (rr[j] - true_rr[j]);
        }
        return tot;
    }
};



