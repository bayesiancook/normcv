#include "Random.hpp"
#include <vector>

class NormalModel   {

    public:

    double theta0;
    double tau0;
    double tau;
    int p;
    int nsite;
    vector<vector<double>> X;
    vector<double> meanx;
    vector<double> meanx2;
    vector<double> postx;
    vector<double> postv;
    vector<double> testmeanx;
    vector<double> testmeanx2;
    int testm;

    NormalModel(int inp, int innsite, double intheta0, double intau, double intau0)    {
        p = inp;
        nsite = innsite;
        theta0 = intheta0;
        tau = intau;
        tau0 = intau0;
        X.assign(nsite, vector<double>(p,0));
        meanx.assign(p,0);
        meanx2.assign(p,0);
        testmeanx.assign(p,0);
        testmeanx2.assign(p,0);
        postx.assign(p,0);
        postv.assign(p,0);

        Draw();
    }

    ~NormalModel() {}

    void Draw() {
        for (int i=0; i<nsite; i++) {
            for (int j=0; j<p; j++) {
                X[i][j] = rnd::GetRandom().sNormal()/sqrt(tau) + theta0;
            }
        }
    }

    void ComputeSuffStats(int n)  {
        vector<double> m1(p,0);
        vector<double> m2(p,0);
        for (int i=0; i<n; i++) {
            for (int j=0; j<p; j++) {
                m1[j] += X[i][j];
                m2[j] += X[i][j]*X[i][j];
            }
        }
        for (int j=0; j<p; j++) {
            meanx[j] = m1[j]/n;
            meanx2[j] = m2[j]/n;
            postx[j] = n*tau / (tau0 + n*tau) * meanx[j];
            postv[j] = n*tau / (tau0 + n*tau) * meanx2[j] - postx[j]*postx[j];
        }
    }

    void ComputeTestSuffStats(int n, int m)   {
        vector<double> m1(p,0);
        vector<double> m2(p,0);
        for (int i=n; i<n+m; i++) {
            for (int j=0; j<p; j++) {
                m1[j] += X[i][j];
                m2[j] += X[i][j]*X[i][j];
            }
        }
        for (int j=0; j<p; j++) {
            testmeanx[j] = m1[j]/m;
            testmeanx2[j] = m2[j]/m;
        }
        testm = m;
    }

    double GetLogMarginal0(int n)   {
        ComputeSuffStats(n);
        double tot = 0;
        for (int j=0; j<p; j++) {
            tot += 0.5 * (-n*log(2*Pi) + n*log(tau) - n*tau*meanx2[j]);
        }
        return tot;
    }

    double GetLogMarginal(int n)    {
        ComputeSuffStats(n);
        double tot = 0;
        for (int j=0; j<p; j++) {
            tot += 0.5 * (-n*log(2*Pi) + n*log(tau) + log(tau0/(tau0 + n*tau)) - (tau0 + n*tau)*postv[j]);
        }
        return tot;
    }

    double GetLogCV(int n, int m)   {
        return GetLogMarginal(n+m) - GetLogMarginal(n);
    }

    double GetLogCV0(int n, int m)  {
        return GetLogMarginal0(n+m) - GetLogMarginal0(n);
    }

    double GetDLogCV(int n, int m)  {
        return GetLogCV(n, m) - GetLogCV0(n, m);
    }

    double GetTestLogProb(vector<double> theta) {
        double tot = 0;
        for (int j=0; j<p; j++) {
            tot += 0.5 * testm * (-log(2*Pi) + log(tau) - tau*(testmeanx2[j] - 2*theta[j]*testmeanx[j] + theta[j]*theta[j]));
        }
        return tot;
    }

    double GetISLogCV(int n, int m, int nsample, double& var, double& ess)    {
        ComputeSuffStats(n);
        ComputeTestSuffStats(n,m);

        double posttau = tau0 + n*tau;
        double sigma = 1.0 / sqrt(posttau);

        vector<double> lnL(nsample,0);

        double max = 0;
        vector<double> theta(p,0);
        for (int i=0; i<nsample; i++)   {
            for (int j=0; j<p; j++) {
                theta[j] = sigma*rnd::GetRandom().sNormal() + postx[j];
            }
            lnL[i] = GetTestLogProb(theta);
            if ((!i) || (max < lnL[i])) {
                max = lnL[i];
            }
        }

        double tot = 0;
        double tot2 = 0;
        for (int i=0; i<nsample; i++)   {
            double tmp = exp(lnL[i] - max);
            tot += tmp;
            tot2 += tmp*tmp;
        }
        tot /= nsample;
        tot2 /= nsample;
        tot2 /= tot*tot;
        var = (tot2 - 1) / (nsample-1);
        ess = nsample/tot2;
        double ret = log(tot) + max;
        return ret;
    }

    double GetISDLogCV(int n, int m, int nsample, double& var, double& ess)   {
        return GetISLogCV(n,m,nsample,var,ess) - GetLogCV0(n,m);
    }

    double GetSteppingLogCV(int n, int m, int nsample, vector<double>& sitecv, double& meanvar, double& meaness)  {
        double tot = 0;
        double var, ess;
        meaness = 0;
        meanvar = 0;
        for (int i=0; i<m; i++)   {
            sitecv[i] = GetISLogCV(n+i, 1, nsample, var, ess);
            meanvar += var;
            meaness += ess;
            tot += sitecv[i];
        }
        meanvar /= m;
        meaness /= m;
        return tot;
    }

    double GetSteppingDLogCV(int n, int m, int nsample, double& meanvar, double& meaness) {
        vector<double> sitecv(m,0);
        return GetSteppingLogCV(n, m, nsample, sitecv, meaness, meanvar) - GetLogCV0(n,m);
    }
};

