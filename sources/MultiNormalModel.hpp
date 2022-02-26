#include "Random.hpp"
#include <vector>

class NormalModel   {

    public:

    double theta0;
    double tau0;
    double tau;
    double overalpha;
    double vintra;
    double vtot;
    int p;
    int nsite;
    vector<vector<double>> X;
    vector<double> meanx;
    vector<double> meanx2;
    vector<double> empmeanx;
    vector<double> empmeanx2;
    vector<double> postx;
    vector<double> postv;
    vector<double> testmeanx;
    vector<double> testmeanx2;
    int testm;

    NormalModel(int inp, int innsite, double intheta0, double intau, double intau0, double inoveralpha = 1.0)    {
        p = inp;
        nsite = innsite;
        theta0 = intheta0;
        tau = intau;
        tau0 = intau0;
        overalpha = inoveralpha;
        X.assign(nsite, vector<double>(p,0));
        meanx.assign(p,0);
        meanx2.assign(p,0);
        empmeanx.assign(p,0);
        empmeanx2.assign(p,0);
        testmeanx.assign(p,0);
        testmeanx2.assign(p,0);
        postx.assign(p,0);
        postv.assign(p,0);

        Draw();
    }

    ~NormalModel() {}

    void DataToStream(ostream& os) {
        os << nsite << '\t' << p << '\n';
        for (int i=0; i<nsite; i++) {
            for (int j=0; j<p; j++) {
                os << X[i][j] << '\t';
            }
            os << '\n';
        }
    }

    void DataFromStream(istream& is)    {
        int nn, pp;
        is >> nn >> pp;
        if (nn != nsite)    {
            cerr << "error when reading data from file: non matching number of sites\n";
            exit(1);
        }
        if (pp != p)    {
            cerr << "error when reading data from file: non matching model dimension\n";
            exit(1);
        }

        for (int i=0; i<nsite; i++) {
            for (int j=0; j<p; j++) {
                is >> X[i][j];
            }
        }
    }


    void Draw() {
        for (int i=0; i<nsite; i++) {
            for (int j=0; j<p; j++) {
                X[i][j] = rnd::GetRandom().sNormal()/sqrt(tau*overalpha) + theta0;
            }
            /*
            double m = 1.0;
            if (overalpha)  {
                m = rnd::GetRandom().Gamma(overalpha, overalpha);
            }
            for (int j=0; j<p; j++) {
                X[i][j] = rnd::GetRandom().sNormal()/sqrt(tau/m) + theta0;
                X[i][j] = rnd::GetRandom().sNormal()/sqrt(tau*overalpha) + theta0;
            }
            */
        }
    }

    void ComputeSuffStats(int n)  {
        if (!n) {
            for (int j=0; j<p; j++) {
                meanx[j] = 0;
                meanx2[j] = 0;
                postx[j] = 0;
                postv[j] = 0;
            }
            vtot = vintra = 0;
        }
        else    {
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
            vtot = 0;
            vintra = 0;
            for (int j=0; j<p; j++) {
                vintra += meanx2[j] - meanx[j]*meanx[j];
                vtot += meanx[j] * meanx[j];
            }
            vintra /= p;
            vtot /= p;
        }
    }

    void ComputeSuffStatsSegmentPlusPoint(int n, int i0)  {
        vector<double> m1(p,0);
        vector<double> m2(p,0);
        for (int i=0; i<n; i++) {
            for (int j=0; j<p; j++) {
                m1[j] += X[i][j];
                m2[j] += X[i][j]*X[i][j];
            }
        }
        for (int j=0; j<p; j++) {
            m1[j] += X[i0][j];
            m2[j] += X[i0][j]*X[i0][j];
        }
        for (int j=0; j<p; j++) {
            meanx[j] = m1[j]/(n+1);
            meanx2[j] = m2[j]/(n+1);
            postx[j] = (n+1)*tau / (tau0 + (n+1)*tau) * meanx[j];
            postv[j] = (n+1)*tau / (tau0 + (n+1)*tau) * meanx2[j] - postx[j]*postx[j];
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

    void ComputeEmpiricalSuffStats(int n)  {
        if (!n) {
            for (int j=0; j<p; j++) {
                empmeanx[j] = 0;
                empmeanx2[j] = 0;
            }
        }
        else    {
            vector<double> m1(p,0);
            vector<double> m2(p,0);
            for (int i=0; i<n; i++) {
                for (int j=0; j<p; j++) {
                    m1[j] += X[i][j];
                    m2[j] += X[i][j]*X[i][j];
                }
            }
            for (int j=0; j<p; j++) {
                empmeanx[j] = m1[j]/n;
                empmeanx2[j] = m2[j]/n;
            }
        }
    }

    double GetSiteLogProb(const vector<double>& theta, int i)   {
        double s2 = 0;
        for (int j=0; j<p; j++) {
            s2 += (X[i][j]-theta[j])*(X[i][j]-theta[j]);
        }
        double ret = 0.5*(p*log(tau/2/Pi) - tau*s2);
        return ret;
    }

    double GetLogProb(const vector<double>& theta, vector<double>& logl)    {
        double tot = 0;
        for (int i=0; i<nsite; i++) {
            logl[i] = GetSiteLogProb(theta, i);
            tot += logl[i];
        }
        return tot;
    }

    double GetLogMarginal0(int n)   {
        ComputeSuffStats(n);
        double tot = 0;
        for (int j=0; j<p; j++) {
            tot += 0.5 * (-n*log(2*Pi) + n*log(tau) - n*tau*meanx2[j]);
        }
        return tot;
    }

    double GetLogMarginal0SegmentPlusPoint(int n, int i0)   {
        ComputeSuffStatsSegmentPlusPoint(n,i0);
        double tot = 0;
        for (int j=0; j<p; j++) {
            tot += 0.5 * (-(n+1)*log(2*Pi) + (n+1)*log(tau) - (n+1)*tau*meanx2[j]);
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

    double GetLogMarginalSegmentPlusPoint(int n, int i0)    {
        ComputeSuffStatsSegmentPlusPoint(n,i0);
        double tot = 0;
        for (int j=0; j<p; j++) {
            tot += 0.5 * (-(n+1)*log(2*Pi) + (n+1)*log(tau) + log(tau0/(tau0 + (n+1)*tau)) - (tau0 + (n+1)*tau)*postv[j]);
        }
        return tot;
    }

    double GetSiteLogCV(int n, int m)   {
        double tot = 0;
        for (int i=0; i<m; i++) {
            tot += GetLogMarginalSegmentPlusPoint(n,n+i) - GetLogMarginal(n);
        }
        return tot;
    }

    double GetSiteLogCV0(int n, int m)  {
        double tot = 0;
        for (int i=0; i<m; i++) {
            tot += GetLogMarginal0SegmentPlusPoint(n,n+i) - GetLogMarginal0(n);
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

    double GetSiteISLogCV(int n, int m, int nsample, double& totvar, double& meaness)    {
        ComputeSuffStats(n);
        // ComputeTestSuffStats(n,m);

        double posttau = tau0 + n*tau;
        double sigma = 1.0 / sqrt(posttau);

        vector<vector<double> > lnL(nsample,vector<double>(m,0));

        vector<double> theta(p,0);
        for (int i=0; i<nsample; i++)   {
            for (int j=0; j<p; j++) {
                theta[j] = sigma*rnd::GetRandom().sNormal() + postx[j];
            }
            for (int j=0; j<m; j++) {
                lnL[i][j] = GetSiteLogProb(theta, n+j);
            }
        }

        double ret = 0;
        totvar = 0;
        meaness = 0;
        for (int j=0; j<m; j++) {
            double max = 0;
            for (int i=0; i<nsample; i++)   {
                if ((!i) || (max < lnL[i][j]))  {
                    max = lnL[i][j];
                }
            }
            double tot = 0;
            double tot2 = 0;
            for (int i=0; i<nsample; i++)   {
                double tmp = exp(lnL[i][j] - max);
                tot += tmp;
                tot2 += tmp*tmp;
            }
            tot /= nsample;
            tot2 /= nsample;
            tot2 /= tot*tot;
            double var = (tot2 - 1) / (nsample-1);
            double ess = nsample/tot2;
            ret += log(tot) + max;
            totvar += var;
            meaness += ess;
        }
        meaness /= m;
        return ret;
    }

    double GetEmpiricalLogPrior(int n, double empfrac, const vector<double>& theta) {
        double emptau0 = tau0 + (1-empfrac)*nsite*tau;
        double tot = 0;
        for (int j=0; j<p; j++) {
            double emppriortheta = (1-empfrac)*nsite*tau*empmeanx[j] / emptau0;
            tot += 0.5 * (log(emptau0/2/Pi) - emptau0*(theta[j] - emppriortheta)*(theta[j]-emppriortheta));
        }
        return tot;
    }

    double GetISLogCVWithEmpPrior(int n, int m, int nsample, double& var, double& ess)  {

        /*
        double empfrac1 = ((double) n) / nsite;
        double empfrac2 = ((double) (n+m)) / nsite;
        */
        double empfrac1 = ((double) 2*n) / nsite;
        if (empfrac1 > 1.0) {
            empfrac1 = 1.0;
        }
        double empfrac2 = ((double) 2*(n+m)) / nsite;
        if (empfrac2 > 1.0) {
            empfrac2 = 1.0;
        }
        ComputeSuffStats(n);
        ComputeTestSuffStats(n,m);

        double emptau0 = tau0 + (1-empfrac1)*nsite*tau;
        double posttau = emptau0 + n*tau;
        double sigma = 1.0 / sqrt(posttau);

        vector<double> lnL(nsample,0);

        double max = 0;
        vector<double> theta(p,0);
        for (int i=0; i<nsample; i++)   {
            for (int j=0; j<p; j++) {
                double emppriortheta = (1-empfrac1)*nsite*tau*empmeanx[j] / emptau0;
                double px = (emptau0*emppriortheta + n*tau*meanx[j]) / posttau;
                theta[j] = sigma*rnd::GetRandom().sNormal() + px;
            }
            lnL[i] = GetTestLogProb(theta);
            lnL[i] += GetEmpiricalLogPrior(n+m, empfrac2, theta) - GetEmpiricalLogPrior(n, empfrac1, theta);
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

    double GetEmpiricalSteppingLogCV(int n, int m, int nsample, vector<double>& sitecv, double& meanvar, double& meaness)  {

        ComputeEmpiricalSuffStats(n+m);

        double tot = 0;
        double var, ess;
        meaness = 0;
        meanvar = 0;
        for (int i=0; i<m; i++)   {
            sitecv[i] = GetISLogCVWithEmpPrior(n+i, 1, nsample, var, ess);
            meanvar += var;
            meaness += ess;
            tot += sitecv[i];
        }
        meanvar /= m;
        meaness /= m;
        return tot;
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

    double GetRMSD0()   {
        return p * theta0 * theta0;
    }

    double GetRMSD(const vector<double>& theta)    {
        double tot = 0;
        for (int j=0; j<p; j++) {
            tot += (theta0 - theta[j]) * (theta0 - theta[j]);
        }
        return tot;
    }

    double GetML_RMSD()  {
        ComputeSuffStats(nsite);
        vector<double> theta(p,0);
        for (int j=0; j<p; j++) {
            theta[j] = meanx[j];
        }
        return GetRMSD(theta);
    }

    double GetMPL_RMSD()    {
        ComputeSuffStats(nsite);
        TunePenalization();
        vector<double> theta(p,0);
        for (int j=0; j<p; j++) {
            theta[j] = meanx[j] * nsite*tau / (tau0 + nsite*tau);
        }
        return GetRMSD(theta);
    }

    double GetAIC() {
        vector<double> sitelogl(nsite,0);
        vector<double> theta(p,0);
        ComputeSuffStats(nsite);
        for (int j=0; j<p; j++) {
            theta[j] = meanx[j];
        }
        double lnl = GetLogProb(theta, sitelogl);
        double aic = lnl - p;
        return aic;
    }

    double GetAIC0() {
        vector<double> sitelogl(nsite,0);
        vector<double> theta(p,0);
        for (int j=0; j<p; j++) {
            theta[j] = 0;
        }
        double aic = GetLogProb(theta, sitelogl);
        return aic;
    }

    double GetBIC0() {
        return GetAIC0();
    }

    double GetBIC() {
        return GetAIC() + p - 0.5*p*log(nsite);
    }

    double GetTraceIJ0() {
        ComputeSuffStats(nsite);
        double tot = 0;
        for (int j=0; j<p; j++) {
            tot += meanx2[j];
        }
        tot *= tau;
        return tot;
    }

    double GetTraceIJ() {
        ComputeSuffStats(nsite);
        double tot = 0;
        for (int j=0; j<p; j++) {
            tot += meanx2[j] - meanx[j]*meanx[j];
        }
        tot *= tau;
        return tot;
    }

    double GetTIC0()    {
        return GetAIC0();
    }

    double GetTIC() {
        return GetAIC() + p - GetTraceIJ();
    }

    double TunePenalization() {
        ComputeSuffStats(nsite);
        double lambda = vintra / nsite / vtot;
        tau0 = nsite * tau * lambda / (1-lambda);
        return 1.0 - lambda;
    }

    double GetRIC0()    {
        return GetAIC0();
    }

    // assumes penalization has already been tuned
    double GetRIC() {
        int n = nsite;
        vector<double> sitelogl(nsite,0);
        vector<double> theta(p,0);
        for (int j=0; j<p; j++) {
            theta[j] = meanx[j] * n*tau / (tau0 + n*tau);
        }
        double lnl = GetLogProb(theta, sitelogl);

        double r = p * vintra * tau;
        double T = n*tau / (tau0 + n*tau) * r;

        return lnl - T;
    }

    double GetMPLLooCV(vector<double>& sitecv)    {
        int n = nsite;
        double tot = 0;
        vector<double> theta(p,0);
        for (int i=0; i<nsite; i++) {
            for (int j=0; j<p; j++) {
                double mx1 = (n*meanx[j] - X[i][j])/(n-1);
                theta[j] = mx1*n*tau / (tau0 + n*tau);
            }
            double score = GetSiteLogProb(theta, i);
            tot += score;
        }
        return tot;
    }

    double GetMLLooCV(vector<double>& sitecv)    {

        ComputeSuffStats(nsite);

        int n = nsite;
        double tot = 0;
        vector<double> theta(p,0);
        for (int i=0; i<nsite; i++) {
            for (int j=0; j<p; j++) {
                double mx1 = (n*meanx[j] - X[i][j])/(n-1);
                theta[j] = mx1;
            }
            double score = GetSiteLogProb(theta, i);
            tot += score;
        }
        return tot;
    }

    double GetLooCV(vector<double>& sitecv)    {

        ComputeSuffStats(nsite);

        double logtotn = GetLogMarginal(nsite);

        int n = nsite;
        double tot = 0;
        for (int i=0; i<nsite; i++) {
            double sitetot = 0;
            for (int j=0; j<p; j++) {
                double mx1 = (n*meanx[j] - X[i][j])/(n-1);
                double mx2 = (n*meanx2[j] - X[i][j]*X[i][j])/(n-1);
                double px = (n-1)*tau / (tau0 + (n-1)*tau) * mx1;
                double pv = (n-1)*tau / (tau0 + (n-1)*tau) * mx2 - px*px;
                sitetot += 0.5 * (-(n-1)*log(2*Pi) + (n-1)*log(tau) + log(tau0/(tau0 + (n-1)*tau)) - (tau0 + (n-1)*tau)*pv);
            }
            sitecv[i] = logtotn - sitetot;
            tot += sitecv[i];
        }
        return tot;
    }

    double GetLooCV0(vector<double>& sitecv)   {
        vector<double> theta(p,0);
        for (int j=0; j<p; j++) {
            theta[j] = 0;
        }
        double logp = GetLogProb(theta, sitecv);
        return logp;
    }

    double GetDMLLooCV(vector<double>& sitecv) {
        vector<double> sitecv0(nsite, 0);
        double logp = GetMLLooCV(sitecv);
        double logp0 = GetLooCV0(sitecv0);
        for (int i=0; i<nsite; i++) {
            sitecv[i] -= sitecv0[i];
        }
        logp -= logp0;
        return logp;
    }

    double GetDLooCV(vector<double>& sitecv) {
        vector<double> sitecv0(nsite, 0);
        double logp = GetLooCV(sitecv);
        double logp0 = GetLooCV0(sitecv0);
        for (int i=0; i<nsite; i++) {
            sitecv[i] -= sitecv0[i];
        }
        logp -= logp0;
        return logp;
    }

    double GetLogCPO(int nsample, vector<double>& sitecv, double& meanvar, double& meaness)  {

        ComputeSuffStats(nsite);

        vector<vector<double>> logl(nsample, vector<double>(nsite,0));
        vector<double> theta(p,0);

        double posttau = tau0 + nsite*tau;
        double sigma = 1.0 / sqrt(posttau);

        for (int rep=0; rep<nsample; rep++) {
            for (int j=0; j<p; j++) {
                theta[j] = sigma*rnd::GetRandom().sNormal() + postx[j];
            }
            GetLogProb(theta, logl[rep]);
        }

        meanvar = 0;
        meaness = 0;
        double cv = 0;
        for (int i=0; i<nsite; i++) {
            double min = 0;
            for (int rep=0; rep<nsample; rep++) {
                if ((!rep) || (min > logl[rep][i])) {
                    min = logl[rep][i];
                }
            }
            double tot = 0;
            double tot2 = 0;
            for (int rep=0; rep<nsample; rep++) {
                double tmp = exp(min - logl[rep][i]);
                tot += tmp;
                tot2 += tmp*tmp;
            }
            tot /= nsample;
            tot2 /= nsample;
            tot2 /= tot*tot;
            double var = (tot2 - 1) / (nsample-1);
            double ess = nsample/tot2;
            double sitecpo = min - log(tot);

            sitecv[i] = sitecpo;
            cv += sitecpo;
            meanvar += var;
            meaness += ess;
        }
        meanvar /= nsite;
        meaness /= nsite;

        return cv;
    }

    double GetDLogCPO(int nsample, vector<double>& sitecv, double& meanvar, double& meaness)  {
        double logp = GetLogCPO(nsample, sitecv, meanvar, meaness);
        vector<double> sitecv0(nsite,0);
        double logp0 = GetLooCV0(sitecv0);
        for (int i=0; i<nsite; i++) {
            sitecv[i] -= sitecv0[i];
        }
        logp -= logp0;
        return logp;
    }

};

