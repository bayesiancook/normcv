#include "Random.hpp"
using namespace std;
#include <vector>

int main(int argc, char* argv[])    {

    int n = atoi(argv[1]);
    int p = atoi(argv[2]);
    int nrep = atoi(argv[3]);
    double min = -10;
    double max = 10;

    double meanself = 0;
    double meancross = 0;

    for (int rep=0; rep<nrep; rep++)    {

        vector<double> counts(p,0.01);
        vector<double> data(n,0);

        for (int i=0; i<n; i++) {
            double x = rnd::GetRandom().sNormal();
            data[i] = x;
            int k = (x-min)/(max-min)*p;
            if (k < 0)  {
                k = 0;
            }
            if (k > p-1) {
                k = p-1;
            }
            counts[k]++;
        }

        for (int k=0; k<p; k++) {
            counts[k] /= (n+0.01*p);
        }

        double loglself = 0;
        double loglcross = 0;
        for (int i=0; i<n; i++) {
            double x = data[i];
            int k = (x-min)/(max-min)*p;
            if (k < 0)  {
                k = 0;
            }
            if (k > p-1) {
                k = p-1;
            }
            loglself += log(counts[k]);

            x = rnd::GetRandom().sNormal();
            k = (x-min)/(max-min)*p;
            if (k < 0)  {
                k = 0;
            }
            if (k > p-1) {
                k = p-1;
            }
            loglcross += log(counts[k]);
        }
        loglself /= n;
        loglcross /= n;

        meanself += loglself;
        meancross += loglcross;
    }
    meanself /= nrep;
    meancross /= nrep;
    cout << n << '\t' << meanself << '\t' << meancross << '\t' << meanself - meancross << '\n';
}

