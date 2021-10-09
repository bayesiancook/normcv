#include "MultiNormalModel.hpp"
#include "Chrono.hpp"
#include <sstream>

int main (int argc, char* argv[])   {

    int p = atoi(argv[1]);
    int n = atoi(argv[2]);
    double theta = atof(argv[4]);
    double tau = atof(argv[5]);
    int nrep = atoi(argv[6]);
    string name = argv[7];

    // not used
    double tau0 = 1.0;

    for (int rep=0; rep<nrep; rep++)    {

        cerr << rep +1 << '/' << nrep << '\n';
        NormalModel model(p, n, theta, tau, tau0);
        ostringstream s;
        s << name << rep << ".data";
        ofstream os(s.str().c_str());
        model.DataToStream(os);
    }
}
