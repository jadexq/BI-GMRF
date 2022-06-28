//
//

// steps:
// 0. Initial initial(constant)

#include "Initial.h"

using namespace std;
using namespace Eigen;

Initial::Initial(Constant constant)
        : alpha(constant.vv),
          a(constant.vv),
          zetaA(constant.vv),
          beta(constant.vv),
          zetaB(constant.vv),
          sigE(constant.vv)
{
    alpha = VectorXd::Constant(constant.vv, 0); // set initial value
    a = VectorXd::Constant(constant.vv, 0);
    zetaA = VectorXd::Constant(constant.vv, 0);
    beta = VectorXd::Constant(constant.vv, 0);
    // beta = constant.beta; // true value as initial
    zetaB = VectorXd::Constant(constant.vv, 0);
    gamma = 0;
    b = 0;
    sigE = VectorXd::Constant(constant.vv, 0.1);
    sigD = constant.sigDFix;
}

// set true value (read from .txt) as initial value for beta
int Initial::setTrueBeta(void)
{
    // read beta true value, col by col vectorized
    // two tmp variables
    string strB; // store a line as str
    vector<double> vecB; // store a line as vec
    // open file
    ifstream fileB;
    fileB.open("../data/betaVec_true.txt");
    // assign
    int i = 0;
    while (getline(fileB, strB))
    {
        strToVec(strB, vecB);
        beta(i) = vecB[0];
        i++;
    } // end while i
    fileB.close();

    //cout << "betaInit2248=" << beta(2247) << endl;

    return 0;
}

