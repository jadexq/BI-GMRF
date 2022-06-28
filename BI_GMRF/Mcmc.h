//
//

#ifndef SIMU1_MCMC_H
#define SIMU1_MCMC_H

#include "Constant.h"
#include "Generate.h"
#include "Initial.h"

#include <Eigen/Dense>
#include <math.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include <string>
#include <cstring>
#include <cstdlib>
#include <sstream>
#include "rnorm.h"
#include "runif.h"
#include "rgamma.h"
#include "strToVec.h"

using namespace std;
using namespace Eigen;

class Mcmc
{
private:
    // parameters to be updated
    VectorXd alpha;
    VectorXd a; // lambda1
    VectorXd zetaA;
    VectorXd beta;
    VectorXd zetaB;
    double gamma;
    double b; // (not in use)
    VectorXd sigE;
    double sigD;

    // data
    VectorXd x;
    VectorXd z;
    MatrixXd m;
    MatrixXd mb;
    VectorXd y;
    MatrixXd nbh;
    VectorXd nNbh;

    // for fast computation
    // fixed
    double sumX2; // sum_{i=1}^n x_i^2, scalar
    VectorXd sumXM; // sum_{i=1}^n m_i(v) x_i, vector, length = V
    VectorXd sumXMb;
    VectorXd sumM2; // sum_{i=1}^n m_i(v)^2, vector, length = V
    VectorXd sumMb2;
    VectorXd sumYM; // sum_{i=1}^n y_i m_i(v), vector, length = V
    VectorXd sumYMb;
    // update together with beta
    VectorXd sumBetaMb; // sum_{v=1}^V beta(v) m_i(v), vector, length = n
    // related to a, b
    double sumZ2; // sum_{i=1}^n z_i^2, scalar
    VectorXd sumZM; // sum_{i=1}^n m_i(v) c1_i, vector, length = V
    VectorXd sumZMb;
    double sumZX; // sum_{i=1}^n x_i c1_i, scalar

public:
    // read data and nbh from .txt
    int readData(Constant constant, Generate generate, int rep);

    // set initial value
    int setInitial(Constant constant, Initial initial);

    // for fast computation
    int fastComp(Constant constant);

    // mcmc updates
    // update alpha, zeta_alpha, lambda1, and sigma^2_epsilon
    int updateA(Constant constant);
    // update beta and zeta_beta
    int updateB(Constant constant);
    // update gamma
    int updateGamma(Constant constant);
    // (not in use) update b
    int updateb(Constant constant);
    // update sigma^2_delta
    int updateSigD(Constant constant);

    // get values
    VectorXd getAlpha(void);
    VectorXd getZetaA(void);
    VectorXd getA(void);
    VectorXd getBeta(void);
    VectorXd getZetaB(void);
    double getGamma(void);
    double getB(void);
    VectorXd getSigE(void);
    double getSigD(void);
    VectorXd getX(void);
    VectorXd getZ(void);
    MatrixXd getM(void);
    VectorXd getY(void);

    // constructor
    Mcmc(Constant constant);
};


#endif //SIMU1_MCMC_H
