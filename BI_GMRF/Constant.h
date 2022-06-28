//
//

#ifndef SIMU1_CONSTANT_H
#define SIMU1_CONSTANT_H

#include <iostream>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

class Constant
{
public:
    // *** setting: 0-simulation, 1-real data/simulation where the data is read from existing .txt files
    int setting = 1;
    // 0-do not apply Cv, 1-appply Cv=sqrt(V)
    int applyCv = 0;

    // *** mcmc
    int nIter = 2000; // #mcmc iterations <<
    int nBurn = 1000; // #mcmc burn-in <<
    int nRep = 1; // #replication <<
    int nPrint = 500; // print "iter" <<

    // *** data
    int n = 500; // sample size <<
    int v1 = 64; // #row of 2D images <<
    int v2 = 64; // #col of 2D images <<
    int v3 = 0;
    int vv = v1*v2;
    double cv = 1/sqrt(vv); // (not in use) the const for sum beta(v)*m(v)
    int mNnbh = 4; // maximum #nbh, =4 for 2d (=6 for 3D)
    double pi = 3.14159;

    // *** (not in use) if generate data use class `Generate`, parameter true values
    VectorXd alpha; // alpha
    double xxx1 = 0.3; // scale of alpha
    VectorXd a; // lambda1 for c1_i
    double yyy1 = 0.2; // scale of lambda1
    VectorXd beta; // beta
    double xxx2 = 0.3; // scale of beta
    double b = 0.2; // lambda2 for c2_i
    double gamma = -0.2; // gamma
    VectorXd sigE; // sig^2_epsilon, = 0.1 (specified in .cpp)
    double sigD = 0.1; // sig^2_delta

    // *** hyper-parameters
    // * binary Ising (zeta^alpha and zeta^beta)
    // p(zeta) ~ exp[ a*zeta + sum_v sum_v' b*I(zeta_v = zeta_v') ]
    // zeta_alpha
    double aIA = -4; // a_{I}^{alpha} <<
    double bIA = 0.55; // b_{I}^{alpha} <<
    // zeta_beta
    double aIB = -2; // a_{I}^{beta} <<
    double bIB = 0.50; // b_{I}^{beta} <<
    // association between zeta_alpha and zeta_beta
    double cI = 0.06; // c_I <<
    // * GMRF (alpha and beta)
    double sig1Alpha = 0.005; // sigma_{alpha 0} <<
    double sig1Beta = 0.06; // sigma_{beta 0} <<
    double sigDFix = 0.96; // sigma_{delta} <<
    // * gamma, normal prior
    double mu0Gamma = 0;
    double sig0Gamma = 100;
    // * sig^2_epsilon(v), IG prior
    double a0sigE = 100;
    double b0sigE = 5;
    // * (not in use) sig^2_delta, IG prior
    double a0sigD = 3;
    double b0sigD = 5;
    // * lambda1, normal prior
    double mu0a = 0;
    double sig0a = 100;
    double mu0b = 0;
    double sig0b = 100;

    // change tuning parameters:
    // alpha: 1-sig1Alpha 2-aIA 3-bIA
    // beta:  4-sigDFix 5-sig1Beta 6-aIB 7-bIB
    // AB: 8-cI
    // #choices
    int nTun = 5;
    int nTun1 = nTun; // alpha
    int nTun2 = nTun;
    int nTun3 = nTun;
    int nTun4 = nTun; // beta
    int nTun5 = nTun;
    int nTun6 = nTun;
    int nTun7 = nTun;
    int nTun8 = nTun; // cI
    // tuning choices, specified in .cpp
    VectorXd tun1;
    VectorXd tun2;
    VectorXd tun3;
    VectorXd tun4;
    VectorXd tun5;
    VectorXd tun6;
    VectorXd tun7;
    VectorXd tun8;

    // update tunings in each replication
    int changeTun(int rep);

    Constant();
};


#endif //SIMU1_CONSTANT_H
