//
//

/*
 how to use this class:
 steps:
 0. Mcmc:mcmc(Constant constant)
 1. (not in use) readData(Constant constant, Generate generate)
 2. setInitial(Constant constant)
 3. fastComp(Constant constant)
 3. update
 */

#include "Mcmc.h"

using namespace std;
using namespace Eigen;

Mcmc::Mcmc(Constant constant)
        : x(constant.n),
          z(constant.n),
          m(constant.vv, constant.n),
          mb(constant.vv, constant.n),
          y(constant.n),
          alpha(constant.vv),
          a(constant.vv),
          zetaA(constant.vv),
          beta(constant.vv),
          zetaB(constant.vv),
          sigE(constant.vv),
          nbh(constant.vv, constant.mNnbh),
          nNbh(constant.vv),
          sumXM(constant.vv),
          sumXMb(constant.vv),
          sumM2(constant.vv),
          sumMb2(constant.vv),
          sumYM(constant.vv),
          sumYMb(constant.vv),
          sumBetaMb(constant.n),
          sumZM(constant.vv),
          sumZMb(constant.vv)
{
    x = VectorXd::Constant(constant.n, 0);
    z = VectorXd::Constant(constant.n, 0);
    m = MatrixXd::Constant(constant.vv, constant.n, 0);
    mb = MatrixXd::Constant(constant.vv, constant.n, 0);
    y = VectorXd::Constant(constant.n, 0);
    alpha = VectorXd::Constant(constant.vv, 0);
    a = VectorXd::Constant(constant.vv, 0);
    zetaA = VectorXd::Constant(constant.vv, 0);
    beta = VectorXd::Constant(constant.vv, 0);
    zetaB = VectorXd::Constant(constant.vv, 0);
    sigE = VectorXd::Constant(constant.vv, 0);
    nbh = MatrixXd::Constant(constant.vv, constant.mNnbh, 0);
    nNbh = VectorXd::Constant(constant.vv, 0);
    sumXM = VectorXd::Constant(constant.vv, 0);
    sumXMb = VectorXd::Constant(constant.vv, 0);
    sumM2 = VectorXd::Constant(constant.vv, 0);
    sumMb2 = VectorXd::Constant(constant.vv, 0);
    sumYM = VectorXd::Constant(constant.vv, 0);
    sumYMb = VectorXd::Constant(constant.vv, 0);
    sumBetaMb = VectorXd::Constant(constant.n, 0);
    sumZM = VectorXd::Constant(constant.vv, 0);
    sumZMb = VectorXd::Constant(constant.vv, 0);
}

// read data (x, m, y) and nbh, nNbh from .txt
int Mcmc::readData(Constant constant, Generate generate, int rep)
{
    // *** (not in use) read from "Generate" class
    // x = generate.getX();
    // z = generate.getZ();
    // m = generate.getM();
    // mb = m;
    // y = generate.getY();
    // nbh = generate.getNbh();
    // nNbh = generate.getNnbh();


    // *** read from file
    if (constant.setting == 0)
    {
        // * x
        // two tmp variables
        string strX; // store a line as str
        vector<double> vecX; // store a line as vec
        // open file
        ifstream fileX;
        fileX.open("../data/x_gen.txt");
        // assign
        int i = 0;
        while (getline(fileX, strX))
        {
            strToVec(strX, vecX);
            x(i) = vecX[0];
            i++;
        } // end while i
        fileX.close();

        // * z
        // two tmp variables
        string strZ; // store a line as str
        vector<double> vecZ; // store a line as vec
        // open file
        ifstream fileZ;
        fileZ.open("../data/c1_gen.txt");
        // assign
        i = 0;
        while (getline(fileZ, strZ))
        {
            strToVec(strZ, vecZ);
            z(i) = vecZ[0];
            i++;
        } // end while i
        fileZ.close();

        // * m/mb
        // two tmp variables
        string strM; // store a line as str
        vector<double> vecM; // store a line as vec
        // open file
        ifstream fileM;
        fileM.open("../data/m_gen.txt");
        // assign
        i = 0;
        while (getline(fileM, strM))
        {
            // convert str to vec
            strToVec(strM, vecM);
            // assign row_i in file to col_i in matrix xi
            for (int j=0; j<constant.vv; j++)
            {
                m(j,i) = vecM[j];
            } // end j
            i++;
        } // end while i
        fileM.close();
        // currently mb = m, No Cv
        mb = m;

        // * y
        // two tmp variables
        string strY; // store a line as str
        vector<double> vecY; // store a line as vec
        // open file
        ifstream fileY;
        fileY.open("../data/y_gen.txt");
        // assign
        i = 0;
        while (getline(fileY, strY))
        {
            strToVec(strY, vecY);
            y(i) = vecY[0];
            i++;
        } // end while i
        fileY.close();

        // * nbh
        // two tmp variables
        string strN; // store a line as str
        vector<double> vecN; // store a line as vec
        // open file
        ifstream fileN;
        fileN.open("../data/nbh.txt");
        // assign
        i = 0;
        while (getline(fileN, strN))
        {
            // convert str to vec
            strToVec(strN, vecN);
            // assign row_i in file to col_i in matrix xi
            for (int j=0; j<constant.mNnbh; j++)
            {
                nbh(i,j) = vecN[j];
            } // end j
            i++;
        } // end while i
        fileN.close();

        // * nNbh
        // two tmp variables
        string strNn; // store a line as str
        vector<double> vecNn; // store a line as vec
        // open file
        ifstream fileNn;
        fileNn.open("../data/nNbh.txt");
        // assign
        i = 0;
        while (getline(fileNn, strNn))
        {
            strToVec(strNn, vecNn);
            nNbh(i) = vecNn[0];
            i++;
        } // end while i
        fileNn.close();
    }

    if (constant.setting == 1)
    {
        // * x
        // two tmp variables
        string strX; // store a line as str
        vector<double> vecX; // store a line as vec
        // open file
        ifstream fileX;
        fileX.open("../data/x_gen.txt");
        // assign
        int i = 0;
        while (getline(fileX, strX))
        {
            //cout << "line=" << i << endl;
            strToVec(strX, vecX);
            if ((i >= constant.n*rep) & (i<constant.n*(rep+1)))
            {
                x(i-constant.n*rep) = vecX[0];
            }
            i++;
        } // end while i
        fileX.close();
        //cout << "xRep[1:5]=" << x.head(5) << endl;

        // * z
        // two tmp variables
        string strZ; // store a line as str
        vector<double> vecZ; // store a line as vec
        // open file
        ifstream fileZ;
        fileZ.open("../data/c1_gen.txt");
        // assign
        i = 0;
        while (getline(fileZ, strZ))
        {
            strToVec(strZ, vecZ);
            if ((i >= constant.n*rep) & (i<constant.n*(rep+1)))
            {
                z(i-constant.n*rep) = vecZ[0];
            }
            i++;
        } // end while i
        fileZ.close();
        //cout << "zRep[1:5]=" << z.head(5) << endl;

        // * m/mb
        // two tmp variables
        string strM; // store a line as str
        vector<double> vecM; // store a line as vec
        // open file
        ifstream fileM;
        fileM.open("../data/m_gen.txt");
        // assign
        i = 0;
        while (getline(fileM, strM))
        {
            // convert str to vec
            strToVec(strM, vecM);
            // assign row_i in file to col_i in matrix xi
            if ((i >= constant.n*rep) & (i<constant.n*(rep+1)))
            {
                for (int j=0; j<constant.vv; j++)
                {
                    m(j,i-constant.n*rep) = vecM[j];
                } // end j
            } // end if i
            i++;
        } // end while i
        fileM.close();
        // currently mb = m, No Cv
        mb = m;
        if (constant.applyCv == 1)
        {
            mb = m * constant.cv;
            cout << "Cv is applied to m" << endl;
        }
        //cout << "mbRep[1:5,1:5]=" << mb.block(0,0,5,5) << endl;

        // * y
        // two tmp variables
        string strY; // store a line as str
        vector<double> vecY; // store a line as vec
        // open file
        ifstream fileY;
        fileY.open("../data/y_gen.txt");
        // assign
        i = 0;
        while (getline(fileY, strY))
        {
            strToVec(strY, vecY);
            if ((i >= constant.n*rep) & (i<constant.n*(rep+1)))
            {
                y(i-constant.n*rep) = vecY[0];
            }
            i++;
        } // end while i
        fileY.close();
        //cout << "yRep[1:5]=" << y.head(5) << endl;

        // * nbh
        // two tmp variables
        string strN; // store a line as str
        vector<double> vecN; // store a line as vec
        // open file
        ifstream fileN;
        fileN.open("../data/nbh.txt");
        // assign
        i = 0;
        while (getline(fileN, strN))
        {
            // convert str to vec
            strToVec(strN, vecN);
            // assign row_i in file to col_i in matrix xi
            for (int j=0; j<constant.mNnbh; j++)
            {
                nbh(i,j) = vecN[j];
            } // end j
            i++;
        } // end while i
        fileN.close();

        // * nNbh
        // two tmp variables
        string strNn; // store a line as str
        vector<double> vecNn; // store a line as vec
        // open file
        ifstream fileNn;
        fileNn.open("../data/nNbh.txt");
        // assign
        i = 0;
        while (getline(fileNn, strNn))
        {
            strToVec(strNn, vecNn);
            nNbh(i) = vecNn[0];
            i++;
        } // end while i
        fileNn.close();
    }

    return 0;
}

// set initial value for mcmc updates
int Mcmc::setInitial(Constant constant, Initial initial)
{
    alpha = initial.alpha;
    a = initial.a;
    zetaA = initial.zetaA;
    beta = initial.beta;
    zetaB = initial.zetaB;
    gamma = initial.gamma;
    b = initial.b;
    sigE = initial.sigE;
    sigD = initial.sigD;

    //cout << "MCMC_ini_beta2147=" << beta(2146) << endl;

    return 0;
}

// for fast computation
int Mcmc::fastComp(Constant constant)
{
    // fixed
    sumX2 = pow(x.array(),2).sum(); // sum_{i=1}^n x_i^2, scalar
    sumXM = ((m.array().rowwise()) * (x.transpose().array())).rowwise().sum(); // sum_{i=1}^n m_i(v) x_i, vector, length = V
    sumXMb = ((mb.array().rowwise()) * (x.transpose().array())).rowwise().sum();
    sumM2 = pow(m.array(), 2).rowwise().sum(); // sum_{i=1}^n m_i(v)^2, vector, length = V
    sumMb2 = pow(mb.array(), 2).rowwise().sum();
    sumYM = ((m.array().rowwise()) * (y.transpose().array())).rowwise().sum();  // sum_{i=1}^n y_i m_i(v), vector, length = V
    sumYMb = ((mb.array().rowwise()) * (y.transpose().array())).rowwise().sum();
    // need to be updated along with beta
    sumBetaMb = beta.transpose() * mb; // sum_{v=1}^V beta(v) m_i(v), vector, length = n
    // related to a, b
    sumZ2 = pow(z.array(),2).sum(); // sum_{i=1}^n z_i^2, scalar
    sumZM = ((m.array().rowwise()) * (z.transpose().array())).rowwise().sum(); // sum_{i=1}^n m_i(v) z_i, vector, length = V
    sumZMb = ((mb.array().rowwise()) * (z.transpose().array())).rowwise().sum();
    sumZX = (x.array() * z.array()).sum(); // sum_{i=1}^n x_i z_i, scalar

    return 0;
}

// get values
VectorXd Mcmc::getAlpha(void)
{
    return alpha;
}

VectorXd Mcmc::getA(void)
{
    return a;
}

VectorXd Mcmc::getZetaA(void)
{
    return zetaA;
}

VectorXd Mcmc::getBeta(void)
{
    //cout << "MCMC_getBeta_beta2147=" << beta(2146) << endl;
    return beta;
}

VectorXd Mcmc::getZetaB(void)
{
    return zetaB;
}

double Mcmc::getGamma(void)
{
    return gamma;
}

double Mcmc::getB(void)
{
    return b;
}

VectorXd Mcmc::getSigE(void)
{
    return sigE;
}

double Mcmc::getSigD(void)
{
    return sigD;
}

VectorXd Mcmc::getX(void)
{
    return x;
}

VectorXd Mcmc::getZ(void)
{
    return z;
}

MatrixXd Mcmc::getM(void)
{
    return m;
}

VectorXd Mcmc::getY(void)
{
    return y;
}

// update alpha, zeta_alpha, a (lambda1), and sigma^2_epsilon
int Mcmc::updateA(Constant constant)
{
    for (int v=0; v<constant.vv; v++)
    {
        // *** update alpha and zeta_alpha
        // generate alpha(v) proposal from alpha(v)|zeta^alpha = 1, all others
        double sig2 = 1/1/(sumX2/sigE(v) + 1/constant.sig1Alpha);
        double barAlphaV = 0;
        double zb01 = 0; // sum((zetaAlpha[nbh[v,]] == 0)-(zetaAlpha[nbh[v,]] == 1), na.rm = T)
        for (int i=0; i<constant.mNnbh; i++)
        {
            if (nbh(v,i) >= 0)
            {
                barAlphaV = barAlphaV + alpha(nbh(v,i));
                zb01 = zb01 + (zetaA(nbh(v,i)) == 0) - (zetaA(nbh(v,i)) == 1);
            }
        }
        barAlphaV = barAlphaV/nNbh(v);
        double mu = sumXM(v)/sigE(v) - a(v)*sumZX/sigE(v) + barAlphaV/constant.sig1Alpha;
        mu = sig2 * mu;
        double alphaVP = mu + rnorm(1)(0)*sqrt(sig2);
        // compute g_v
        // where g_v = P(alpha(v)=0, zeta(v)=0|all) / P(alpha(v)=alpha(v)P, zeta(v)=1|all)
        double gv = 0;
        gv = gv -0.5/sigE(v)*(2*alphaVP*sumXM(v)-a(v)*sumZX-pow(alphaVP,2)*sumX2) + 0.5/constant.sig1Alpha*pow(alphaVP - barAlphaV,2);
        //gv = gv - constant.aIA + 2*constant.bIA*zb01 + constant.cI*((zetaB(v) == 0) - (zetaB(v) == 1));
        gv = gv - constant.aIA + 2*constant.bIA*zb01 - constant.cI*zetaB(v);
        gv = sqrt(2*constant.pi*constant.sig1Alpha)*exp(gv);
        // prob of alpha(v)=alphaVP, zeta(v)=1| all = 1/(g_v + 1)
        double prob = 1/(gv+1);
        // ac
        double rand = runif();
        if (rand < prob)
        {
            alpha(v) = alphaVP;
            zetaA(v) = 1;
        }
        else
        {
            alpha(v) = 0;
            zetaA(v) = 0;
        }

        //alpha(v) = 0; // fix alpha=0 for debugging
        //zetaA(v) = 0;

        // update a(v)
        double sigStar = 1/(sumZ2/sigE(v) + 1/constant.sig0a);
        double muStar = (sumZM(v) - alpha(v)*sumZX)/sigE(v) + constant.mu0a/constant.sig0a;
        muStar = sigStar * muStar;
        a(v) = muStar + rnorm(1)(0)*sqrt(sigStar);

        // a(v) = 0; // fix a = 0 for debugging

        // *** update sigE(v)
        double aS = constant.a0sigE + 0.5*constant.n;
        double bS = constant.b0sigE + 0.5*pow(m.row(v).array()-alpha(v)*x.transpose().array()-a(v)*z.transpose().array(), 2).sum();
        sigE(v) = 1/rgamma(1, aS, 1/bS)(0); // NOTE rgamma parameterization
    }

    return 0;
}

// update beta and zeta_beta
int Mcmc::updateB(Constant constant)
{
    for (int v=0; v<constant.vv; v++)
    {
        // sum_{v1!=v} beta(v1) m_i(v1)
        if (beta(v) != 0)
        {
            sumBetaMb = sumBetaMb - beta(v) * mb.row(v).transpose();
        }
        // generate beta(v) proposal from beta(v)|zeta^beta = 1, all others
        double sig2 = 1/(sumMb2(v)/sigD + 1/constant.sig1Beta);
        double barBetaV = 0;
        double zb01 = 0; // sum((zetaBeta[nbh[v,]] == 0)-(zetaBeta[nbh[v,]] == 1), na.rm = T)
        for (int i=0; i<constant.mNnbh; i++)
        {
            if (nbh(v,i) >= 0)
            {
                barBetaV = barBetaV + beta(nbh(v,i));
                zb01 = zb01 + (zetaB(nbh(v,i)) == 0) - (zetaB(nbh(v,i)) == 1);
            }
        }
        barBetaV = barBetaV/nNbh(v);
        double mu = 0;
        mu = (sumYMb(v) - gamma*sumXMb(v) - b*sumZMb(v))/sigD - (mb.row(v).transpose().array()*sumBetaMb.array()).sum() - barBetaV/constant.sig1Beta;
        mu = sig2 * mu;
        double betaVP = mu + rnorm(1)(0)*sqrt(sig2);
        // compute g_v
        // where g_v = P(beta(v)=0, zeta(v)=0|all) / P(beta(v)=beta(v)P, zeta(v)=1|all)
        double gv = 0;
        gv = -0.5/sigD * (2*betaVP*(sumYMb(v) - (mb.row(v).transpose().array()*sumBetaMb.array()).sum() - gamma*sumXMb(v) - b*sumZMb(v)) - pow(betaVP,2)*sumMb2(v));
        gv = gv + 0.5/constant.sig1Beta * pow(betaVP - barBetaV,2);
        gv = gv - constant.aIB + 2*constant.bIB*zb01 - constant.cI*zetaA(v);
        gv = sqrt(2*constant.pi*constant.sig1Beta) * exp(gv);
        // prob of beta(v)=betaVP, zeta(v)=1| all = 1/(g_v + 1)
        double prob = 1/(gv+1);
        // ac
        double rand = runif();
        if (rand < prob)
        {
            beta(v) = betaVP;
            zetaB(v) = 1;
            // update sumBetaM = sum_{v=1}^V beta(v) m_i(v)
            sumBetaMb = sumBetaMb + beta(v) * mb.row(v).transpose();
        }
        else
        {
            beta(v) = 0;
            zetaB(v) = 0;
        }
    }

    return 0;
}

// update gamma
int Mcmc::updateGamma(Constant constant)
{
    double sigStar = 1/(sumX2/sigD + 1/constant.sig0Gamma);
    double muStar = ((y - sumBetaMb - b*z).array()*x.array()).sum()/sigD + constant.mu0Gamma/constant.sig0Gamma;
    muStar = sigStar * muStar;
    gamma = muStar + rnorm(1)(0)*sqrt(sigStar);

    //cout << "MCMC_beta2147=" << beta(2146) << endl;

    return 0;
}

// update b
int Mcmc::updateb(Constant constant)
{
    double sigStar = 1/(sumZ2/sigD + 1/constant.sig0b);
    double muStar = ((y - sumBetaMb - gamma*x).array()*z.array()).sum()/sigD + constant.mu0b/constant.sig0b;
    muStar = sigStar * muStar;
    b = muStar + rnorm(1)(0)*sqrt(sigStar);
    //cout << "muStar=" << muStar << " sigStar=" << sigStar << " b=" << b << endl;
    //cout << "sumZ2=" << sumZ2 << endl;

    return 0;
}

// update sigD
int Mcmc::updateSigD(Constant constant)
{
    double aS = constant.a0sigD + 0.5*constant.n;
    double bS = constant.b0sigD + pow((y - sumBetaMb - gamma*x - b*z).array(),2).sum();
    sigD = 1/rgamma(1, aS, 1/bS)(0); // NOTE rgamma parameterization

    return 0;
}

