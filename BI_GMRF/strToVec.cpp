//
//

// change a string to vector<double>
#include "strToVec.h"
using namespace std;

void strToVec(string str, vector<double>& fea)
{
    stringstream ss(str);
    string buf;
    vector<double> vec;
    //vector<double> vec(1024);
    while (ss >> buf)
        vec.push_back(atof(buf.c_str()));

    fea = vec;
}



