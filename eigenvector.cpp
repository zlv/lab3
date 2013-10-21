#include <iostream>
#include <string>
#include <sstream>
#include <string.h>
#include <cmath>
#include "PolStr.h"
#include <stdexcept>
using namespace std;
double** findP(int,double**);
double* findLa(int,double**);
void mul(int,double**,double**,double**);
int main(int argc, char **argv) {
        int type;
        int n;
        double **matrix;
        double **b;
        cin >>type >> n;
        matrix=new double*[n];
        double eps=1e-3;
        cout.precision(-log10(eps));
        cout << "matrix input:\n";
        for(int i=0;i<n;i++)
        {
                matrix[i]= new double[n];
                for(int j=0;j<n;j++) {
                        cin >> matrix[i][j];
                        cout << matrix[i][j] << ' ';
                }
                cout << endl;
        }
        double epsf=0;
        double *epsv = new double[n];
        try {
            double **result;
            if(type==1) {
               double** P=findP(n,matrix);
               findLa(n,P);
            }
            for (int i=0; i<n; i++) {
                delete[] matrix[i];
            }
            delete[] matrix;
        }
        catch (invalid_argument e) {
            cerr << e.what() << endl;
        }
                
        return 0;
}
double** findP(int n, double** a) {
    double **acur = new double*[n];
    double **m = new double*[n];
    double **minv = new double*[n];
    double **e = new double*[n];
    double **temp = new double*[n];
    for (int i=0; i<n; i++) {
        acur[i] = new double[n];
        m[i] = new double[n];
        minv[i] = new double[n];
        e[i] = new double[n];
        temp[i] = new double[n];
        for (int j=0; j<n; j++) {
            acur[i][j] = a[i][j];
            e[i][j] = i==j;
        }
    }
    for (int k=0; k<n-1; k++) {
        int realk = n-k-2;
        for (int i=0; i<n; i++) {
            for (int j=0; j<n; j++) {
                if (i!=realk) {
                    m[i][j] = minv[i][j] = e[i][j];
                }
                else {
                    minv[i][j] = acur[realk+1][j];
                    if (j!=realk) {
                        m[i][j] = -acur[realk+1][j]/acur[realk+1][realk];
                    }
                    else {
                        m[i][j] = 1/acur[realk+1][realk];
                    }
                }
            }
        }
        mul(n,temp,minv,acur);
        mul(n,acur,temp,m);
        cout << "current matrix" << k << endl;
        for(int i=0;i<n;i++) {
            for(int j=0;j<n;j++) {
                cout << acur[i][j] << ' ';
            }
            cout << endl;
        }
    }
    for (int i=0; i<n; i++) {
        delete[] m[i];
        delete[] minv[i];
        delete[] e[i];
        delete[] temp[i];
    }
    delete[] m;
    delete[] minv;
    delete[] e;
    delete[] temp;
    return acur;
}

double* findLa(int n, double** a) {
        std::stringstream ss;
        string str;
        int h=-1*n%2;
        ss<<"x^"<<n<<" + ";
        for(int i=0;i<n;i++)
                ss<<a[0][i]<<" * x^"<<n-i-1<<" + ";
        ss<<"0";
//        ss>>str;

        char* t = new char[1024+1];
              str="x^4 + 3 * x^3 + 3.9 * x^2 + 18.66 * x^1 + 14.08 * x^0 + 0";
        strcpy(t,str.c_str());
       cout<<str<<'\n';
       
        t=CreatePolStr(t,0);
         cout<<t<<'\n';
      
        double aa=1;
        cout<<EvalPolStr(t,aa)<<'\n';
  
//        cout<<str;
}
void mul(int n, double** res, double** a, double** b) {
    for(int i=0;i<n;i++) {
        for(int j=0;j<n;j++) {
            res[i][j] = 0;
            for(int k=0;k<n;k++) {
                res[i][j]+=a[i][k]*b[k][j];
            }
        }
    }
}
