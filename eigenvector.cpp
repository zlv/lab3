#include <iostream>
#include <string>
#include <sstream>
#include <string.h>
#include <cmath>
#include "PolStr.h"
#include <stdexcept>
using namespace std;
double** findP(int,double**);
double* findLa(int,double**,const double,double**);
void mul(int,double**,double**,double**);
void mulV(int, double*, double**, double*);
int main(int argc, char **argv) {
        int type;
        int n;
        double **matrix;
        double **b;
        cin >>type >> n;
        matrix=new double*[n];
        double eps=1e-5;
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
            double **P;
            if(type==1) {
               P=findP(n,matrix);
               findLa(n,P,eps,matrix);
               
            }
            for (int i=0; i<n; i++) {
                delete[] matrix[i];
                delete[] P[i];
            }
            delete[] matrix;
            delete[] P;
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

double* findLa(int n, double** a, const double eps,double** A) 
{
    std::stringstream ss;
    string str;
    bool h=n%2;
    ss<<(h?" ":" - ")<<"x^"<<n;
    for(int i=0;i<n;i++)
        ss<< (a[0][i]<0&&h?" + ":" - ")<<fabs(a[0][i])<<" * x^"<<n-i-1;
    char* t = new char[1024+1];
    strcpy(t,ss.str().c_str());
    cout<<ss.str() << endl << str<<'\n';
    t=CreatePolStr(t,0);
    cout<<t<<'\n';
      
    const double delta = eps*1e2;
    const double minx = -1;
    const double maxx = 4;
    double *xi = new double[n];
    double *ximult = new double[n];
    int lastxi = -1;
    int multsum = 0;
    for (double x=minx; x<maxx; x+=eps) {
        if (fabs(EvalPolStr(t,x))<delta) {
            if (lastxi==-1 || fabs(xi[lastxi]-x)>2*eps) {
                //cout << x << endl;
                xi[++lastxi] = x;
                ximult[multsum++] = 1;
            }
            else if (lastxi!=-1) {
                //cout << "same: " << x << endl;
                xi[lastxi] = x;
            }
        }
    }
    int i = 0;
    while (multsum<n) {
        int derIndex = 1;
        while (fabs(EvalPolStr(t,xi[i],derIndex++))<delta) {
            ximult[i++]++;
            multsum++;
            if (multsum>=n) break;
        }
    }
    double **y = new double*[n];

    for (int i=0; i<=lastxi; i++)
    {
        cout <<">>"<<endl<< xi[i] << ' ' << ximult[i] << endl<<"============"<<endl;
        for(int j=0;j<n;j++)
        {
                y[i]=new double[n];
                y[i][j]=pow(xi[i],j);
                cout<<y[i][j]<<endl;
        }
    }
    cout <<"========"<< endl;
    double* r=new double[n];
    mulV(n,r,A,y[0]);
    for(int i=0;i<n;i++)
        cout<<r[i]<<'\n';
    delete[] xi;
    delete[] ximult;
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
void mulV(int n, double* res, double** a, double* b) {
    for(int i=0;i<n;i++) {
            res[i] = 0;
            for(int k=0;k<n;k++) {
                res[i]+=a[i][k]*b[k];
            }
    }
}
