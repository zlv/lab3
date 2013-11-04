#include <iostream>
#include <string>
#include <sstream>
#include <string.h>
#include <cmath>
#include <stdexcept>
using namespace std;
double** findP(int,double**,double ***s);
double* findLa(int,double**,const double,double**,double**);
void mul(int,double**,double**,double**);
void mulV(int, double*, double**, double*);
void mulC(int, double*, double*, double);
void mulC(int, double**, double**, double);
double determ(double**, int);
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
            double **s;
            if(type==1) {
               P=findP(n,matrix,&s);
               findLa(n,P,eps,matrix,s);
               
            }
            for (int i=0; i<n; i++) {
                delete[] matrix[i];
                delete[] s[i];
                delete[] P[i];
            }
            delete[] matrix;
            delete[] P;
            delete[] s;
        }
        catch (invalid_argument &e) {
            cerr << e.what() << endl;
        }
                
        return 0;
}
double** findP(int n, double** a, double ***s) {
    double **acur = new double*[n];
    double **m = new double*[n];
    double **minv = new double*[n];
    double **e = new double*[n];
    double **temp = new double*[n];
    *s = new double*[n];
    for (int i=0; i<n; i++) {
        acur[i] = new double[n];
        m[i] = new double[n];
        minv[i] = new double[n];
        e[i] = new double[n];
        temp[i] = new double[n];
        (*s)[i] = new double[n];
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
        if (k==0) {
            for(int i=0;i<n;i++) {
                for(int j=0;j<n;j++) {
                    (*s)[i][j] = m[i][j];
                }
            }
        }
        else {
            mul(n,temp,*s,m);
            for(int i=0;i<n;i++) {
                for(int j=0;j<n;j++) {
                    (*s)[i][j] = temp[i][j];
                }
            }
        }
    }
    cout << "Frobenius matrix" << endl;
    for(int i=0;i<n;i++) {
        for(int j=0;j<n;j++) {
            cout << acur[i][j] << ' ';
        }
        cout << endl;
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

double findDerMult(int n, int iDerivative) {
    double advancedDerMult = 1;
    for (int i=0; i<iDerivative; i++)
        advancedDerMult *= n-i;
    return advancedDerMult;
}

double findSumElement(double x, double a, int powIndex, bool h, double advancedDerMult) {
    return a*pow(x,powIndex)*(h?1:-1)*advancedDerMult;
}

double eval(double* a, int n, double x, int iDerivative=0) {
    bool h=n%2;
    double sum = findSumElement(x,-1,n-iDerivative,h,findDerMult(n,iDerivative));
    for (int i=0; i<n; i++) {
        if (n-i-1-iDerivative<0)
            break;
        sum += findSumElement(x,a[i],n-i-1-iDerivative,h,findDerMult(n-i-1,iDerivative));
    }
    return sum;
}

double* findLa(int n, double** a, const double eps,double** A,double** S) 
{
    std::stringstream ss;
    string str;
    bool h=n%2;
    ss<<(!h?"":" - ")<<"x^"<<n;
    for(int i=0;i<n;i++)
        ss<< (a[0][i]<0^h?" + ":" - ")<<fabs(a[0][i])<<" * x^"<<n-i-1;
    cout<<ss.str() << endl;
      
    const double delta = eps*1e2;
    const double minx = -3;
    const double maxx = 4;
    double *xi = new double[n];
    double *ximult = new double[n];
    int lastxi = -1;
    int multsum = 0;
    double dxil = minx;
    for (int i=-1; i<5; i++)
        cerr << "x = " << i << " y = " << eval(a[0],n,i) << endl;
    bool foundValue = 0;
    for (double x=minx; x<maxx; x+=eps) {
        if (eval(a[0],n,x)<delta) {
            if (lastxi==-1 || fabs(dxil-x)>2*eps) {
                xi[++lastxi] = x;
                ximult[multsum++] = 1;
                foundValue = 1;
            }
            dxil = x;
        }
        else if (foundValue && fabs(xi[lastxi]-dxil)>2*eps) {
            xi[lastxi] = (xi[lastxi]+dxil)/2;
            foundValue = 0;
        }
    }
    int i = 0;
    while (multsum<n) {
        int derIndex = 1;
        while (fabs(eval(a[0],n,xi[i],derIndex++))<delta) {
            cerr << xi[i] << endl;
            ximult[i++]++;
            multsum++;
            if (multsum>=n) break;
        }
    }
    double **y = new double*[n];

    cout << "eigennumber : ";
    for (int i=0; i<=lastxi; i++) {
        cout << xi[i] << ' ';
    }
    cout << "\neigen multiplicity : ";
    for (int i=0; i<=lastxi; i++) {
        cout << ximult[i] << ' ';
    }
    cout << endl;
    for (int i=0; i<=lastxi; i++)
    {
        y[i]=new double[n];
        for(int j=0;j<n;j++)
        {
                y[i][j]=pow(xi[i],n-1-j);
        }
    }
    double** tt=new double*[n];
    double** r=new double*[n];
    for(int i=0;i<n;i++)
    {
            tt[i]=new double[n];
            r[i]=new double[n];
            mulV(n,tt[i],S,y[i]);
            mulV(n,r[i],A,tt[i]);
            mulC(n,tt[i],r[i],xi[i]);
    }
    
    double **e = new double*[n];
    for (int i=0; i<n; i++) {
        e[i] = new double[n];
        for (int j=0; j<n; j++) {
            e[i][j] = i==j;
        }
    }
    for(int i=0;i<n;i++)
    {
    
        double **eTmp = new double*[n];
        for(int j=0;j<n;j++)
                eTmp[j]=new double[n];
        mulC(n,eTmp,e,xi[i]);
        cout <<"|A-eig" << i << "*E|: ";
        for(int k=0;k<n;k++) {
                for(int l=0;l<n;l++) {
                    
                    eTmp[k][l]=A[k][l]-eTmp[k][l];
                }
        }
        
        cout<<determ(eTmp,n) << endl;

        
    }
    for(int i=0;i<n;i++)
    {
        cout <<"A*x" << i << "-eig" << i << "*x" << i << ":\n";
        for(int j=0;j<n;j++)
                cout<<r[i][j]-tt[i][j]<<'\n';
    }
    for(int i=0;i<n;i++)
        delete[] y[i];
    delete[] r;
    delete[] y;
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
void mulC(int n, double* res, double* a, double b) {
    for(int i=0;i<n;i++) {
            res[i] = a[i]*b;
    }
}
void mulC(int n, double** res, double** a, double b) {
    for(int i=0;i<n;i++) {
        for(int j=0;j<n;j++)
            res[i][j] = a[i][j]*b;
    }
}
 
double determ(double** A, int n)
{
        int i,j;
        double det=0;
        double** matr;
        if(n==1)
        {
                det=A[0][0];
        }
        else if(n==2)
        {
                det=A[0][0]*A[1][1]-A[0][1]*A[1][0];
        }
        else
        {
                matr=new double*[n-1];
                for(i=0;i<n;++i)
                {
                        for(j=0;j<n-1;++j)
                        {
                                if(j<i) 
                                        matr[j]=A[j];
                                else
                                        matr[j]=A[j+1];
                        }
                        det+=pow((double)-1, (i+j))*determ(matr, n-1)*A[i][n-1];
                }
                delete[] matr;
        }
        return det;
}
