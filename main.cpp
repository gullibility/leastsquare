//

#include "stdafx.h"
#include <math.h>
#include <iostream>

using namespace std;

void Matrix_T(double *K, int m, int n, double *KT)//返回矩阵K的转置KT.k[m][n]
{
    int i, j, a, b;
    for (i = 0, a = 0; i<m; i++)
    {
        for (j = 0, b = 0; j<n; j++)
        {
            KT[b + i] = K[a + j];
            b += m;
        }
        a += n;
    }
}

bool Matrix_LU(double *K, int n, double *L, double *U)//对方阵K进行LU分解.分解失败返回False.成功返回True以及分解得到的L与U
{
    int i, j, a, b, c, d;
    double temp;
    for (i = 0, a = 0; i<n; i++)
    {
        for (j = 0; j<n; j++)
        {
            L[a + j] = U[a + j] = 0;
        }
        U[a + i] = 1;
        a += n;
    }
    for (j = 0, d = 0; j<n; j++)
    {
        for (i = j, b = d; i<n; i++)
        {
            temp = 0;
            a = 0, c = j;
            while (a<j)
            {
                temp += L[b + a] * U[c];
                c += n;
                a++;
            }
            L[b + j] = K[b + j] - temp;
            b += n;
        }
        i = j + 1;
        while (i<n)
        {
            temp = 0;
            a = 0, c = i;
            while (a<j)
            {
                temp += L[d + a] * U[c];
                a++;
                c += n;
            }
            if (L[d + j] == 0)
            {
                return false;
            }
            U[d + i] = (K[d + i] - temp) / L[d + j];
            i++;
        }
        d += n;
    }
    return true;
}

bool Matrix_Inv(double *K, int n, double *InvK)//采用LU分解方法求方阵K的逆InvK,K[n][n]
{
    if (1 == n)
    {
        if (K[0] == 0)
        {
            return false;
        }
        else
        {
            InvK[0] = 1 / K[0];
        }
    }
    else if (n<1)
    {
        return false;
    }
    else
    {
        int i, j, a, b;
        double *L, *U, *d, *x, *e, temp;
        a = n*n;
        L = new double[a];
        U = new double[a];
        if (Matrix_LU(K, n, L, U))
        {
            d = new double[n];
            x = new double[n];
            e = new double[n];
            for (i = 0; i<n; i++)
            {
                x[i] = d[i] = 0;
            }
            for (i = 0; i<n; i++)
            {
                for (j = 0; j<n; j++)
                {
                    e[j] = 0;
                }
                e[i] = 1;
                j = 0;
                b = 0;
                while (j<n)
                {
                    temp = 0;
                    a = 0;
                    while (a<j)
                    {
                        temp += d[a] * L[b + a];
                        a++;
                    }
                    d[j] = e[j] - temp;
                    d[j] /= L[b + j];
                    j++;
                    b += n;
                }
                j = n - 1;
                b -= n;
                while (j>-1)
                {
                    temp = 0;
                    a = j + 1;
                    while (a<n)
                    {
                        temp += U[b + a] * x[a];
                        a++;
                    }
                    x[j] = d[j] - temp;
                    x[j] /= U[b + j];
                    j--;
                    b -= n;
                }
                for (j = 0, b = i; j<n; j++)
                {
                    InvK[b] = x[j];
                    b += n;
                }
            }
            delete[]d;
            delete[]x;
            delete[]e;
        }
        else
        {
            delete[]L;
            delete[]U;
            return false;
        }
        delete[]L;
        delete[]U;
    }
    return true;
}

void Matrix_Mul(double *Mul1, int Mul1_m, double *Mul2, int Mul2_n, int nm, double *Mul)
{
    //Mul1[Mul1_m][nm]*Mul2[nm][Mul2_n]=Mul即矩阵的乘法
    int i, j, k, a, b, c, d;
    for (i = 0, a = 0, c = 0; i<Mul1_m; i++)
    {
        for (j = 0; j<Mul2_n; j++)
        {
            b = a + j;
            Mul[b] = 0;
            for (k = 0, d = 0; k<nm; k++)
            {
                Mul[b] += Mul1[c + k] * Mul2[d + j];
                d += Mul2_n;
            }
        }
        c += nm;
        a += Mul2_n;
    }
}

bool Matrix_Solve(double *K, double *B, int m, int n, double *x)//Kx=B求解x。K[m][n]。其结果返回最小二乘解,B[m][1]
{
    double *KT, *Kmul, *Kb, *Kinv;
    int i;
    i = n*n;
    KT = new double[m*n];
    Kmul = new double[i];
    Kinv = new double[i];
    Kb = new double[n];
    Matrix_T(K, m, n, KT);
    Matrix_Mul(KT, n, K, n, m, Kmul);
    Matrix_Mul(KT, n, B, 1, m, Kb);
    if (Matrix_Inv(Kmul, n, Kinv))
    {
        Matrix_Mul(Kinv, n, Kb, 1, n, x);
        delete[]KT;
        delete[]Kmul;
        delete[]Kinv;
        delete[]Kb;
        return true;
    }
    else
    {
        delete[]KT;
        delete[]Kmul;
        delete[]Kinv;
        delete[]Kb;
        return false;
    }
}

bool Matrix_PolyFit(double *x, double *y, int X_Y_Number, int Fit_N, double *ks)
{
    /*
    对X_Y_Number组数据x,y进行Fit_N次多项式拟合,拟合返回多项式的系数ks
    x[X_Y_Number],y[X_Y_Number],ks[Fit_N+1]
    */
    if (Fit_N>X_Y_Number || X_Y_Number<1)
    {
        return false;
    }
    int i, j, index, n;
    double temp, *x2, *y2;
    Fit_N++;
    y2 = new double[Fit_N];
    x2 = new double[Fit_N*Fit_N];
    for (i = 0, index = 0; i<Fit_N; i++)
    {
        y2[i] = 0;
        for (j = 0; j<Fit_N; j++)
        {
            x2[index + j] = 0;
        }
        index += Fit_N;
    }
    x2[0] = X_Y_Number;
    for (i = 0; i<Fit_N; i++)
    {
        for (j = i + 1; j<Fit_N; j++)
        {
            temp = 0;
            n = i + j;
            for (index = 0; index<X_Y_Number; index++)
            {
                temp += pow(x[index], n);
            }
            index = j;
            for (n = i; n<Fit_N; n++)
            {
                if (index >= 0)
                {
                    x2[n*Fit_N + index] = temp;
                }
                index--;
            }
        }
    }
    n = Fit_N + Fit_N - 2;
    temp = 0;
    for (i = 0; i<X_Y_Number; i++)
    {
        temp += pow(x[i], n);
    }
    x2[Fit_N*Fit_N - 1] = temp;
    for (i = 0; i<Fit_N; i++)
    {
        temp = 0;
        for (j = 0; j<X_Y_Number; j++)
        {
            temp += y[j] * pow(x[j], i);
        }
        y2[i] = temp;
    }

    if (Matrix_Solve(x2, y2, Fit_N, Fit_N, ks))
    {
        delete[]y2;
        delete[]x2;
        return true;
    }
    else
    {
        delete[]y2;
        delete[]x2;
        return false;
    }

}

int main()
{
    double x[] = { 1,2,3,4,5,6,7,8,9,10 };
    double y[] = { 9,31,83,177,325,539,831,1213,1697,2295 };
    double ks[4];
    Matrix_PolyFit(x, y, 10, 3, ks);
    cout << ks[0] << "," << ks[1] << "," << ks[2] << "," << ks[3] << endl;

    double x2[] = { 1,2,3,4,5,6 };
    double y2[] = { 4,12,32,70,132,224};
    double ks2[4];
    Matrix_PolyFit(x2, y2, 6, 3, ks2);
    cout << ks2[0] << "," << ks2[1] << "," << ks2[2] << "," << ks[3] << endl;
    system("pause");
    return 0;
}
