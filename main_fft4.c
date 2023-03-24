#include <stdio.h>
#include <math.h>
#include <time.h>

#define PI 3.1415926535


void kfft(double pr[], double pi[], int n, int k, double mag[])
{
    int it, m, is, i, j, nv, l0;
    double p, q, s, vr, vi, poddr, poddi;

    double fr[1024] = {0.0};
    double fi[1024] = {0.0};

    for (it = 0; it < n; it++) //将pr[0]和pi[0]循环赋值给fr[]和fi[]
    {
        m = it;
        is = 0;
        for (i = 0; i < k; i++)
        {
            j = m / 2;
            is = 2 * is + (m - 2 * j);
            m = j;
        }
        fr[it] = pr[is];
        //fi[it] = pi[is];
    }

    pr[0] = 1.0;
    pi[0] = 0.0;
    p = 6.283185306 / n;
    pr[1] = cos(p); //将w=e^-j2pi/n用欧拉公式表示
    pi[1] = -sin(p);

    for (i = 2; i < n; i++) //计算pr[]
    {
        p = pr[i - 1] * pr[1];
        q = pi[i - 1] * pi[1];
        s = (pr[i - 1] + pi[i - 1]) * (pr[1] + pi[1]);
        pr[i] = p - q;
        pi[i] = s - p - q;
    }

    for (it = 0; it <= n - 2; it = it + 2)
    {
        vr = fr[it];
        //vi = fi[it];
        fr[it] = vr + fr[it + 1];
        //fi[it] = vi + fi[it + 1];
        fr[it + 1] = vr - fr[it + 1];
        //fi[it + 1] = vi - fi[it + 1];
    }

    m = n / 2;
    nv = 2;
    for (l0 = k - 2; l0 >= 0; l0--) //蝴蝶操作
    {
        m = m / 2;
        nv = 2 * nv;
        for (it = 0; it <= (m - 1) * nv; it = it + nv)
            for (j = 0; j <= (nv / 2) - 1; j++)
            {
                p = pr[m * j] * fr[it + j + nv / 2];
                q = pi[m * j] * fi[it + j + nv / 2];
                s = (pr[m * j] + pi[m * j]) * (fr[it + j + nv / 2] + fi[it + j + nv / 2]);
                poddr = p - q;
                poddi = s - p - q;
                fr[it + j + nv / 2] = fr[it + j] - poddr;
                fi[it + j + nv / 2] = fi[it + j] - poddi;
                fr[it + j] = fr[it + j] + poddr;
                fi[it + j] = fi[it + j] + poddi;
            }
    }

    for (i = 0; i < (n>>1); i++)
    {
        mag[i] = sqrt(fr[i] * fr[i] + fi[i] * fi[i]); //幅值计算
    }

    return;
}


void kfft4(double pr[], double pi[], int n, int p, double mag[])
{

    int i, j, k, k1, kw;
    int l, m, t;
    int jm, lm;

    // prepare (W_1024)^(0 ~ 1023)
    double Wreal[1024] = {0.0};
    double Wimag[1024] = {0.0};
    double PI2n = 3.141592653589793 * 2 / n;
    double atr, ati, pt, qt, st;

    Wreal[0] = 1.0;
    Wimag[0] = 0.0;
    Wreal[1] = cos(PI2n);
    Wimag[1] = -sin(PI2n);

    for (i = 2; i < n; i++) {
        pt = Wreal[i - 1] * Wreal[1];
        qt = Wimag[i - 1] * Wimag[1];
        st = (Wreal[i - 1] + Wimag[i - 1]) * (Wreal[1] + Wimag[1]);
        Wreal[i] = pt - qt;
        Wimag[i] = st - pt - qt;
    }

    double c0r, c1r, c2r, c3r;
    double c0i, c1i, c2i, c3i;
    double d0r, d1r, d2r, d3r;
    double d0i, d1i, d2i, d3i;

    double *pX0r, *pX0i;
    double *pX1r, *pX1i;
    double *pTemp;

    double fr[1024] = {0.0};
    double fi[1024] = {0.0};

    pX0r = pr;
    pX0i = pi;

    pX1r = fr;
    pX1i = fi;

    l = (n >> 2);
    m = 1;

    for (t = 0; t < (p>>1); t++) {
        for (j = 0; j < l; j++) {
            for (k = 0; k < m; k++) {
                jm = j * m;
                lm = l * m;

                k1 = k + jm;
                c0r = pX0r[k1];
                c0i = pX0i[k1];

                k1 += lm;
                c1r = pX0r[k1];
                c1i = pX0i[k1];

                k1 += lm;
                c2r = pX0r[k1];
                c2i = pX0i[k1];

                k1 += lm;
                c3r = pX0r[k1];
                c3i = pX0i[k1];

                d0r = c0r + c2r;
                d0i = c0i + c2i;

                d1r = c0r - c2r;
                d1i = c0i - c2i;

                d2r = c1r + c3r;
                d2i = c1i + c3i;

                d3r = c1i - c3i;
                d3i = c3r - c1r;

                k1 = k + (jm<<2);
                pX1r[k1] = d0r + d2r;
                pX1i[k1] = d0i + d2i;

                k1 += m;
                kw = jm;
                atr = d1r + d3r;
                ati = d1i + d3i;
                pt = atr * Wreal[kw];
                qt = ati * Wimag[kw];
                st = (atr + ati) * (Wreal[kw] + Wimag[kw]);
                pX1r[k1] = pt - qt;
                pX1i[k1] = st - pt - qt;

                k1 += m;
                kw += jm;
                atr = d0r - d2r;
                ati = d0i - d2i;
                pt = atr * Wreal[kw];
                qt = ati * Wimag[kw];
                st = (atr + ati) * (Wreal[kw] + Wimag[kw]);
                pX1r[k1] = pt - qt;
                pX1i[k1] = st - pt - qt;

                k1 += m;
                kw += jm;
                atr = d1r - d3r;
                ati = d1i - d3i;
                pt = atr * Wreal[kw];
                qt = ati * Wimag[kw];
                st = (atr + ati) * (Wreal[kw] + Wimag[kw]);
                pX1r[k1] = pt - qt;
                pX1i[k1] = st - pt - qt;
            }
        }
        l >>= 2;
        m <<= 2;

        pTemp = pX0r;
        pX0r = pX1r;
        pX1r = pTemp;

        pTemp = pX0i;
        pX0i = pX1i;
        pX1i = pTemp;
    }

    for (i = 0; i < (n>>1); i++)
    {
        mag[i] = sqrt(pX0r[i] * pX0r[i] + pX0i[i] * pX0i[i]); //幅值计算
    }

    return;
}



int main()
{
    clock_t a, b;
	int i,j;
    double pr[1024],pi[1024],mag[1024],t[1024];

    for (i=0; i<1024; i++)  //生成输入信号
    {
		t[i] = i*0.001;
		pr[i]=1.2+2.7*cos(2*PI*33*t[i])+5*cos(2*PI*200*t[i]+PI/2);
		pi[i]=0.0;
	}

    a = clock();
    kfft4(pr,pi,1024,10,mag);  //调用FFT函数
    b = clock();
    printf("kfft_time = %ld\n", (b-a));


	for (i=0; i<4; i++)
    {
        printf("%d\t%lf\n",i,mag[i]); //输出结果
    }


    return 0;
}
