#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#define MIN(a, b) (((a) < (b)) ? (a) : (b))


typedef void (*dFuncType)(double*, double*, double);

int ode45t(dFuncType dfp, double xv[], double beta, double rtol, double atol);
void dxdv(double xv[], double xvdot[], double beta);
void vdp(double xv[], double xvdot[], double beta);

int main()
{
    int ret = 0;

    //double beta = 3.0;
    //double xv[2] = {0.0, beta};
    //ret = ode45t(dxdv, xv, beta, 1E-6, 1E-8);

    double beta = 1.0;
    double xv[2] = {2.0, 0.0};
    ret = ode45t(vdp, xv, beta, 1E-6, 1E-8);

    printf("Return of ode45t is %d\n", ret);

    return 0;
}


int ode45t(dFuncType dfp, double xv[], double beta, double rtol, double atol)
{
    int ret = 0;

    /*** coefficients ***/
    double pw = 1.0/5;

    double a21 = 1.0/5;

    double a31 = 3.0/40;
    double a32 = 9.0/40;

    double a41 = 44.0/45;
    double a42 = -56.0/15;
    double a43 = 32.0/9;

    double a51 = 19372.0/6561;
    double a52 = -25360.0/2187;
    double a53 = 64448.0/6561;
    double a54 = -212.0/729;

    double a61 = 9017.0/3168;
    double a62 = -355.0/33;
    double a63 = 46732.0/5247;
    double a64 = 49.0/176;
    double a65 = -5103.0/18656;

    double b1 = 35.0/384;
    double b3 = 500.0/1113;
    double b4 = 125.0/192;
    double b5 = -2187.0/6784;
    double b6 = 11.0/84;

    double e1 = 71.0/57600;
    double e3 = -71.0/16695;
    double e4 = 71.0/1920;
    double e5 = -17253.0/339200;
    double e6 = 22.0/525;
    double e7 = -1.0/40;

    //*** ------Initialize-------***//
    double k1[2], k2[2], k3[2], k4[2], k5[2], k6[2], k7[2];
    double xc[2] = {xv[0], xv[1]};
    double xtemp[2] = {0, 0};
    double xnew[2] = {0, 0};
    double xdot[2] = {0, 0};
    double ts = 0.0;

    double fE[2] = {0, 0};
    double err = 0.0;
    double max_temp[2] = {0, 0};
    double temp = 0.0;

    double threshold = atol / rtol;

    //***  hmin and hmax  ***//
    double hmin = 1E-14;
    double hmax = 0.1;

    // initialize time step: h
    dfp(xc, xdot, beta);
    double rh0 = fabs(xdot[0]) / MAX(fabs(xc[0]), threshold);
    double rh1 = fabs(xdot[1]) / MAX(fabs(xc[1]), threshold);

    double h = 0.8 * pow(rtol, pw) / MAX(rh0, rh1);
    h = MAX(h, hmin);
    h = MIN(h, hmax);


    //*** print for debug ***//
    FILE *gFPlog;
    gFPlog = fopen("ode_log.txt", "w");

    fprintf(gFPlog, "%.15f  %.15f  %.15f\n", xc[0], xc[1], ts);
    //***  -----end-----  ***//

    for (int ii = 0; ii < 200; ii++) {

        dfp(xc, xdot, beta);
        k1[0] = xdot[0];
        k1[1] = xdot[1];

        for (;;) {
            xtemp[0] = xc[0] + h*a21*k1[0];
            xtemp[1] = xc[1] + h*a21*k1[1];
            dfp(xtemp, xdot, beta);
            k2[0] = xdot[0];
            k2[1] = xdot[1];

            xtemp[0] = xc[0] + h * (a31*k1[0] + a32*k2[0]);
            xtemp[1] = xc[1] + h * (a31*k1[1] + a32*k2[1]);
            dfp(xtemp, xdot, beta);
            k3[0] = xdot[0];
            k3[1] = xdot[1];

            xtemp[0] = xc[0] + h * (a41*k1[0] + a42*k2[0] + a43*k3[0]);
            xtemp[1] = xc[1] + h * (a41*k1[1] + a42*k2[1] + a43*k3[1]);
            dfp(xtemp, xdot, beta);
            k4[0] = xdot[0];
            k4[1] = xdot[1];

            xtemp[0] = xc[0] + h * (a51*k1[0] + a52*k2[0] + a53*k3[0] + a54*k4[0]);
            xtemp[1] = xc[1] + h * (a51*k1[1] + a52*k2[1] + a53*k3[1] + a54*k4[1]);
            dfp(xtemp, xdot, beta);
            k5[0] = xdot[0];
            k5[1] = xdot[1];

            xtemp[0] = xc[0] + h * (a61*k1[0] + a62*k2[0] + a63*k3[0] + a64*k4[0] + a65*k5[0]);
            xtemp[1] = xc[1] + h * (a61*k1[1] + a62*k2[1] + a63*k3[1] + a64*k4[1] + a65*k5[1]);
            dfp(xtemp, xdot, beta);
            k6[0] = xdot[0];
            k6[1] = xdot[1];

            //xnew
            xnew[0] = xc[0] + h * (b1*k1[0] + b3*k3[0] + b4*k4[0] + b5*k5[0] + b6*k6[0]);
            xnew[1] = xc[1] + h * (b1*k1[1] + b3*k3[1] + b4*k4[1] + b5*k5[1] + b6*k6[1]);
            dfp(xnew, xdot, beta);
            k7[0] = xdot[0];
            k7[1] = xdot[1];

            //err
            fE[0] = e1*k1[0] + e3*k3[0] + e4*k4[0] + e5*k5[0] + e6*k6[0] + e7*k7[0];
            fE[1] = e1*k1[1] + e3*k3[1] + e4*k4[1] + e5*k5[1] + e6*k6[1] + e7*k7[1];

            max_temp[0] = MAX(fabs(xc[0]), fabs(xnew[0]));
            max_temp[0] = MAX(max_temp[0], threshold);
            max_temp[0] = fabs(fE[0]) / max_temp[0];

            max_temp[1] = MAX(fabs(xc[1]), fabs(xnew[1]));
            max_temp[1] = MAX(max_temp[1], threshold);
            max_temp[1] = fabs(fE[1]) / max_temp[1];

            err = h * MAX(max_temp[0], max_temp[1]);

            if (err > rtol) {
                temp = 0.8 * pow(rtol/err, pw);
                temp = h * MAX(0.1, temp);

                if (temp > hmin) {
                    h = temp;
                }
                else {
                    // cannot shrink h below hmin
                    if (ret == 2) {
                        return ret;
                    }

                    h = hmin;
                    ret = 2;
                }
            }
            else {
                ts += h;

                //***  print for debug  ***//
                fprintf(gFPlog, "%.15f  %.15f  %.15f\n", xnew[0], xnew[1], ts);

                // stop event
                if (ts > 33.0) {

                    //***  print for debug  ***//
                    fclose(gFPlog);

                    return 0;
                }

                // refresh xc[]
                xc[0] = xnew[0];
                xc[1] = xnew[1];                

                // loosen h
                temp = 1.25 * pow(err/rtol, pw);
                if (temp > 0.2) {
                    h = h / temp;
                }
                else {
                    h = 5.0 * h;
                }

                h = MIN(h, hmax);

                break;
            }
        }

    }

    //***  print for debug  ***//
    fclose(gFPlog);

    return 1;
}

void dxdv(double xv[], double xvdot[], double beta)
{
    xvdot[0] = beta * xv[1];
    xvdot[1] = -beta * xv[0];
}

void vdp(double xv[], double xvdot[], double beta)
{
    xvdot[0] = xv[1];
    xvdot[1] = beta * (1 - xv[0]*xv[0]) * xv[1] - xv[0];
}


