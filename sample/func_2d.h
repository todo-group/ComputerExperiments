#ifndef FUNC_2D_H
#define FUNC_2D_H

double f(const double* x)
{
    const double r2 = x[0] * x[0] + x[1] * x[1];
    return -10 * r2 + r2 * r2 - 2 * (x[0] + x[1]);
}

void df(const double* x, double* grad)
{
    const double r2 = x[0] * x[0] + x[1] * x[1];
    grad[0] = -2 - 20 * x[0] + 4 * x[0] * r2;
    grad[1] = -2 - 20 * x[1] + 4 * x[1] * r2;
}

const double minx = 1.6289485065147173539;
const double miny = 1.6289485065147173539;

/*
  r2[x0_,x1_] := x0^2 + x1^2
  f[x0_,x1_] := -10 * r2[x0,x1] + r2[x0,x1]^2 - 2 * (x0+x1)
  d1[x0_,x1_] := D[f[x0,x1],x0]
  d2[x0_,x1_] := D[f[x0,x1],x1]
  NSolve[{d1[x,y]==0,d2[x,y]==0},{x,y},WorkingPrecision->20]
*/

#endif // FUNC_2D_H