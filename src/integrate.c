/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "emc.h"

double oooFunction (double e, const_t constants)
{
    double back = sqrt(2) * pow(9.11e-31, 1.5) * constants.pi / (constants.h * constants.h * constants.h) * 
                  ( sqrt(e) / ( exp( (e - constants.ef) / (constants.kt)) + 1) );
    return back;
}

double oooIntegrateLeftrect(double from, double to, double n, const_t constants, double (*func)())
{
   double h = (to-from)/n;
   if (h == 0) return 0;
   double sum = 0.0, x;
   for(x=from; x <= (to-h); x += h)
      sum += func(x, constants);
   return h*sum;
}

double oooIntegrateRightrect(double from, double to, double n, const_t constants, double (*func)())
{
   double h = (to-from)/n;
   if (h == 0) return 0;
   double sum = 0.0, x;
   for(x=from; x <= (to-h); x += h)
     sum += func(x+h, constants);
   return h*sum;
}

double oooIntegrateMidrect(double from, double to, double n, const_t constants, double (*func)())
{
   double h = (to-from)/n;
   if (h == 0) return 0;
   double sum = 0.0, x;
   for(x=from; x <= (to-h); x += h)
     sum += func(x+h/2.0, constants);
   return h*sum;
}

double oooIntegrateTrapezium(double from, double to, double n, const_t constants, double (*func)())
{
   double h = (to - from) / n;
   if (h == 0) return 0;
   double sum = func(from, constants) + func(to, constants);
   int i;
   for(i = 1;i < n;i++)
       sum += 2.0*func(from + i * h, constants);
   return  h * sum / 2.0;
}

double oooIntegrateSimpson(double from, double to, double n, const_t constants, double (*func)())
{
   double h = (to - from) / n;
   if (h == 0) return 0;
   double sum1 = 0.0;
   double sum2 = 0.0;
   int i;
 
   double x;
 
   for(i = 0;i < n;i++)
      sum1 += func(from + h * i + h / 2.0, constants);

   for(i = 1;i < n;i++)
      sum2 += func(from + h * i, constants);
 
   double back = h / 6.0 * (func(from, constants) + func(to, constants) + 4.0 * sum1 + 2.0 * sum2);

   return back;
}

double oooIntegrate(const_t constants, double narrow)
{
    double from = (constants.ef - 0.5*narrow)*constants.kt;
    double to = (constants.ef + 0.5*narrow)*constants.kt;
    double n = 50;
    return oooIntegrateSimpson(from, to, n, constants, oooFunction);
}