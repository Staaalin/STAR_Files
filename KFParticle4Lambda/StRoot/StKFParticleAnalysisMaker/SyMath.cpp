#include "KaonPID.h"
#include <cmath>

#include "NMin.h"

float NMin(float (*func)(float),
           float Input_X_min,
           float Input_X_max,
           int   MAX_Itr)
    {
    
    float X_min = 1.0*Input_X_min;
    float X_max = 1.0*Input_X_max;
    float X_mid = 0.5*(X_min + X_max);
    float X_Tem;
    float Y_Left,Y_Right,Y_Mid;
    Y_Left  = (*func)(X_min);
    Y_Mid   = (*func)(X_mid);
    Y_Right = (*func)(X_max);
    for(int Itr=1;Itr <= MAX_Itr;Itr++){

        if (Y_Left > Y_Mid && Y_Right > Y_Mid){
            X_Tem = - ;
        }
        
        Y_Left  = (*func)(X_min);
        Y_Mid   = (*func)(X_mid);
        Y_Right = (*func)(X_max);
    }

    }