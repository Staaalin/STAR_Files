#include "KaonPID.h"
#include <cmath>

#include "NMin.h"

float NMin(float (*func)(float),
           float Input_X_min,
           float Input_X_max,
           char  Method,
           int MAX_Itr)
    {
    
    float X_min = 1.0*Input_X_min;
    float X_max = 1.0*Input_X_max;
    if (Method == 'BI' || Method == 'BISECTION'){
        float X_mid,Y_mid;
        float Y_Left  = (*func)(X_min);
        float Y_Right = (*func)(X_max);
        for(int Itr=1;Itr <= MAX_Itr;Itr++){
            Y_Left  = (*func)(X_min);
            Y_Right = (*func)(X_max);

            if (Y_Left * Y_Right < 0){
                X_mid = 0.5*(X_min + X_max);
                Y_mid = (*func)(X_mid);
                if(Y_Left * Y_mid < 0){
                    X_max = X_mid;
                    Y_Right = Y_mid;
                }
            }
        }
    }

    }