#include "math.h"
#include "string.h"
#include <iostream>
#include <map>
#include <vector>
// #include "SyMath.h"
using namespace std;
pair<std::vector<float> , std::vector<float>> NMin(float (*func)(float),
                                                   float Input_X_min,
                                                   float Input_X_max,
                                                   float Err,
                                                   int   MAX_Itr,
                                                   bool  IfGlobal = false,
                                                   int   Global_First_Bin = 10,
                                                   bool  IfExtend = false)
{
    
    vector<float> result,result_err;
    float X_min = 1.0*Input_X_min;
    float X_max = 1.0*Input_X_max;
    float X_mid = 0.5*(X_min + X_max);
    float X_Tem,Y_Tem;
    float Y_Left,Y_Right,Y_Mid;
    float X_Array[3],Y_Array[3];
    vector<float> TX_Array,TY_Array;
    vector<float> VX0_Array,VX1_Array,VX2_Array;
    vector<float> VY0_Array,VY1_Array,VY2_Array;
    result.resize(0);result_err.resize(0);
    VX0_Array.resize(0);VY0_Array.resize(0);
    VX1_Array.resize(0);VY1_Array.resize(0);
    VX2_Array.resize(0);VY2_Array.resize(0);

    if (!IfGlobal) {
        Global_First_Bin = 3;
    }
    float Scan_Step = (X_max - X_min)/(Global_First_Bin - 1);
    TX_Array.resize(Global_First_Bin);TY_Array.resize(Global_First_Bin);
    for(int Itr=0;Itr < Global_First_Bin;Itr++){
        TX_Array[Itr] = X_min + Scan_Step * Itr;
        TY_Array[Itr] = (*func)(TX_Array[Itr]);
    }
    for(int Itr=1;Itr < Global_First_Bin - 1;Itr++){
        if (TY_Array[Itr - 1] >= TY_Array[Itr] && TY_Array[Itr + 1] > TY_Array[Itr]){
            VX0_Array.emplace_back(TX_Array[Itr - 1]);VY0_Array.emplace_back(TY_Array[Itr - 1]);
            VX1_Array.emplace_back(TX_Array[Itr]);    VY1_Array.emplace_back(TY_Array[Itr]);    
            VX2_Array.emplace_back(TX_Array[Itr + 1]);VY2_Array.emplace_back(TY_Array[Itr + 1]);
        }
    }
    if (VX0_Array.size() == 0 && !IfGlobal && IfExtend){
        for(int Itr=1;Itr <= 4;Itr++){
            int Size_Of_TX_Array = TX_Array.size();
            TX_Array.resize(2*Size_Of_TX_Array - 1);TY_Array.resize(2*Size_Of_TX_Array - 1);
            for(int Jtr=Size_Of_TX_Array-1;Jtr >= 1;Jtr--){
                TX_Array[2*Jtr] = TX_Array[Jtr];
                TY_Array[2*Jtr] = TY_Array[Jtr];
            }
            for(int Jtr=0;Jtr < Size_Of_TX_Array - 1;Jtr++){
                TX_Array[2*Jtr+1] = 0.5*(TX_Array[2*Jtr] + TX_Array[2*Jtr+2]);
                TY_Array[2*Jtr+1] = (*func)(TX_Array[2*Jtr+1]);
                if (TY_Array[2*Jtr] >= TY_Array[2*Jtr+1] && TY_Array[2*Jtr+2] > TY_Array[2*Jtr+1]){
                    VX0_Array.emplace_back(TX_Array[2*Jtr])  ;VY0_Array.emplace_back(TY_Array[2*Jtr]);
                    VX1_Array.emplace_back(TX_Array[2*Jtr+1]);VY1_Array.emplace_back(TY_Array[2*Jtr+1]);    
                    VX2_Array.emplace_back(TX_Array[2*Jtr+2]);VY2_Array.emplace_back(TY_Array[2*Jtr+2]);
                    break;
                }
            }
            if(VX0_Array.size() != 0){break;}
            for(int Jtr=1;Jtr < Size_Of_TX_Array - 1;Jtr++){
                if (TY_Array[2*Jtr] <= TY_Array[2*Jtr-1] && TY_Array[2*Jtr] < TY_Array[2*Jtr+1]){
                    VX0_Array.emplace_back(TX_Array[2*Jtr-1]);VY0_Array.emplace_back(TY_Array[2*Jtr-1]);
                    VX1_Array.emplace_back(TX_Array[2*Jtr])  ;VY1_Array.emplace_back(TY_Array[2*Jtr]);    
                    VX2_Array.emplace_back(TX_Array[2*Jtr+1]);VY2_Array.emplace_back(TY_Array[2*Jtr+1]);
                    break;
                }
            }
            if(VX0_Array.size() != 0){break;}
        }
    }

    for(int Jtr=0;Jtr < VX0_Array.size();Jtr++){
        X_Array[0] = VX0_Array[Jtr];
        Y_Array[0] = VY0_Array[Jtr];
        X_Array[1] = VX1_Array[Jtr];
        Y_Array[1] = VY1_Array[Jtr];
        X_Array[2] = VX2_Array[Jtr];
        Y_Array[2] = VY2_Array[Jtr];
        for(int Itr=1;Itr <= MAX_Itr;Itr++){

            if (X_Array[2] - X_Array[0] < fabs(Err) || 
            (X_Tem < X_Array[0] || X_Tem > X_Array[2]) ||
            Itr == MAX_Itr)
            {
                result.emplace_back(0.5*(X_Array[0] + X_Array[2]));
                result_err.emplace_back(0.5*(X_Array[2] - X_Array[0]));
                cout<<0.5*(X_Array[0] + X_Array[2])<<endl;
                break;
            }

            X_Tem = (X_Array[2]*X_Array[2]*(Y_Array[0] - Y_Array[1]) + X_Array[0]*X_Array[0]*(Y_Array[1] - Y_Array[2]) + X_Array[1]*X_Array[1]*(-Y_Array[0] + Y_Array[2]))/(2*(X_Array[2]*(Y_Array[0] - Y_Array[1]) + X_Array[0]*(Y_Array[1] - Y_Array[2]) + X_Array[1]*(-Y_Array[0] + Y_Array[2])));
            if (X_Tem == X_Array[1]){X_Tem = 0.5*(X_Array[0] + X_Array[1]);}
            if (X_Tem > X_Array[2]) {X_Tem = 0.5*(X_Array[1] + X_Array[2]);}
            if (X_Tem < X_Array[0]) {X_Tem = 0.5*(X_Array[0] + X_Array[1]);}
            Y_Tem = (*func)(X_Tem);

            // Y_Array[0] , Y_Array[1] , Y_Array[2]
            // X_Array[0] , X_Array[1] , X_Array[2]
            cout<<"["<<X_Array[0]<<" , "<<X_Array[1]<<" , "<<X_Array[2]<<"]"<<endl;
            

            if (Y_Tem <= Y_Array[1]){
                if (X_Tem > X_Array[1]){X_Array[0] = X_Array[1];Y_Array[0] = Y_Array[1];}
                else                   {X_Array[2] = X_Array[1];Y_Array[2] = Y_Array[1];}
                X_Array[1] = X_Tem;Y_Array[1] = Y_Tem;
            }else{
                if (X_Tem > X_Array[1]){X_Array[2] = X_Tem;Y_Array[2] = Y_Tem;}
                else                   {X_Array[0] = X_Tem;Y_Array[0] = Y_Tem;}
            }
            
        }
    }
    return std::make_pair(result,result_err);

}

float TestFunc(float x){
    float a = 14;
    return 2*(x - a)*(x - a)*(x - a)*(x - a) + 12*(x - a)*(x - a) + 8;
}

void SyMath(){
    cout<<"Here"<<endl;
    pair<std::vector<float> , std::vector<float>>RV = NMin(TestFunc,0.0,20.0,0.001,40);
    std::cout<<RV.first.size()<<endl;
    std::cout<<RV.second.size()<<endl;
}