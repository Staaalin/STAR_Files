#include "math.h"
#include "string.h"
#include <iostream>
#include <map>
#include <vector>
// #include "SyMath.h"
using namespace std;

inline bool IfYType(double Y0,double Y1,double Y2){
    if (Y1 <= Y0 && Y1 < Y2){
        return true;
    }
    if (Y1 < Y0 && Y1 <= Y2){
        return true;
    }
    return false;
}

pair<std::vector<double> , std::vector<double>> NMin(double (*func)(double),
                                                     double Input_X_min,
                                                     double Input_X_max,
                                                     double Err,
                                                     int    MAX_Itr,
                                                     bool   IfGlobal = false,
                                                     int    Global_First_Bin = 10,
                                                     bool   IfExtend = false)
{
    
    vector<double> result,result_err;
    double X_min = 1.0*Input_X_min;
    double X_max = 1.0*Input_X_max;
    double X_mid = 0.5*(X_min + X_max);
    double X_Tem,Y_Tem;
    double Y_Left,Y_Right,Y_Mid;
    double X_Array[3],Y_Array[3];
    vector<double> TX_Array,TY_Array;
    vector<double> VX0_Array,VX1_Array,VX2_Array;
    vector<double> VY0_Array,VY1_Array,VY2_Array;
    result.resize(0);result_err.resize(0);
    VX0_Array.resize(0);VY0_Array.resize(0);
    VX1_Array.resize(0);VY1_Array.resize(0);
    VX2_Array.resize(0);VY2_Array.resize(0);

    if (!IfGlobal) {
        Global_First_Bin = 3;
    }
    double Scan_Step = (X_max - X_min)/(Global_First_Bin - 1);
    TX_Array.resize(Global_First_Bin);TY_Array.resize(Global_First_Bin);
    for(int Itr=0;Itr < Global_First_Bin;Itr++){
        TX_Array[Itr] = X_min + Scan_Step * Itr;
        TY_Array[Itr] = (*func)(TX_Array[Itr]);
    }
    for(int Itr=1;Itr < Global_First_Bin - 1;Itr++){
        if (IfYType(TY_Array[Itr - 1],TY_Array[Itr],TY_Array[Itr + 1])){
            VX0_Array.emplace_back(TX_Array[Itr - 1]);VY0_Array.emplace_back(TY_Array[Itr - 1]);
            VX1_Array.emplace_back(TX_Array[Itr]);    VY1_Array.emplace_back(TY_Array[Itr]);    
            VX2_Array.emplace_back(TX_Array[Itr + 1]);VY2_Array.emplace_back(TY_Array[Itr + 1]);
        }
    }
    if (VX0_Array.size() == 0 && IfExtend){
        for(int Itr=1;Itr <= 4;Itr++){
            VX0_Array.resize(0);VY0_Array.resize(0);
            VX1_Array.resize(0);VY1_Array.resize(0);
            VX2_Array.resize(0);VY2_Array.resize(0);
            int Size_Of_TX_Array = TX_Array.size();
            TX_Array.resize(2*Size_Of_TX_Array - 1);TY_Array.resize(2*Size_Of_TX_Array - 1);
            for(int Jtr=Size_Of_TX_Array-1;Jtr >= 1;Jtr--){
                TX_Array[2*Jtr] = TX_Array[Jtr];
                TY_Array[2*Jtr] = TY_Array[Jtr];
            }
            for(int Jtr=0;Jtr < Size_Of_TX_Array - 1;Jtr++){
                TX_Array[2*Jtr+1] = 0.5*(TX_Array[2*Jtr] + TX_Array[2*Jtr+2]);
                TY_Array[2*Jtr+1] = (*func)(TX_Array[2*Jtr+1]);
                if (IfYType(TY_Array[2*Jtr],TY_Array[2*Jtr+1],TY_Array[2*Jtr+2])){
                    VX0_Array.emplace_back(TX_Array[2*Jtr])  ;VY0_Array.emplace_back(TY_Array[2*Jtr]);
                    VX1_Array.emplace_back(TX_Array[2*Jtr+1]);VY1_Array.emplace_back(TY_Array[2*Jtr+1]);    
                    VX2_Array.emplace_back(TX_Array[2*Jtr+2]);VY2_Array.emplace_back(TY_Array[2*Jtr+2]);
                    break;
                }
            }
            if(VX0_Array.size() != 0 && !IfGlobal){break;}
            for(int Jtr=1;Jtr < Size_Of_TX_Array - 1;Jtr++){
                if (IfYType(TY_Array[2*Jtr-1],TY_Array[2*Jtr],TY_Array[2*Jtr+1])){
                    VX0_Array.emplace_back(TX_Array[2*Jtr-1]);VY0_Array.emplace_back(TY_Array[2*Jtr-1]);
                    VX1_Array.emplace_back(TX_Array[2*Jtr])  ;VY1_Array.emplace_back(TY_Array[2*Jtr]);    
                    VX2_Array.emplace_back(TX_Array[2*Jtr+1]);VY2_Array.emplace_back(TY_Array[2*Jtr+1]);
                    break;
                }
            }
            if(VX0_Array.size() != 0 && !IfGlobal){break;}
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
                Itr == MAX_Itr)
            {
                result.emplace_back(0.5*(X_Array[0] + X_Array[2]));
                result_err.emplace_back(0.5*(X_Array[2] - X_Array[0]));
                // cout<<0.5*(X_Array[0] + X_Array[2])<<endl;
                break;
            }

            X_Tem = (X_Array[2]*X_Array[2]*(Y_Array[0] - Y_Array[1]) + X_Array[0]*X_Array[0]*(Y_Array[1] - Y_Array[2]) + X_Array[1]*X_Array[1]*(-Y_Array[0] + Y_Array[2]))/(2*(X_Array[2]*(Y_Array[0] - Y_Array[1]) + X_Array[0]*(Y_Array[1] - Y_Array[2]) + X_Array[1]*(-Y_Array[0] + Y_Array[2])));
            if (X_Tem == X_Array[1]){X_Tem = 0.5*(X_Array[0] + X_Array[1]);}
            if (X_Tem > X_Array[2]) {X_Tem = 0.5*(X_Array[1] + X_Array[2]);}
            if (X_Tem < X_Array[0]) {X_Tem = 0.5*(X_Array[0] + X_Array[1]);}
            Y_Tem = (*func)(X_Tem);

            // Y_Array[0] , Y_Array[1] , Y_Array[2]
            // X_Array[0] , X_Array[1] , X_Array[2]
            // cout<<"Y_Array = ["<<Y_Array[0]<<" , "<<Y_Array[1]<<" , "<<Y_Array[2]<<"]"<<endl;
            // cout<<"X_Array = ["<<X_Array[0]<<" , "<<X_Array[1]<<" , "<<X_Array[2]<<"]"<<endl;
            

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

pair<std::vector<double> , std::vector<double>> NMin(double (*func)(double),
                                                     double Input_X_Array[3],
                                                     double Input_Y_Array[3],
                                                     double Err,
                                                     int    MAX_Itr)
{
    
    vector<double> result,result_err;
    double X_Tem,Y_Tem;
    double X_Array[3],Y_Array[3];
    for(int Itr = 0;Itr < 3;Itr++){
        X_Array[Itr] = Input_X_Array[Itr];
        Y_Array[Itr] = Input_Y_Array[Itr];
    }
    result.resize(0);result_err.resize(0);

    for(int Itr=1;Itr <= MAX_Itr;Itr++){

        if (X_Array[2] - X_Array[0] < fabs(Err) || 
            Itr == MAX_Itr)
        {
            result.emplace_back(0.5*(X_Array[0] + X_Array[2]));
            result_err.emplace_back(0.5*(X_Array[2] - X_Array[0]));
            // cout<<0.5*(X_Array[0] + X_Array[2])<<endl;
            break;
        }

        X_Tem = (X_Array[2]*X_Array[2]*(Y_Array[0] - Y_Array[1]) + X_Array[0]*X_Array[0]*(Y_Array[1] - Y_Array[2]) + X_Array[1]*X_Array[1]*(-Y_Array[0] + Y_Array[2]))/(2*(X_Array[2]*(Y_Array[0] - Y_Array[1]) + X_Array[0]*(Y_Array[1] - Y_Array[2]) + X_Array[1]*(-Y_Array[0] + Y_Array[2])));
        if (X_Tem == X_Array[1]){X_Tem = 0.5*(X_Array[0] + X_Array[1]);}
        if (X_Tem > X_Array[2]) {X_Tem = 0.5*(X_Array[1] + X_Array[2]);}
        if (X_Tem < X_Array[0]) {X_Tem = 0.5*(X_Array[0] + X_Array[1]);}
        Y_Tem = (*func)(X_Tem);

        // Y_Array[0] , Y_Array[1] , Y_Array[2]
        // X_Array[0] , X_Array[1] , X_Array[2]
        // cout<<"Y_Array = ["<<Y_Array[0]<<" , "<<Y_Array[1]<<" , "<<Y_Array[2]<<"]"<<endl;
        // cout<<"X_Array = ["<<X_Array[0]<<" , "<<X_Array[1]<<" , "<<X_Array[2]<<"]"<<endl;
        

        if (Y_Tem <= Y_Array[1]){
            if (X_Tem > X_Array[1]){X_Array[0] = X_Array[1];Y_Array[0] = Y_Array[1];}
            else                   {X_Array[2] = X_Array[1];Y_Array[2] = Y_Array[1];}
            X_Array[1] = X_Tem;Y_Array[1] = Y_Tem;
        }else{
            if (X_Tem > X_Array[1]){X_Array[2] = X_Tem;Y_Array[2] = Y_Tem;}
            else                   {X_Array[0] = X_Tem;Y_Array[0] = Y_Tem;}
        }
        
    }
    return std::make_pair(result,result_err);

}

pair<std::vector<double> , std::vector<double>> NMin(double (*func)(double),
                                                     vector<double> Input_X_Array,
                                                     vector<double> Input_Y_Array,
                                                     double Err,
                                                     int    MAX_Itr,
                                                     bool   IfGlobal = false,
                                                     int    Global_First_Bin = 10,
                                                     bool   IfExtend = false)
{
    
    vector<double> result,result_err;
    result.resize(0);result_err.resize(0);

    int Input_X_Array_Size = Input_X_Array.size();
    int Itr = 0;
    while(Itr < Input_X_Array_Size - 2){
        Itr++;
        if (IfYType(Input_Y_Array[Itr-1],Input_Y_Array[Itr],Input_Y_Array[Itr+1])){
            double X_Array[3],Y_Array[3];
            X_Array[0] = Input_X_Array[Itr-1];X_Array[1] = Input_X_Array[Itr];X_Array[2] = Input_X_Array[Itr+1];
            Y_Array[0] = Input_Y_Array[Itr-1];Y_Array[1] = Input_Y_Array[Itr];Y_Array[2] = Input_Y_Array[Itr+1];
            pair<std::vector<double> , std::vector<double>>RV = NMin(func,X_Array,Y_Array,Err,MAX_Itr);
            result.emplace_back(RV.first[0]);result_err.emplace_back(RV.second[0]);
            Itr++;
            cout<<Y_Array[0]<<" , "<<Y_Array[1]<<" , "<<Y_Array[2]<<endl;
            cout<<X_Array[0]<<" , "<<X_Array[1]<<" , "<<X_Array[2]<<endl;
        }
    }
    return std::make_pair(result,result_err);

}

double TestFunc(double x){
    double a = 20.0;
    // return 2*(x - a)*(x - a)*(x - a)*(x - a) + 12*(x - a)*(x - a) + 8;
    return -cos(2*3.1415926*(x - a));
}

void SyMath(){
    // pair<std::vector<double> , std::vector<double>>RV = NMin(TestFunc,0.0,21.0,0.001,40,false,10,false);

    // double Input_X_Array[3] = {0.0 , 19.0 , 21.0},Input_Y_Array[3];
    // for(int Itr = 0;Itr < 3;Itr++){
    //     Input_Y_Array[Itr] = TestFunc(Input_X_Array[Itr]);
    // }
    // pair<std::vector<double> , std::vector<double>>RV = NMin(TestFunc,Input_X_Array,Input_Y_Array,0.001,40);

    vector<double> Input_X_Array,Input_Y_Array;
    Input_X_Array.emplace_back(0.0);Input_X_Array.emplace_back(19.0);Input_X_Array.emplace_back(21.0);Input_X_Array.emplace_back(22.0);Input_X_Array.emplace_back(22.2);Input_X_Array.emplace_back(45.0);
    cout<<"Input_X_Array = ";
    for(int Itr = 0;Itr < Input_X_Array.size();Itr++){
        cout<<Input_X_Array[Itr]<<" , ";
    }
    cout<<endl;
    for(int Itr = 0;Itr < Input_X_Array.size();Itr++){
        Input_Y_Array.emplace_back(TestFunc(Input_X_Array[Itr]));
    }
    cout<<"Input_Y_Array = ";
    for(int Itr = 0;Itr < Input_X_Array.size();Itr++){
        cout<<Input_Y_Array[Itr]<<" , ";
    }
    cout<<endl;
    pair<std::vector<double> , std::vector<double>>RV = NMin(TestFunc,Input_X_Array,Input_Y_Array,0.001,40,true,10,false);

    std::cout<<"FOUND NUMBER OF MIN: "<<RV.first.size()<<endl;
    if (RV.first.size()!=0){
        for(int i=0;i < RV.first.size();i++){
            std::cout<<RV.first[i]<<endl;
        }
    }
}