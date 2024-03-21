// #include "math.h"
// #include "string.h"
// #include <iostream>
// #include <map>
// #include <vector>
pair<std::vector<double> , std::vector<double>> NMin(double (*func)(double),
                                                     double Input_X_min,
                                                     double Input_X_max,
                                                     double Err,
                                                     int    MAX_Itr,
                                                     bool   IfGlobal = false,
                                                     int    Global_First_Bin = 10,
                                                     bool   IfExtend = false) const;
pair<std::vector<double> , std::vector<double>> NMin(double (*func)(double),
                                                     double Input_X_Array[3],
                                                     double Input_Y_Array[3],
                                                     double Err,
                                                     int    MAX_Itr,
                                                     bool   IfGlobal = false,
                                                     bool   IfExtend = false)