pair<std::vector<float> , std::vector<float>> NMin(float (*func)(float),
                                                   float Input_X_min,
                                                   float Input_X_max,
                                                   float Err,
                                                   int   MAX_Itr,
                                                   bool  IfGlobal = false,
                                                   int   Global_First_Bin = 10,
                                                   bool  IfExtend = false) const;