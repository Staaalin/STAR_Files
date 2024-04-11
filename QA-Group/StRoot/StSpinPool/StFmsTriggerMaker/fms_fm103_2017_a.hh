#ifndef FMS_FM103_2017_A_HH
#define FMS_FM103_2017_A_HH

#include "Board.hh"

void fms_fm103_2017_a(Board& fm103, int t=MAXPP, int simdat=0);

int getFM103_2017a_BS3(int out);
int getFM103_2017a_BS2(int out);
int getFM103_2017a_BS1T(int out);
int getFM103_2017a_BS1M(int out);
int getFM103_2017a_BS1B(int out);
int getFM103_2017a_JpT(int out);
int getFM103_2017a_JpM(int out);
int getFM103_2017a_JpB(int out);

#endif	// FMS_FM103_2017_A_HH
