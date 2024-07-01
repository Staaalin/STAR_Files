double a0=-2.47672886142121, a1=1.21087977860322, a2=-0.00108715908482163, a3=5.55017591702583e-06, a4=-1.26922203954713e-08;
double b0=25.670640560464, b1=1.8212690904889, b2=-0.00420979541903353, b3=1.44007073803242e-05, b4=-2.27741752869024e-08;
double c0=-14.4809591959354, c1=0.814028899556729, c2=0.000851480786781806, c3=3.07563604740816e-08, c4=-6.31230659685062e-09;

double paraPileup[15]={
    a0, a1, a2, a3, a4,
    b0, b1, b2, b3, b4,
    c0, c1, c2, c3, c4,
};

bool isPileUp(UShort_t refMult, UShort_t btofMatched, int VzBin) 
{
  //double refmultcutmode = m_a0 + m_a1*(btofMatched) + m_a2*pow(btofMatched,2) + m_a3*pow(btofMatched,3) + m_a4*pow(btofMatched,4);
  //double refmultcutmax  = m_b0 + m_b1*(btofMatched) + m_b2*pow(btofMatched,2) + m_b3*pow(btofMatched,3) + m_b4*pow(btofMatched,4);
  //double refmultcutmin  = m_c0 + m_c1*(btofMatched) + m_c2*pow(btofMatched,2) + m_c3*pow(btofMatched,3) + m_c4*pow(btofMatched,4);

  double refmultcutmode = paraPileup[0] + paraPileup[1]*(btofMatched) + paraPileup[2]*pow(btofMatched,2) + paraPileup[3]*pow(btofMatched,3) + paraPileup[4]*pow(btofMatched,4);
  double refmultcutmax  = paraPileup[5] + paraPileup[6]*(btofMatched) + paraPileup[7]*pow(btofMatched,2) + paraPileup[8]*pow(btofMatched,3) + paraPileup[9]*pow(btofMatched,4);
  double refmultcutmin  = paraPileup[10] + paraPileup[11]*(btofMatched) + paraPileup[12]*pow(btofMatched,2) + paraPileup[13]*pow(btofMatched,3) + paraPileup[14]*pow(btofMatched,4);

  return ( refMult > refmultcutmax || refMult < refmultcutmin );
}