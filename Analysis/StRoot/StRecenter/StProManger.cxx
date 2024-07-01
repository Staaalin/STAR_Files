#include "StProManger.h"
#include "StConstants.h"
#include "TProfile2D.h"
#include "TProfile.h"
#include "TString.h"
#include "TMath.h"
#include "TH1F.h"

ClassImp(StProManger)

TString StProManger::mVStr[2] = {"pos","neg"};

//---------------------------------------------------------------------------------

StProManger::StProManger()
{
}

//---------------------------------------------------------------------------------

StProManger::~StProManger()
{
}


//----------------------------------------------------------------------------

void StProManger::InitShift()
{
    TString ProName;
    for(Int_t i = 0; i < 2; i++) // vertex pos/neg
    {
        for(Int_t j = 0; j < 4; j++) // eta_gap
        {
            for(Int_t k = 0; k < 5; k++) // Shift Order
            {
                ProName = Form("CosPsi2_Vertex_%s_EtaGap_%d_Order_%d_East_EP",mVStr[i].Data(),j,k);
                p_mcos2_East_EP[i][j][k] = new TProfile2D(ProName.Data(),ProName.Data(),3000,-0.5,2999.5,9,-0.5,8.5); // x axis is RunIndex, y axis is Centrality
                ProName = Form("SinPsi2_Vertex_%s_EtaGap_%d_Order_%d_East_EP",mVStr[i].Data(),j,k);
                p_msin2_East_EP[i][j][k] = new TProfile2D(ProName.Data(),ProName.Data(),3000,-0.5,2999.5,9,-0.5,8.5); // x axis is RunIndex, y axis is Centrality
                ProName = Form("CosPsi2_Vertex_%s_EtaGap_%d_Order_%d_West_EP",mVStr[i].Data(),j,k);
                p_mcos2_West_EP[i][j][k] = new TProfile2D(ProName.Data(),ProName.Data(),3000,-0.5,2999.5,9,-0.5,8.5); // x axis is RunIndex, y axis is Centrality
                ProName = Form("SinPsi2_Vertex_%s_EtaGap_%d_Order_%d_West_EP",mVStr[i].Data(),j,k);
                p_msin2_West_EP[i][j][k] = new TProfile2D(ProName.Data(),ProName.Data(),3000,-0.5,2999.5,9,-0.5,8.5); // x axis is RunIndex, y axis is Centrality

                ProName = Form("CosPsi3_Vertex_%s_EtaGap_%d_Order_%d_East_EP",mVStr[i].Data(),j,k);
                p_mcos3_East_EP[i][j][k] = new TProfile2D(ProName.Data(),ProName.Data(),3000,-0.5,2999.5,9,-0.5,8.5); // x axis is RunIndex, y axis is Centrality
                ProName = Form("SinPsi3_Vertex_%s_EtaGap_%d_Order_%d_East_EP",mVStr[i].Data(),j,k);
                p_msin3_East_EP[i][j][k] = new TProfile2D(ProName.Data(),ProName.Data(),3000,-0.5,2999.5,9,-0.5,8.5); // x axis is RunIndex, y axis is Centrality
                ProName = Form("CosPsi3_Vertex_%s_EtaGap_%d_Order_%d_West_EP",mVStr[i].Data(),j,k);
                p_mcos3_West_EP[i][j][k] = new TProfile2D(ProName.Data(),ProName.Data(),3000,-0.5,2999.5,9,-0.5,8.5); // x axis is RunIndex, y axis is Centrality
                ProName = Form("SinPsi3_Vertex_%s_EtaGap_%d_Order_%d_West_EP",mVStr[i].Data(),j,k);
                p_msin3_West_EP[i][j][k] = new TProfile2D(ProName.Data(),ProName.Data(),3000,-0.5,2999.5,9,-0.5,8.5); // x axis is RunIndex, y axis is Centrality

            }
        }

        for(Int_t k = 0; k < 5; k++) // Shift Order
        {
            ProName = Form("CosPsi2_Vertex_%s_Order_%d_Full_EP",mVStr[i].Data(),k);
            p_mcos2_Full_EP[i][k] = new TProfile2D(ProName.Data(),ProName.Data(),3000,-0.5,2999.5,9,-0.5,8.5); // x axis is RunIndex, y axis is Centrality
            ProName = Form("SinPsi2_Vertex_%s_Order_%d_Full_EP",mVStr[i].Data(),k);
            p_msin2_Full_EP[i][k] = new TProfile2D(ProName.Data(),ProName.Data(),3000,-0.5,2999.5,9,-0.5,8.5); // x axis is RunIndex, y axis is Centrality

            ProName = Form("CosPsi3_Vertex_%s_Order_%d_Full_EP",mVStr[i].Data(),k);
            p_mcos3_Full_EP[i][k] = new TProfile2D(ProName.Data(),ProName.Data(),3000,-0.5,2999.5,9,-0.5,8.5); // x axis is RunIndex, y axis is Centrality
            ProName = Form("SinPsi3_Vertex_%s_Order_%d_Full_EP",mVStr[i].Data(),k);
            p_msin3_Full_EP[i][k] = new TProfile2D(ProName.Data(),ProName.Data(),3000,-0.5,2999.5,9,-0.5,8.5); // x axis is RunIndex, y axis is Centrality

        }
    }
}


//----------------------------------------------------------------------------
// Event Plane method
void StProManger::FillEventEast_EP(TVector2 Psi2Vector, TVector2 Psi3Vector, Int_t Cent9, Int_t RunIndex, Int_t i, Int_t j, Int_t k) // i = vertex pos/neg, j = eta_gap, k = ShiftOrder
{
    const Float_t cos2 = Psi2Vector.X();
    const Float_t sin2 = Psi2Vector.Y();
    const Float_t cos3 = Psi3Vector.X();
    const Float_t sin3 = Psi3Vector.Y();

    p_mcos2_East_EP[i][j][k]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)cos2);
    p_msin2_East_EP[i][j][k]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)sin2);

    p_mcos3_East_EP[i][j][k]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)cos3);
    p_msin3_East_EP[i][j][k]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)sin3);
}

void StProManger::FillEventWest_EP(TVector2 Psi2Vector, TVector2 Psi3Vector, Int_t Cent9, Int_t RunIndex, Int_t i, Int_t j, Int_t k) // i = vertex pos/neg, j = eta_gap, k = ShiftOrder
{
    const Float_t cos2 = Psi2Vector.X();
    const Float_t sin2 = Psi2Vector.Y();
    const Float_t cos3 = Psi3Vector.X();
    const Float_t sin3 = Psi3Vector.Y();

    p_mcos2_West_EP[i][j][k]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)cos2);
    p_msin2_West_EP[i][j][k]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)sin2);

    p_mcos3_West_EP[i][j][k]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)cos3);
    p_msin3_West_EP[i][j][k]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)sin3);
}

void StProManger::FillEventFull_EP(TVector2 Psi2Vector, TVector2 Psi3Vector, Int_t Cent9, Int_t RunIndex, Int_t i, Int_t k) // i = vertex pos/neg, k = ShiftOrder
{
    const Float_t cos2 = Psi2Vector.X();
    const Float_t sin2 = Psi2Vector.Y();
    const Float_t cos3 = Psi3Vector.X();
    const Float_t sin3 = Psi3Vector.Y();

    p_mcos2_Full_EP[i][k]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)cos2);
    p_msin2_Full_EP[i][k]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)sin2);

    p_mcos3_Full_EP[i][k]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)cos3);
    p_msin3_Full_EP[i][k]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)sin3);
}

//----------------------------------------------------------------------------

void StProManger::WriteShift()
{
    for(Int_t i = 0; i < 2; i++) // vertex pos/neg
    {
        for(Int_t j = 0; j < 4; j++) // eta_gap
        {
            for(Int_t k = 0; k < 5; k++) // Shift Order
            {
                // Event Plane method
                p_mcos2_East_EP[i][j][k]->Write();
                p_msin2_East_EP[i][j][k]->Write();
                p_mcos2_West_EP[i][j][k]->Write();
                p_msin2_West_EP[i][j][k]->Write();

                p_mcos3_East_EP[i][j][k]->Write();
                p_msin3_East_EP[i][j][k]->Write();
                p_mcos3_West_EP[i][j][k]->Write();
                p_msin3_West_EP[i][j][k]->Write();

            }
        }
        for(Int_t k = 0; k < 5; k++) // Shift Order
        {
            // Event Plane method
            p_mcos2_Full_EP[i][k]->Write();
            p_msin2_Full_EP[i][k]->Write();
            p_mcos3_Full_EP[i][k]->Write();
            p_msin3_Full_EP[i][k]->Write();

        }
    }
}

//----------------------------------------------------------------------------

