#include "ProtonPID.h"
#include <cmath>

const float ProtonPID::num_std = 3.5;
const float ProtonPID::nSigma_mean_proton[] = {0.3481744713145836, 0.42627956049030336, 0.4845762573372773, 0.4541828502014521, 0.4794322463583554, 0.5372607806606937, 0.5824024490081556, 0.613983210979091, 0.6364227284494539};
const float ProtonPID::nSigma_std_proton[] = {0.9839644270333957, 1.0478942656011285, 1.035594551262602, 1.043522745858929, 1.0543203307016884, 1.0439088493541977, 1.031475195705601, 1.0194402693815687, 1.0060704355761803};
const float ProtonPID::zTOF_mean_proton[] = {-0.0009490698008325743, -0.0014660377290297262, -0.0005466357052841915, -0.00011984648153294772, 6.561066064719087e-05, 0.00020529856707863408, 0.00030822910331527434, 0.0003937925608661451, 0.0003466630630411491};
const float ProtonPID::zTOF_std_proton[] = {0.04335564914608217, 0.03447789310393692, 0.021913191955711918, 0.018359466855459226, 0.01666521744699522, 0.015585142215884383, 0.01481883206539533, 0.014206826973834847, 0.013845388788371804};
const float ProtonPID::nSigma_mean_antiproton[] = {0.2624094716707337, 0.34535459577820904, 0.42087515830020217, 0.38612545989726466, 0.41143281280576977, 0.4737568543085128, 0.5244760639924474, 0.5615107791443341, 0.5902118855744951};
const float ProtonPID::nSigma_std_antiproton[] = {0.8875031320990929, 1.082392326161313, 1.0137465029293287, 1.023095737803164, 1.0334968869103895, 1.0188464313816605, 0.9995677262229304, 0.9757308616886522, 0.9499793013788619};
const float ProtonPID::zTOF_mean_antiproton[] = {-0.0019131323791304747, -0.005280054435136126, -0.0028508679855806583, -0.0023643062170443465, -0.001962437837163167, -0.001681865750799462, -0.0013014430549617014, -0.0012212218155660392, -0.0012945494559558684};
const float ProtonPID::zTOF_std_antiproton[] = {0.040577429368338865, 0.03835035865727878, 0.021104695136188633, 0.01804861104602098, 0.016244532287523598, 0.015076870878868956, 0.014100537930229218, 0.013273514449286745, 0.012325901894896357};

// bool ProtonPID::IsProton()
// // deprecated
// {
//     int pTbin = static_cast<int>(floor(pT / 0.2));
//     if (pTbin < 1 || pTbin > 8) return false;

//     // rectangular 2D cut
//     if (zTOF > zTOF_mean[pTbin-1] + num_std * zTOF_std[pTbin-1] || zTOF < zTOF_mean[pTbin-1] - num_std * zTOF_std[pTbin-1]) return false;
//     if (nSigma > nSigma_mean[pTbin-1] + num_std * nSigma_std[pTbin-1] || nSigma < nSigma_mean[pTbin-1] - num_std * nSigma_std[pTbin-1]) return false;

//     // decision boundary cut
//     float x = nSigma;
//     if (pTbin == 2)
//     {
//         if (zTOF < -0.00740428142528834*x - 0.309396675573395*sqrt(4.89619390913469e-5*x*x - 0.0608378198738696*x + 1) + 0.237735090905258) return false;
//     }
//     if (pTbin == 3)
//     {
//         if (zTOF < -0.0281648248295252*x - 0.962886729723417*sqrt(0.000443038401028118*x*x - 0.0733453661219122*x + 1) + 0.900787586767528) return false;
//     } 
//     if (pTbin == 4)
//     {
//         if (zTOF > -0.0032406934598125*x + 0.13698443471894*sqrt(0.00280927391482182*x*x + 0.00221847237193595*x + 1) - 0.0716263807099521) return false;
//         if (zTOF < -0.00336998016087225*x - 0.068836431881989*sqrt(-0.00207507565876224*x*x - 0.209026087342191*x + 1) + 0.036153467837767) return false;
//     }
//     if (pTbin == 5)
//     {
//         if (zTOF > -0.00213967993133456*x + 0.0676440825750894*sqrt(0.00250094845330029*x*x - 0.00751326905926489*x + 1) - 0.0243909416305087) return false;
//         if (zTOF < -0.0038428148130518*x - 0.0747002505703554*sqrt(-0.000587438544893784*x*x - 0.207476923610012*x + 1) + 0.0460501491083069) return false;
//     }
//     if (pTbin == 6)
//     {
//         if (zTOF > -0.00175411118868176*x + 0.045714123736958*sqrt(0.000712016627218395*x*x + 0.0104467802787154*x + 1) - 0.0156521506962398) return false;
//         if (zTOF < -0.00388033959336782*x - 0.0748278805316802*sqrt(0.000117474696868625*x*x - 0.213332869998102*x + 1) + 0.0499585395905155) return false;
//     }
//     if (pTbin == 7)
//     {
//         if (zTOF > -0.00135823151446061*x + 0.0393009793559247*sqrt(-8.1076042377777e-6*x*x + 0.0375948173484308*x + 1) - 0.0158373005570447) return false;
//         if (zTOF < -0.00311068719528727*x - 0.0700715902989559*sqrt(-0.000765659313779474*x*x - 0.217735809837499*x + 1) + 0.0482913994238915) return false;
//     }
//     if (pTbin == 8)
//     {
//         if (zTOF > -0.000810505839619179*x + 0.035206338530538*sqrt(-0.00195961733178985*x*x + 0.0541320905456119*x + 1) - 0.0164687958929413) return false;
//         if (zTOF < -0.00191783901912651*x - 0.0639157056402445*sqrt(-0.00304299323779905*x*x - 0.221235277690126*x + 1) + 0.0443158469801602) return false;
//     }

//     return true;
// }

bool ProtonPID::IsProtonSimple(float nSigmaCut, int charge)
{
    int pTbin = static_cast<int>(floor(pT / 0.2));
    if (pTbin < 1 || pTbin > 9) return false;
    if (std::abs(charge) != 1) return false;

    // loose nSigma cut
    if (charge > 0)
    {
        if (fabs(nSigma-nSigma_mean_proton[pTbin-1])*1.0 / nSigma_std_proton[pTbin-1] > nSigmaCut) return false;
        return true;
    }
    else
    {
        if (fabs(nSigma-nSigma_mean_antiproton[pTbin-1])*1.0 / nSigma_std_antiproton[pTbin-1] > nSigmaCut) return false;
        return true;
    }
}