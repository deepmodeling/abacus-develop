/**
 * Intro
 * -----
 * This file stores XC dependent DFT-D3 parameters for Grimme-D3
 * dispersion correction.
 * 
 * Supported forms:
 * 
 * DFT-D3(0): zero-damping
 * DFT-D3(BJ): Becke-Johnson damping
 * DFT-D3M(0): zero-damping with modified damping function
 * DFT-D3M(BJ): Becke-Johnson damping with modified damping function
 * 
 * A detailed introduction of undamped, and BJ damping, the modified 
 * damping can be found in DFT-D3 software manual, see:
 * https://www.chemie.uni-bonn.de/grimme/de/software/dft-d3/man.pdf
 *
 * Other excellent learning materials (where you can find expression
 * of both DFT-D2 and DFT-D3):
 * DFT-D2: https://www.vasp.at/wiki/index.php/DFT-D2
 * DFT-D3: https://www.vasp.at/wiki/index.php/DFT-D3
 * 
 * Usage
 * -----
 * call function DFTD3::search(xc, method, param) to get the DFT-D3 parameters
 * for the given XC functional. The obtained param should be a std::vector<double>,
 * in which the first 9 elements are the DFT-D3 parameters:
 * 's6', 'sr6', 'a1', 's8', 'sr8', 'a2', 's9', 'alp', 'bet'
 * 
 * ParamNotFoundError
 * ------------------
 * If the requested D3 parameters of XC are not found, then the ABACUS will 
 * WARNING_QUIT with the message "DFT-D3 parameters for XC not found".
 * 
 * Other dispersion correction
 * ---------------------------
 * there are other kinds of dispersion correction, such as the xc VV09, VV10,
 * and rVV10, and the vdw-DF family nonlocal dispersion correction. They will
 * be mixed directly with the correlation and exchange part, which act 
 * differently from the DFT-D2 and D3 methods.
 * 
 * Special: Omega-B97 family
 * -------------------------
 * (thanks for help and discussion with @hhebrewsnabla and @moolawooda)
 * wB97 XC family is special, their DFT-D3 supports are quite complicated.
 * 
 * wB97          long-range exx with B97
 * wB97X         wB97 with additional short-range exx
 * wB97X-D       wB97X_D from libXC with DFTD2, not in DFTD3 framework
 * wB97X-D3      wB97X_D3 from libXC with DFTD3(0)
 * wB97X-D3(BJ)  wB97X_V from libXC with DFTD3(BJ)
 * wB97X-V       with VV10, not in DFTD3 framework
 * wB97M-V       with VV10, not in DFTD3 framework
 * 
 * Recommended: http://bbs.keinsci.com/thread-19076-1-1.html
 * Related information from Pyscf Github repo: 
 * https://github.com/pyscf/pyscf/issues/2069
 * 
 */
#include <map>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <cassert>
#include <algorithm>
#include "module_base/formatter.h"
#include "module_base/tool_quit.h"
#include "module_hamilt_general/module_vdw/vdwd3_parameters.h"

// DFT-D3(BJ)
const std::map<std::string, std::vector<double>> bj = {
    {"bp", {1.0, 0.3946, 0.3946, 3.2822, 4.8516, 4.8516, 1.0, 14.0, 0.0}},
    {"blyp", {1.0, 0.4298, 0.4298, 2.6996, 4.2359, 4.2359, 1.0, 14.0, 0.0}},
    {"revpbe", {1.0, 0.5238, 0.5238, 2.355, 3.5016, 3.5016, 1.0, 14.0, 0.0}},
    {"rpbe", {1.0, 0.182, 0.182, 0.8318, 4.0094, 4.0094, 1.0, 14.0, 0.0}},
    {"b97_d", {1.0, 0.5545, 0.5545, 2.2609, 3.2297, 3.2297, 1.0, 14.0, 0.0}},
    {"b973c", {1.0, 0.37, 0.37, 1.5, 4.1, 4.1, 1.0, 14.0, 0.0}},
    {"pbe", {1.0, 0.4289, 0.4289, 0.7875, 4.4407, 4.4407, 1.0, 14.0, 0.0}},
    {"rpw86pbe", {1.0, 0.4613, 0.4613, 1.3845, 4.5062, 4.5062, 1.0, 14.0, 0.0}},
    {"b3lyp", {1.0, 0.3981, 0.3981, 1.9889, 4.4211, 4.4211, 1.0, 14.0, 0.0}},
    {"tpss", {1.0, 0.4535, 0.4535, 1.9435, 4.4752, 4.4752, 1.0, 14.0, 0.0}},
    {"hf", {1.0, 0.3385, 0.3385, 0.9171, 2.883, 2.883, 1.0, 14.0, 0.0}},
    {"tpss0", {1.0, 0.3768, 0.3768, 1.2576, 4.5865, 4.5865, 1.0, 14.0, 0.0}},
    {"pbe0", {1.0, 0.4145, 0.4145, 1.2177, 4.8593, 4.8593, 1.0, 14.0, 0.0}},
    {"hse06", {1.0, 0.383, 0.383, 2.31, 5.685, 5.685, 1.0, 14.0, 0.0}},
    {"hse", {1.0, 0.383, 0.383, 2.31, 5.685, 5.685, 1.0, 14.0, 0.0}}, // ABACUS implements HSE06 as HSE
    {"revpbe38", {1.0, 0.4309, 0.4309, 1.476, 3.9446, 3.9446, 1.0, 14.0, 0.0}},
    {"pw6b95", {1.0, 0.2076, 0.2076, 0.7257, 6.375, 6.375, 1.0, 14.0, 0.0}},
    {"b2plyp", {0.64, 0.3065, 0.3065, 0.9147, 5.057, 5.057, 1.0, 14.0, 0.0}},
    {"dsdblyp", {0.5, 0.0, 0.0, 0.213, 6.0519, 6.0519, 1.0, 14.0, 0.0}},
    {"dsdblypfc", {0.5, 0.0009, 0.0009, 0.2112, 5.9807, 5.9807, 1.0, 14.0, 0.0}},
    {"dodscan66", {0.3152, 0.0, 0.0, 0.0, 5.75, 5.75, 1.0, 14.0, 0.0}},
    {"revdsdblyp", {0.5451, 0.0, 0.0, 0.0, 5.2, 5.2, 1.0, 14.0, 0.0}},
    {"revdsdpbep86", {0.4377, 0.0, 0.0, 0.0, 5.5, 5.5, 1.0, 14.0, 0.0}},
    {"revdsdpbeb95", {0.3686, 0.0, 0.0, 0.0, 5.5, 5.5, 1.0, 14.0, 0.0}},
    {"revdsdpbe", {0.5746, 0.0, 0.0, 0.0, 5.5, 5.5, 1.0, 14.0, 0.0}},
    {"revdodblyp", {0.6145, 0.0, 0.0, 0.0, 5.2, 5.2, 1.0, 14.0, 0.0}},
    {"revdodpbep86", {0.477, 0.0, 0.0, 0.0, 5.5, 5.5, 1.0, 14.0, 0.0}},
    {"revdodpbeb95", {0.4107, 0.0, 0.0, 0.0, 5.5, 5.5, 1.0, 14.0, 0.0}},
    {"revdodpbe", {0.6067, 0.0, 0.0, 0.0, 5.5, 5.5, 1.0, 14.0, 0.0}},
    {"bop", {1.0, 0.487, 0.487, 3.295, 3.5043, 3.5043, 1.0, 14.0, 0.0}},
    {"mpwlyp", {1.0, 0.4831, 0.4831, 2.0077, 4.5323, 4.5323, 1.0, 14.0, 0.0}},
    {"olyp", {1.0, 0.5299, 0.5299, 2.6205, 2.8065, 2.8065, 1.0, 14.0, 0.0}},
    {"pbesol", {1.0, 0.4466, 0.4466, 2.9491, 6.1742, 6.1742, 1.0, 14.0, 0.0}},
    {"bpbe", {1.0, 0.4567, 0.4567, 4.0728, 4.3908, 4.3908, 1.0, 14.0, 0.0}},
    {"opbe", {1.0, 0.5512, 0.5512, 3.3816, 2.9444, 2.9444, 1.0, 14.0, 0.0}},
    {"ssb", {1.0, -0.0952, -0.0952, -0.1744, 5.217, 5.217, 1.0, 14.0, 0.0}},
    {"revssb", {1.0, 0.472, 0.472, 0.4389, 4.0986, 4.0986, 1.0, 14.0, 0.0}},
    {"otpss", {1.0, 0.4634, 0.4634, 2.7495, 4.3153, 4.3153, 1.0, 14.0, 0.0}},
    {"b3pw91", {1.0, 0.4312, 0.4312, 2.8524, 4.4693, 4.4693, 1.0, 14.0, 0.0}},
    {"bhlyp", {1.0, 0.2793, 0.2793, 1.0354, 4.9615, 4.9615, 1.0, 14.0, 0.0}},
    {"revpbe0", {1.0, 0.4679, 0.4679, 1.7588, 3.7619, 3.7619, 1.0, 14.0, 0.0}},
    {"tpssh", {1.0, 0.4529, 0.4529, 2.2382, 4.655, 4.655, 1.0, 14.0, 0.0}},
    {"mpw1b95", {1.0, 0.1955, 0.1955, 1.0508, 6.4177, 6.4177, 1.0, 14.0, 0.0}},
    {"pwb6k", {1.0, 0.1805, 0.1805, 0.9383, 7.7627, 7.7627, 1.0, 14.0, 0.0}},
    {"b1b95", {1.0, 0.2092, 0.2092, 1.4507, 5.5545, 5.5545, 1.0, 14.0, 0.0}},
    {"bmk", {1.0, 0.194, 0.194, 2.086, 5.9197, 5.9197, 1.0, 14.0, 0.0}},
    {"camb3lyp", {1.0, 0.3708, 0.3708, 2.0674, 5.4743, 5.4743, 1.0, 14.0, 0.0}},
    {"lcwpbe", {1.0, 0.3919, 0.3919, 1.8541, 5.0897, 5.0897, 1.0, 14.0, 0.0}},
    {"b2gpplyp", {0.56, 0.0, 0.0, 0.2597, 6.3332, 6.3332, 1.0, 14.0, 0.0}},
    {"ptpss", {0.75, 0.0, 0.0, 0.2804, 6.5745, 6.5745, 1.0, 14.0, 0.0}},
    {"pwpb95", {0.82, 0.0, 0.0, 0.2904, 7.3141, 7.3141, 1.0, 14.0, 0.0}},
    {"hf_mixed", {1.0, 0.5607, 0.5607, 3.9027, 4.5622, 4.5622, 1.0, 14.0, 0.0}},
    {"hf_sv", {1.0, 0.4249, 0.4249, 2.1849, 4.2783, 4.2783, 1.0, 14.0, 0.0}},
    {"hf_minis", {1.0, 0.1702, 0.1702, 0.9841, 3.8506, 3.8506, 1.0, 14.0, 0.0}},
    {"b3lyp_631gd", {1.0, 0.5014, 0.5014, 4.0672, 4.8409, 4.8409, 1.0, 14.0, 0.0}},
    {"hcth120", {1.0, 0.3563, 0.3563, 1.0821, 4.3359, 4.3359, 1.0, 14.0, 0.0}},
    {"dftb3", {1.0, 0.5719, 0.5719, 0.5883, 3.6017, 3.6017, 1.0, 14.0, 0.0}},
    {"pw1pw", {1.0, 0.3807, 0.3807, 2.3363, 5.8844, 5.8844, 1.0, 14.0, 0.0}},
    {"pwgga", {1.0, 0.2211, 0.2211, 2.691, 6.7278, 6.7278, 1.0, 14.0, 0.0}},
    {"hsesol", {1.0, 0.465, 0.465, 2.9215, 6.2003, 6.2003, 1.0, 14.0, 0.0}},
    {"hf3c", {1.0, 0.4171, 0.4171, 0.8777, 2.9149, 2.9149, 1.0, 14.0, 0.0}},
    {"hf3cv", {1.0, 0.3063, 0.3063, 0.5022, 3.9856, 3.9856, 1.0, 14.0, 0.0}},
    {"pbeh3c", {1.0, 0.486, 0.486, 0.0, 4.5, 4.5, 1.0, 14.0, 0.0}},
    {"scan", {1.0, 0.538, 0.538, 0.0, 5.42, 5.42, 1.0, 14.0, 0.0}},
    {"rscan", {1.0, 0.47023427, 0.47023427, 1.08859014, 5.73408312, 5.73408312, 1.0, 14.0, 0.0}},
    {"r2scan", {1.0, 0.49484001, 0.49484001, 0.78981345, 5.73083694, 5.73083694, 1.0, 14.0, 0.0}},
    {"r2scanh", {1.0, 0.4709, 0.4709, 1.1236, 5.9157, 5.9157, 1.0, 14.0, 0.0}},
    {"r2scan0", {1.0, 0.4534, 0.4534, 1.1846, 5.8972, 5.8972, 1.0, 14.0, 0.0}},
    {"r2scan50", {1.0, 0.4311, 0.4311, 1.3294, 5.924, 5.924, 1.0, 14.0, 0.0}},
    {"wb97x_v", {1.0, 0.0, 0.0, 0.2641, 5.4959, 5.4959, 1.0, 14.0, 0.0}}, 
    // NOTE: the key `wb97x_v` directly corresonding to HYB_GGA_XC_WB97X_V, which can be further
    // employed to construct either wB97X-V with VV10, or wB97X-D3BJ with D3BJ. Here it is the D3BJ
    // parameter of wB97X-D3BJ, instead of those of wB97X-V.
    {"wb97m", {1.0, 0.566, 0.566, 0.3908, 3.128, 3.128, 1.0, 14.0, 0.0}},
    {"b97m", {1.0, -0.078, -0.078, 0.1384, 5.5946, 5.5946, 1.0, 14.0, 0.0}},
    {"pbehpbe", {1.0, 0.0, 0.0, 1.1152, 6.7184, 6.7184, 1.0, 14.0, 0.0}},
    {"xlyp", {1.0, 0.0809, 0.0809, 1.5669, 5.3166, 5.3166, 1.0, 14.0, 0.0}},
    {"mpwpw", {1.0, 0.3168, 0.3168, 1.7974, 4.7732, 4.7732, 1.0, 14.0, 0.0}},
    {"hcth407", {1.0, 0.0, 0.0, 0.649, 4.8162, 4.8162, 1.0, 14.0, 0.0}},
    {"revtpss", {1.0, 0.4326, 0.4326, 1.4023, 4.4723, 4.4723, 1.0, 14.0, 0.0}},
    {"tauhcth", {1.0, 0.0, 0.0, 1.2626, 5.6162, 5.6162, 1.0, 14.0, 0.0}},
    {"b3p", {1.0, 0.4601, 0.4601, 3.3211, 4.9858, 4.9858, 1.0, 14.0, 0.0}},
    {"b1p", {1.0, 0.4724, 0.4724, 3.5681, 4.9858, 4.9858, 1.0, 14.0, 0.0}},
    {"b1lyp", {1.0, 0.1986, 0.1986, 2.1167, 5.3875, 5.3875, 1.0, 14.0, 0.0}},
    {"mpwb1k", {1.0, 0.1474, 0.1474, 0.9499, 6.6223, 6.6223, 1.0, 14.0, 0.0}},
    {"mpw1pw", {1.0, 0.3342, 0.3342, 1.8744, 4.9819, 4.9819, 1.0, 14.0, 0.0}},
    {"mpw1kcis", {1.0, 0.0576, 0.0576, 1.0893, 5.5314, 5.5314, 1.0, 14.0, 0.0}},
    {"pbeh1pbe", {1.0, 0.0, 0.0, 1.4877, 7.0385, 7.0385, 1.0, 14.0, 0.0}},
    {"pbe1kcis", {1.0, 0.0, 0.0, 0.7688, 6.2794, 6.2794, 1.0, 14.0, 0.0}},
    {"x3lyp", {1.0, 0.2022, 0.2022, 1.5744, 5.4184, 5.4184, 1.0, 14.0, 0.0}},
    {"o3lyp", {1.0, 0.0963, 0.0963, 1.8171, 5.994, 5.994, 1.0, 14.0, 0.0}},
    {"b97_1", {1.0, 0.0, 0.0, 0.4814, 6.2279, 6.2279, 1.0, 14.0, 0.0}},
    {"b97_2", {1.0, 0.0, 0.0, 0.9448, 5.994, 5.994, 1.0, 14.0, 0.0}},
    {"b98", {1.0, 0.0, 0.0, 0.7086, 6.0672, 6.0672, 1.0, 14.0, 0.0}},
    {"hiss", {1.0, 0.0, 0.0, 1.6112, 7.3539, 7.3539, 1.0, 14.0, 0.0}},
    {"hse03", {1.0, 0.0, 0.0, 1.1243, 6.8889, 6.8889, 1.0, 14.0, 0.0}},
    {"revtpssh", {1.0, 0.266, 0.266, 1.4076, 5.3761, 5.3761, 1.0, 14.0, 0.0}},
    {"revtpss0", {1.0, 0.2218, 0.2218, 1.6151, 5.7985, 5.7985, 1.0, 14.0, 0.0}},
    {"tpss1kcis", {1.0, 0.0, 0.0, 1.0542, 6.0201, 6.0201, 1.0, 14.0, 0.0}},
    {"tauhcthhyb", {1.0, 0.0, 0.0, 0.9585, 10.1389, 10.1389, 1.0, 14.0, 0.0}},
    {"m11", {1.0, 0.0, 0.0, 2.8112, 10.1389, 10.1389, 1.0, 14.0, 0.0}},
    {"sogga11x", {1.0, 0.133, 0.133, 1.1426, 5.7381, 5.7381, 1.0, 14.0, 0.0}},
    {"n12sx", {1.0, 0.3283, 0.3283, 2.49, 5.7898, 5.7898, 1.0, 14.0, 0.0}},
    {"mn12sx", {1.0, 0.0983, 0.0983, 1.1674, 8.0259, 8.0259, 1.0, 14.0, 0.0}},
    {"mn12l", {1.0, 0.0, 0.0, 2.2674, 9.1494, 9.1494, 1.0, 14.0, 0.0}},
    {"mn15", {1.0, 2.0971, 2.0971, 0.7862, 7.5923, 7.5923, 1.0, 14.0, 0.0}},
    {"lc_whpbe", {1.0, 0.2746, 0.2746, 1.1908, 5.3157, 5.3157, 1.0, 14.0, 0.0}},
    {"mpw2plyp", {0.66, 0.4105, 0.4105, 0.6223, 5.0136, 5.0136, 1.0, 14.0, 0.0}},
    {"pw91", {1.0, 0.6319, 0.6319, 1.9598, 4.5718, 4.5718, 1.0, 14.0, 0.0}},
    {"drpa75", {0.3754, 0.0, 0.0, 0.0, 4.5048, 4.5048, 1.0, 14.0, 0.0}},
    {"scsdrpa75", {0.2528, 0.0, 0.0, 0.0, 4.505, 4.505, 1.0, 14.0, 0.0}},
    {"optscsdrpa75", {0.2546, 0.0, 0.0, 0.0, 4.505, 4.505, 1.0, 14.0, 0.0}},
    {"dsdpbedrpa75", {0.3223, 0.0, 0.0, 0.0, 4.505, 4.505, 1.0, 14.0, 0.0}},
    {"dsdpbep86drpa75", {0.3012, 0.0, 0.0, 0.0, 4.505, 4.505, 1.0, 14.0, 0.0}},
    {"dsdpbep86_2011", {0.418, 0.0, 0.0, 0.0, 5.65, 5.65, 1.0, 14.0, 0.0}},
    {"dsdsvwn5", {0.46, 0.0, 0.0, 0.0, 5.6, 5.6, 1.0, 14.0, 0.0}},
    {"dsdsp86", {0.3, 0.0, 0.0, 0.0, 5.8, 5.8, 1.0, 14.0, 0.0}},
    {"dsdslyp", {0.3, 0.0, 0.0, 0.0, 5.6, 5.6, 1.0, 14.0, 0.0}},
    {"dsdspbe", {0.4, 0.0, 0.0, 0.0, 6.0, 6.0, 1.0, 14.0, 0.0}},
    {"dsdbvwn5", {0.61, 0.0, 0.0, 0.0, 5.2, 5.2, 1.0, 14.0, 0.0}},
    {"dsdblyp_2013", {0.57, 0.0, 0.0, 0.0, 5.4, 5.4, 1.0, 14.0, 0.0}},
    {"dsdbpbe", {1.22, 0.0, 0.0, 0.0, 6.6, 6.6, 1.0, 14.0, 0.0}},
    {"dsdbp86", {0.76, 0.0, 0.0, 0.0, 6.0, 6.0, 1.0, 14.0, 0.0}},
    {"dsdbpw91", {1.14, 0.0, 0.0, 0.0, 6.5, 6.5, 1.0, 14.0, 0.0}},
    {"dsdbb95", {1.02, 0.0, 0.0, 0.0, 6.8, 6.8, 1.0, 14.0, 0.0}},
    {"dsdpbevwn5", {0.54, 0.0, 0.0, 0.0, 5.1, 5.1, 1.0, 14.0, 0.0}},
    {"dsdpbelyp", {0.43, 0.0, 0.0, 0.0, 5.2, 5.2, 1.0, 14.0, 0.0}},
    {"dsdpbe", {0.78, 0.0, 0.0, 0.0, 6.1, 6.1, 1.0, 14.0, 0.0}},
    {"dsdpbep86", {0.48, 0.0, 0.0, 0.0, 5.6, 5.6, 1.0, 14.0, 0.0}},
    {"dsdpbepw91", {0.73, 0.0, 0.0, 0.0, 6.0, 6.0, 1.0, 14.0, 0.0}},
    {"dsdpbeb95", {0.61, 0.0, 0.0, 0.0, 6.2, 6.2, 1.0, 14.0, 0.0}},
    {"dsdpbehb95", {0.58, 0.0, 0.0, 0.0, 6.2, 6.2, 1.0, 14.0, 0.0}},
    {"dsdpbehp86", {0.46, 0.0, 0.0, 0.0, 5.6, 5.6, 1.0, 14.0, 0.0}},
    {"dsdmpwlyp", {0.48, 0.0, 0.0, 0.0, 5.3, 5.3, 1.0, 14.0, 0.0}},
    {"dsdmpwpw91", {0.9, 0.0, 0.0, 0.0, 6.2, 6.2, 1.0, 14.0, 0.0}},
    {"dsdmpwp86", {0.59, 0.0, 0.0, 0.0, 5.8, 5.8, 1.0, 14.0, 0.0}},
    {"dsdmpwpbe", {0.96, 0.0, 0.0, 0.0, 6.3, 6.3, 1.0, 14.0, 0.0}},
    {"dsdmpwb95", {0.82, 0.0, 0.0, 0.0, 6.6, 6.6, 1.0, 14.0, 0.0}},
    {"dsdhsepbe", {0.79, 0.0, 0.0, 0.0, 6.1, 6.1, 1.0, 14.0, 0.0}},
    {"dsdhsepw91", {0.74, 0.0, 0.0, 0.0, 6.0, 6.0, 1.0, 14.0, 0.0}},
    {"dsdhsep86", {0.46, 0.0, 0.0, 0.0, 5.6, 5.6, 1.0, 14.0, 0.0}},
    {"dsdhselyp", {0.4, 0.0, 0.0, 0.0, 5.2, 5.2, 1.0, 14.0, 0.0}},
    {"dsdtpss", {0.72, 0.0, 0.0, 0.0, 6.5, 6.5, 1.0, 14.0, 0.0}},
    {"dsdtpssb95", {0.91, 0.0, 0.0, 0.0, 7.9, 7.9, 1.0, 14.0, 0.0}},
    {"dsdolyp", {0.93, 0.0, 0.0, 0.0, 5.8, 5.8, 1.0, 14.0, 0.0}},
    {"dsdxlyp", {0.51, 0.0, 0.0, 0.0, 5.3, 5.3, 1.0, 14.0, 0.0}},
    {"dsdxb95", {0.92, 0.0, 0.0, 0.0, 6.7, 6.7, 1.0, 14.0, 0.0}},
    {"dsdb98", {0.07, 0.0, 0.0, 0.0, 3.7, 3.7, 1.0, 14.0, 0.0}},
    {"dsdbmk", {0.17, 0.0, 0.0, 0.0, 3.9, 3.9, 1.0, 14.0, 0.0}},
    {"dsdthcth", {0.39, 0.0, 0.0, 0.0, 4.8, 4.8, 1.0, 14.0, 0.0}},
    {"dsdhcth407", {0.53, 0.0, 0.0, 0.0, 5.0, 5.0, 1.0, 14.0, 0.0}},
    {"dodsvwn5", {0.57, 0.0, 0.0, 0.0, 5.6, 5.6, 1.0, 14.0, 0.0}},
    {"dodblyp", {0.96, 0.0, 0.0, 0.0, 5.1, 5.1, 1.0, 14.0, 0.0}},
    {"dodpbe", {0.91, 0.0, 0.0, 0.0, 5.9, 5.9, 1.0, 14.0, 0.0}},
    {"dodpbep86", {0.72, 0.0, 0.0, 0.0, 5.4, 5.4, 1.0, 14.0, 0.0}},
    {"dodpbeb95", {0.71, 0.0, 0.0, 0.0, 6.0, 6.0, 1.0, 14.0, 0.0}},
    {"dodhsep86", {0.69, 0.0, 0.0, 0.0, 5.4, 5.4, 1.0, 14.0, 0.0}},
    {"dodpbehb95", {0.67, 0.0, 0.0, 0.0, 6.0, 6.0, 1.0, 14.0, 0.0}},
};
// DFT-D3(0)
const std::map<std::string, std::vector<double>> zero = {
    {"slaterdirac", {1.0, 0.999, 0.999, -1.957, 0.697, 0.697, 1.0, 14.0, 0.0}},
    {"bp", {1.0, 1.139, 1.139, 1.683, 1.0, 1.0, 1.0, 14.0, 0.0}},
    {"blyp", {1.0, 1.094, 1.094, 1.682, 1.0, 1.0, 1.0, 14.0, 0.0}},
    {"revpbe", {1.0, 0.923, 0.923, 1.01, 1.0, 1.0, 1.0, 14.0, 0.0}},
    {"rpbe", {1.0, 0.872, 0.872, 0.514, 1.0, 1.0, 1.0, 14.0, 0.0}},
    {"b97_d", {1.0, 0.892, 0.892, 0.909, 1.0, 1.0, 1.0, 14.0, 0.0}},
    {"b973c", {1.0, 1.06, 1.06, 1.5, 1.0, 1.0, 1.0, 14.0, 0.0}},
    {"pbe", {1.0, 1.217, 1.217, 0.722, 1.0, 1.0, 1.0, 14.0, 0.0}},
    {"rpw86pbe", {1.0, 1.224, 1.224, 0.901, 1.0, 1.0, 1.0, 14.0, 0.0}},
    {"b3lyp", {1.0, 1.261, 1.261, 1.703, 1.0, 1.0, 1.0, 14.0, 0.0}},
    {"tpss", {1.0, 1.166, 1.166, 1.105, 1.0, 1.0, 1.0, 14.0, 0.0}},
    {"hf", {1.0, 1.158, 1.158, 1.746, 1.0, 1.0, 1.0, 14.0, 0.0}},
    {"tpss0", {1.0, 1.252, 1.252, 1.242, 1.0, 1.0, 1.0, 14.0, 0.0}},
    {"pbe0", {1.0, 1.287, 1.287, 0.928, 1.0, 1.0, 1.0, 14.0, 0.0}},
    {"hse06", {1.0, 1.129, 1.129, 0.109, 1.0, 1.0, 1.0, 14.0, 0.0}},
    {"hse", {1.0, 1.129, 1.129, 0.109, 1.0, 1.0, 1.0, 14.0, 0.0}}, // ABACUS implements HSE06 as HSE
    {"revpbe38", {1.0, 1.021, 1.021, 0.862, 1.0, 1.0, 1.0, 14.0, 0.0}},
    {"pw6b95", {1.0, 1.532, 1.532, 0.862, 1.0, 1.0, 1.0, 14.0, 0.0}},
    {"b2plyp", {0.64, 1.427, 1.427, 1.022, 1.0, 1.0, 1.0, 14.0, 0.0}},
    {"dsdblyp", {0.5, 1.569, 1.569, 0.705, 1.0, 1.0, 1.0, 14.0, 0.0}},
    {"mpwlyp", {1.0, 1.239, 1.239, 1.098, 1.0, 1.0, 1.0, 14.0, 0.0}},
    {"olyp", {1.0, 0.806, 0.806, 1.764, 1.0, 1.0, 1.0, 14.0, 0.0}},
    {"bpbe", {1.0, 1.087, 1.087, 2.033, 1.0, 1.0, 1.0, 14.0, 0.0}},
    {"opbe", {1.0, 0.837, 0.837, 2.033, 1.0, 1.0, 1.0, 14.0, 0.0}},
    {"ssb", {1.0, 1.215, 1.215, 0.663, 1.0, 1.0, 1.0, 14.0, 0.0}},
    {"revssb", {1.0, 1.221, 1.221, 0.56, 1.0, 1.0, 1.0, 14.0, 0.0}},
    {"otpss", {1.0, 1.128, 1.128, 1.494, 1.0, 1.0, 1.0, 14.0, 0.0}},
    {"b3pw91", {1.0, 1.176, 1.176, 1.775, 1.0, 1.0, 1.0, 14.0, 0.0}},
    {"bhlyp", {1.0, 1.37, 1.37, 1.442, 1.0, 1.0, 1.0, 14.0, 0.0}},
    {"tpssh", {1.0, 1.223, 1.223, 1.219, 1.0, 1.0, 1.0, 14.0, 0.0}},
    {"mpw1b95", {1.0, 1.605, 1.605, 1.118, 1.0, 1.0, 1.0, 14.0, 0.0}},
    {"pwb6k", {1.0, 1.66, 1.66, 0.55, 1.0, 1.0, 1.0, 14.0, 0.0}},
    {"b1b95", {1.0, 1.613, 1.613, 1.868, 1.0, 1.0, 1.0, 14.0, 0.0}},
    {"bmk", {1.0, 1.931, 1.931, 2.168, 1.0, 1.0, 1.0, 14.0, 0.0}},
    {"camb3lyp", {1.0, 1.378, 1.378, 1.217, 1.0, 1.0, 1.0, 14.0, 0.0}},
    {"lcwpbe", {1.0, 1.355, 1.355, 1.279, 1.0, 1.0, 1.0, 14.0, 0.0}},
    {"b2gpplyp", {0.56, 1.586, 1.586, 0.76, 1.0, 1.0, 1.0, 14.0, 0.0}},
    {"ptpss", {0.75, 1.541, 1.541, 0.879, 1.0, 1.0, 1.0, 14.0, 0.0}},
    {"pwpb95", {0.82, 1.557, 1.557, 0.705, 1.0, 1.0, 1.0, 14.0, 0.0}},
    {"pw1pw", {1.0, 1.4968, 1.4968, 1.1786, 1.0, 1.0, 1.0, 14.0, 0.0}},
    {"scan", {1.0, 1.324, 1.324, 0.0, 1.0, 1.0, 1.0, 14.0, 0.0}},
    {"wb97x_d3", {1.0, 1.281, 1.281, 1.0, 1.094, 1.094, 1.0, 14.0, 0.0}}, 
    // NOTE: simple-dftd3 assign the D3(0) parameters of functional wB97X-D3
    // to a key `wb97x`, but the functional wB97X itself does not own these params.
    // instead, there is a XC in libxc really names HYB_GGA_WB97X_D3
    {"pbehpbe", {1.0, 1.5703, 1.5703, 1.401, 1.0, 1.0, 1.0, 14.0, 0.0}},
    {"xlyp", {1.0, 0.9384, 0.9384, 0.7447, 1.0, 1.0, 1.0, 14.0, 0.0}},
    {"mpwpw", {1.0, 1.3725, 1.3725, 1.9467, 1.0, 1.0, 1.0, 14.0, 0.0}},
    {"hcth407", {1.0, 4.0426, 4.0426, 2.7694, 1.0, 1.0, 1.0, 14.0, 0.0}},
    {"revtpss", {1.0, 1.3491, 1.3491, 1.3666, 1.0, 1.0, 1.0, 14.0, 0.0}},
    {"tauhcth", {1.0, 0.932, 0.932, 0.5662, 1.0, 1.0, 1.0, 14.0, 0.0}},
    {"b3p", {1.0, 1.1897, 1.1897, 1.1961, 1.0, 1.0, 1.0, 14.0, 0.0}},
    {"b1p", {1.0, 1.1815, 1.1815, 1.1209, 1.0, 1.0, 1.0, 14.0, 0.0}},
    {"b1lyp", {1.0, 1.3725, 1.3725, 1.9467, 1.0, 1.0, 1.0, 14.0, 0.0}},
    {"mpwb1k", {1.0, 1.671, 1.671, 1.061, 1.0, 1.0, 1.0, 14.0, 0.0}},
    {"mpw1lyp", {1.0, 2.0512, 2.0512, 1.9529, 1.0, 1.0, 1.0, 14.0, 0.0}},
    {"mpw1pw", {1.0, 1.2892, 1.2892, 1.4758, 1.0, 1.0, 1.0, 14.0, 0.0}},
    {"mpw1kcis", {1.0, 1.7231, 1.7231, 2.2917, 1.0, 1.0, 1.0, 14.0, 0.0}},
    {"mpwkcis1k", {1.0, 1.4853, 1.4853, 1.7553, 1.0, 1.0, 1.0, 14.0, 0.0}},
    {"pbeh1pbe", {1.0, 1.3719, 1.3719, 1.043, 1.0, 1.0, 1.0, 14.0, 0.0}},
    {"pbe1kcis", {1.0, 3.6355, 3.6355, 1.7934, 1.0, 1.0, 1.0, 14.0, 0.0}},
    {"x3lyp", {1.0, 1.0, 1.0, 0.299, 1.0, 1.0, 1.0, 14.0, 0.0}},
    {"o3lyp", {1.0, 1.406, 1.406, 1.8058, 1.0, 1.0, 1.0, 14.0, 0.0}},
    {"b97_1", {1.0, 3.7924, 3.7924, 1.6418, 1.0, 1.0, 1.0, 14.0, 0.0}},
    {"b97_2", {1.0, 1.7066, 1.7066, 1.6418, 1.0, 1.0, 1.0, 14.0, 0.0}},
    {"b98", {1.0, 2.6895, 2.6895, 1.9078, 1.0, 1.0, 1.0, 14.0, 0.0}},
    {"hiss", {1.0, 1.3338, 1.3338, 0.7615, 1.0, 1.0, 1.0, 14.0, 0.0}},
    {"hse03", {1.0, 1.3944, 1.3944, 1.0156, 1.0, 1.0, 1.0, 14.0, 0.0}},
    {"revtpssh", {1.0, 1.3224, 1.3224, 1.2504, 1.0, 1.0, 1.0, 14.0, 0.0}},
    {"revtpss0", {1.0, 1.2881, 1.2881, 1.0649, 1.0, 1.0, 1.0, 14.0, 0.0}},
    {"tpss1kcis", {1.0, 1.7729, 1.7729, 2.0902, 1.0, 1.0, 1.0, 14.0, 0.0}},
    {"tauhcthhyb", {1.0, 1.5001, 1.5001, 1.6302, 1.0, 1.0, 1.0, 14.0, 0.0}},
    {"pkzb", {1.0, 0.6327, 0.6327, 0.0, 1.0, 1.0, 1.0, 14.0, 0.0}},
    {"n12", {1.0, 1.3493, 1.3493, 2.3916, 1.0, 1.0, 1.0, 14.0, 0.0}},
    {"mpw2plyp", {0.66, 1.5527, 1.5527, 0.7529, 1.0, 1.0, 1.0, 14.0, 0.0}},
    {"m05", {1.0, 1.373, 1.373, 0.595, 1.0, 1.0, 1.0, 14.0, 0.0}},
    {"m052x", {1.0, 1.417, 1.417, 0.0, 1.0, 1.0, 1.0, 14.0, 0.0}},
    {"m06l", {1.0, 1.581, 1.581, 0.0, 1.0, 1.0, 1.0, 14.0, 0.0}},
    {"m06", {1.0, 1.325, 1.325, 0.0, 1.0, 1.0, 1.0, 14.0, 0.0}},
    {"m062x", {1.0, 1.619, 1.619, 0.0, 1.0, 1.0, 1.0, 14.0, 0.0}},
    {"m08hx", {1.0, 1.6247, 1.6247, 0.0, 1.0, 1.0, 1.0, 14.0, 0.0}},
    {"m11l", {1.0, 2.3933, 2.3933, 1.1129, 1.0, 1.0, 1.0, 14.0, 0.0}},
    {"mn15l", {1.0, 3.3388, 3.3388, 0.0, 1.0, 1.0, 1.0, 14.0, 0.0}},
    {"pwp", {1.0, 2.104, 2.104, 0.8747, 1.0, 1.0, 1.0, 14.0, 0.0}},
};
// DFT-D3M(BJ): not implemented for beta
const std::map<std::string, std::vector<double>> bjm = {
    {"bp", {1.0, 0.82185, 0.82185, 3.140281, 2.728151, 2.728151, 1.0, 14.0, 0.0}},
    {"blyp", {1.0, 0.448486, 0.448486, 1.875007, 3.610679, 3.610679, 1.0, 14.0, 0.0}},
    {"b97_d", {1.0, 0.240184, 0.240184, 1.206988, 3.864426, 3.864426, 1.0, 14.0, 0.0}},
    {"pbe", {1.0, 0.012092, 0.012092, 0.35894, 5.938951, 5.938951, 1.0, 14.0, 0.0}},
    {"b3lyp", {1.0, 0.278672, 0.278672, 1.466677, 4.606311, 4.606311, 1.0, 14.0, 0.0}},
    {"pbe0", {1.0, 0.007912, 0.007912, 0.528823, 6.162326, 6.162326, 1.0, 14.0, 0.0}},
    {"b2plyp", {0.64, 0.486434, 0.486434, 0.67282, 3.656466, 3.656466, 1.0, 14.0, 0.0}},
    {"lcwpbe", {1.0, 0.563761, 0.563761, 0.906564, 3.59368, 3.59368, 1.0, 14.0, 0.0}},
};
// DFT-D3M(0): not implemented for beta
const std::map<std::string, std::vector<double>> zerom = {
    {"bp", {1.0, 1.23346, 1.23346, 1.945174, 1.0, 1.0, 1.0, 14.0, 0.0}},
    {"blyp", {1.0, 1.279637, 1.279637, 1.841686, 1.0, 1.0, 1.0, 14.0, 0.01437}},
    {"b97_d", {1.0, 1.151808, 1.151808, 1.020078, 1.0, 1.0, 1.0, 14.0, 0.035964}},
    {"pbe", {1.0, 2.340218, 2.340218, 0.0, 1.0, 1.0, 1.0, 14.0, 0.129434}},
    {"b3lyp", {1.0, 1.338153, 1.338153, 1.532981, 1.0, 1.0, 1.0, 14.0, 0.013988}},
    {"pbe0", {1.0, 2.077949, 2.077949, 8.1e-05, 1.0, 1.0, 1.0, 14.0, 0.116755}},
    {"b2plyp", {0.64, 1.313134, 1.313134, 0.717543, 1.0, 1.0, 1.0, 14.0, 0.016035}},
    {"lcwpbe", {1.0, 1.366361, 1.366361, 1.280619, 1.0, 1.0, 1.0, 14.0, 0.00316}},
};
// DFT-D3(OptimizedPower)
const std::map<std::string, std::vector<double>> op = {
    {"blyp", {1.0, 0.425, 0.425, 1.31867, 3.5, 3.5, 1.0, 14.0, 2.0}},
    {"revpbe", {1.0, 0.6, 0.6, 1.44765, 2.5, 2.5, 1.0, 14.0, 0.0}},
    {"b97_d", {1.0, 0.6, 0.6, 1.46861, 2.5, 2.5, 1.0, 14.0, 0.0}},
    {"pbe", {0.91826, 0.2, 0.2, 0.0, 4.75, 4.75, 1.0, 14.0, 6.0}},
    {"b3lyp", {1.0, 0.3, 0.3, 0.78311, 4.25, 4.25, 1.0, 14.0, 4.0}},
    {"tpss", {1.0, 0.575, 0.575, 0.51581, 3.0, 3.0, 1.0, 14.0, 8.0}},
    {"pbe0", {0.8829, 0.15, 0.15, 0.0, 4.75, 4.75, 1.0, 14.0, 6.0}},
    {"revpbe0", {1.0, 0.725, 0.725, 1.25684, 2.25, 2.25, 1.0, 14.0, 0.0}},
    {"tpssh", {1.0, 0.575, 0.575, 0.43185, 3.0, 3.0, 1.0, 14.0, 8.0}},
    {"revtpss", {1.0, 0.7, 0.7, 0.27632, 2.5, 2.5, 1.0, 14.0, 8.0}},
    {"b97_1", {0.97388, 0.15, 0.15, 0.0, 4.25, 4.25, 1.0, 14.0, 6.0}},
    {"revtpssh", {1.0, 0.575, 0.575, 0.12467, 3.0, 3.0, 1.0, 14.0, 10.0}},
    {"ms2", {1.0, 0.7, 0.7, 0.90743, 4.0, 4.0, 1.0, 14.0, 2.0}},
    {"ms2h", {1.0, 0.65, 0.65, 1.69464, 4.75, 4.75, 1.0, 14.0, 0.0}},
};
    

// 's6', 'rs6', 'a1', 's8', 'rs8', 'a2', 's9', 'alp', 'bet'
/**
 * @brief Get the dftd3 params object. 
 * dftd3 method fall back: xc-bjm -> xc-bj -> pbe-bj
 *                         xc-zerom -> xc-zero -> pbe-zero
 * 
 * @param xc the functional name
 * @param d3method the d3 method, can be "bj", "zero-damping", "bj-modified", "zero-damping-modified", "op"
 * @param param the dftd3 parameters, ALL_KEYS = {'s6', 'rs6', 'a1', 's8', 'rs8', 'a2', 's9', 'alp', 'bet'}
 */
void _search(const std::string& xc, 
                const std::string& method, 
                std::vector<double>& param)
{
    const std::string xc_lowercase = FmtCore::lower(xc);
    const std::vector<std::string> allowed_ = { "bj", "zero", "bjm", "zerom", "op" };
    assert(std::find(allowed_.begin(), allowed_.end(), method) != allowed_.end());
    if (method == "op")
    {
        if (op.find(xc_lowercase) != op.end())
        {
            param = op.at(xc_lowercase);
        }
        else
        {
            ModuleBase::WARNING_QUIT("ModuleHamiltGeneral::ModuleVDW::DFTD3::_search",
                                        "XC (`" + xc + "`)'s DFT-D3(OP) parameters not found");
        }
    }
    else if (method == "bjm")
    {
        if (bjm.find(xc_lowercase) != bjm.end())
        {
            param = bjm.at(xc_lowercase);
        }
        else
        {
            ModuleBase::WARNING_QUIT("ModuleHamiltGeneral::ModuleVDW::DFTD3::_search",
                                        "XC (`" + xc + "`)'s DFT-D3M(BJ) parameters not found");
        }
    }
    else if (method == "bj")
    {
        if (bj.find(xc_lowercase) != bj.end())
        {
            param = bj.at(xc_lowercase);
        }
        else
        {
            ModuleBase::WARNING_QUIT("ModuleHamiltGeneral::ModuleVDW::DFTD3::_search",
                                        "XC (`" + xc + "`)'s DFT-D3(BJ) parameters not found");
        }
    }
    else if (method == "zerom")
    {
        if (zerom.find(xc_lowercase) != zerom.end())
        {
            param = zerom.at(xc_lowercase);
        }
        else
        {
            ModuleBase::WARNING_QUIT("ModuleHamiltGeneral::ModuleVDW::DFTD3::_search",
                                        "XC (`" + xc + "`)'s DFT-D3M(0) parameters not found");
        }
    }
    else if (method == "zero")
    {
        if (zero.find(xc_lowercase) != zero.end())
        {
            param = zero.at(xc_lowercase);
        }
        else
        {
            ModuleBase::WARNING_QUIT("ModuleHamiltGeneral::ModuleVDW::DFTD3::_search",
                                        "XC (`" + xc + "`)'s DFT-D3(0) parameters not found");
        }
    }
    else // should not reach here
    {
        ModuleBase::WARNING_QUIT("ModuleHamiltGeneral::ModuleVDW::DFTD3::_search",
                                    "Unknown DFT-D3 method: " + method);
    }
}

/**
 * @brief Get DFT-D3 parameters. If if there are parameters defined,
 * then it will overwrite the search result. If all parameters are
 * defined already by user, then search will not performed. 
 * 
 * @param xc XC functional name
 * @param d3method can be "d3_0" or "d3_bj"
 * @param s6_in user defined s6, default is "default"
 * @param s8_in user defined s8, default is "default"
 * @param a1_in user defined a1, default is "default"
 * @param a2_in user defined a2, default is "default"
 * @param s6 [out] s6 parameter
 * @param s8 [out] s8 parameter
 * @param a1 [out] a1 parameter
 * @param a2 [out] a2 parameter
 */
void vdw::Vdwd3Parameters::_vdwd3_autoset_xcparam(const std::string& xc_in,
                                                  const std::string& d3method,
                                                  const std::string& s6_in,
                                                  const std::string& s8_in,
                                                  const std::string& a1_in,
                                                  const std::string& a2_in,
                                                  double& s6,
                                                  double& s8,
                                                  double& a1,
                                                  double& a2,
                                                  std::ofstream* plog)
{
    const std::map<std::string, std::string> param_map = {
        {"d3_bj", "bj"}, {"d3_0", "zero"}, {"d3_bjm", "bjm"}, {"d3_0m", "zerom"},
        {"op", "op"}};

    const std::vector<std::string> flag = {s6_in, s8_in, a1_in, a2_in};
    const bool autoset = std::any_of(flag.begin(), flag.end(), [](const std::string& s) { return s == "default"; });
    if (!autoset) // all parameters are defined
    {
        s6 = std::stod(s6_in);
        s8 = std::stod(s8_in);
        a1 = std::stod(a1_in);
        a2 = std::stod(a2_in);
    }
    else
    {
        std::vector<double> param;
        const std::string xc = _vdwd3_xcname(xc_in);
        _search(xc, param_map.at(d3method), param);
        s6 = (s6_in == "default") ? param[0] : std::stod(s6_in);
        s8 = (s8_in == "default") ? param[3] : std::stod(s8_in);
        a1 = (a1_in == "default") ? param[2] : std::stod(a1_in);
        a2 = (a2_in == "default") ? param[5] : std::stod(a2_in);
        if (plog != nullptr) // logging the autoset
        {
            param = {s6, s8, a1, a2};
            FmtTable vdwd3tab({"Parameters", "Original", "Autoset"}, 4, {"%10s", "%10s", "%10.4f"});
            const std::vector<std::string> items = {"s6", "s8", "a1", "a2"};
            vdwd3tab << items << flag << param;
            (*plog) << "\nDFT-D3 Dispersion correction parameters autoset\n" << vdwd3tab.str()
                    << "XC functional: " << xc_in << std::endl;
        }

    }
}


/*
'''
dftd3 parameters from 
https://github.com/dftd3/simple-dftd3/blob/main/assets/parameters.toml

this script is to convert the toml file to c++ map
'''

import toml

def load(fn):
    with open(fn, 'r') as f:
        data = toml.load(f)
    return data

def xc_indexing(data):
    out = {'bj': {}, 'zero': {}, 'bjm': {}, 'zerom': {}, 'op': {}}
    for xc, param in data['parameter'].items():
        for vdw_method, value in param['d3'].items():
            out[vdw_method][xc] = {k: v for k, v in value.items() if k != 'doi'}
    return out

def complete(vdw_method, value):
    '''
    for each functional, the zero damping version must be provided
    for each vdw method, all parameters including 
    s6, rs6/a1, s8, rs8/a2, s9, alp, bet must be provided, otherwise
    use the default value 
    '''
    DEFAULT = {
        'bj': {'s6': 1.0, 's9': 1.0, 'alp': 14.0},
        'zero': {'s6': 1.0, 's9': 1.0, 'rs8': 1.0, 'alp': 14.0},
        'bjm': {'s6': 1.0, 's9': 1.0, 'alp': 14.0},
        'zerom': {'s6': 1.0, 's9': 1.0, 'rs8': 1.0, 'alp': 14.0},
        'op': {'s9': 1.0, 'alp': 14.0}
    }
    ALL_KEYS = {'s6', 'rs6', 'a1', 's8', 'rs8', 'a2', 's9', 'alp', 'bet'}
    EQUIVALENT = {'rs6': 'a1', 'a1': 'rs6', 'rs8': 'a2', 'a2': 'rs8'}
    out = value.copy()
    for k in ALL_KEYS:
        equilk = EQUIVALENT.get(k, k)
        val = [out.get(k), out.get(equilk), 
               DEFAULT[vdw_method].get(k), DEFAULT[vdw_method].get(equilk)]
        val = [v for v in val if v is not None]
        val = [0.0] if not val else val
        out[k] = val[0]
        out[equilk] = out[k]
    # equivalent? 
    # according to 
    # abacus-develop/source/module_hamilt_general/module_vdw/vdwd3_parameters.cpp
    # https://abacus.deepmodeling.com/en/latest/advanced/input_files/input-main.html

    return out

def make_stdmap(data):
    for vdw_method, param in data.items():
        print(f'std::map<std::string, std::vector<double>> {vdw_method} = {{')
        for xc, value in param.items():
            print(f'    {{\"{xc}\", {{', end='')
            print(', '.join([f'{v}' for v in value.values()]), end='')
            print('}},')
        print('};')
    
if __name__ == '__main__':
    fn = 'dftd3.toml'
    data = load(fn)
    data = xc_indexing(data)
    for vdw_method, param in data.items():
        for xc, value in param.items():
            raw = complete(vdw_method, value)
            data[vdw_method][xc] = {k: raw[k] 
            for k in ['s6', 'rs6', 'a1', 's8', 'rs8', 'a2', 's9', 'alp', 'bet']}
    make_stdmap(data)
*/

