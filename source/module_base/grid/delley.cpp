#include "module_base/grid/delley.h"

/**
 * @brief Delley's table for quadrature on the unit sphere.
 *
 * Reference:
 * Delley, B. (1996). High order integration schemes on the unit sphere.
 * Journal of computational chemistry, 17(9), 1152-1155.
 *
 */
struct Table {
    const int lmax_;
    const int ngrid_;
    const int ntype_[6];
    const std::vector<double> data_;
};

const std::vector<Table> table = {
    {
        17, 110, {1, 1, 0, 3, 1, 0},
        {
            0.00000000000000000, 0.00000000000000000, 0.0038282704949371616,
            0.97735026918962576, 0.57735026918962576, 0.0097937375124875125,
            0.18511563534473617, 0.18511563534473617, 0.0082117372831911110,
            0.39568947305594191, 0.39568947305594191, 0.0095954713360709628,
            0.69042104838229218, 0.21595729184584883, 0.0099428148911781033,
            0.47836902881215020, 0.00000000000000000, 0.0096949963616630283,
        }
    },
    {
        23, 194, {1, 1, 1, 4, 1, 1},
        {
            0.00000000000000000, 0.00000000000000000, 0.0017823404472446112,
            0.57735026918962576, 0.57735026918962576, 0.0055733831788487380,
            0.70710678118654752, 0.00000000000000000, 0.0057169059499771019,
            0.44469331787174373, 0.44469331787174373, 0.0055187714672736137,
            0.28924656275754386, 0.28924656275754386, 0.0051582377118053831,
            0.67129734426952263, 0.31419699418258608, 0.0056087040825879968,
            0.12993354476500669, 0.12993354476500669, 0.0041067770281693941,
            0.34577021976112827, 0.00000000000000000, 0.0050518460646148085,
            0.52511857244364202, 0.15904171053835295, 0.0055302489162330937,
        }
    },
    {
        29, 302, {1, 1, 0, 6, 2, 2},
        {
            0.00000000000000000, 0.00000000000000000, 0.0008545911725128148,
            0.57735026918962576, 0.57735026918962576, 0.0035991192850255715,
            0.70117664160895449, 0.12923867271051493, 0.0036500458076772554,
            0.65663294102196118, 0.37103417838482119, 0.0036048226014198817,
            0.47290541325810046, 0.47290541325810046, 0.0035767296617433671,
            0.35156403455701051, 0.35156403455701051, 0.0034497884243058833,
            0.22196452362941784, 0.22196452362941784, 0.0031089531224136753,
            0.09618308522614784, 0.09618308522614784, 0.0023521014136891644,
            0.97189558918789607, 0.00000000000000000, 0.0036008209322164603,
            0.26441528870606625, 0.00000000000000000, 0.0029823449631718039,
            0.04486773725807738, 0.25100347517704651, 0.0035715405542733871,
            0.41277240831685310, 0.12335485325833274, 0.0033923122050061702,
        }
    },
    {
        35, 434, {1, 1, 1, 7, 2, 4},
        {
            0.00000000000000000, 0.00000000000000000, 0.0005265897968224436, 
            0.57735026918962576, 0.57735026918962576, 0.0025123174189273072,
            0.70710678118654752, 0.00000000000600000, 0.0025482199720026072,
            0.69093463075091106, 0.21264682470755207, 0.0025304038011863550,
            0.64566647074242561, 0.40771266489776951, 0.0025132671745975644,
            0.49143426377847465, 0.49143426377847465, 0.0025017251684029361,
            0.39272597633680022, 0.39272597633680022, 0.0024453734373129800,
            0.28612890103076384, 0.28612890103076384, 0.0023026947822274158,
            0.17748360546091578, 0.17748360546091578, 0.0020142790209185282,
            0.07568084367178018, 0.07568084367178018, 0.0014624956215946138,
            0.21027252285730696, 0.00000000000000000, 0.0019109512821795323,
            0.47159869115131592, 0.00000000000000000, 0.0024174423756389808,
            0.33443631453434549, 0.09921769636429237, 0.0022366077604378487,
            0.45023303825826254, 0.20548236964030437, 0.0024169300443247753,
            0.55501523610768072, 0.31042840351665415, 0.0024966440545530860,
            0.99051570489252711, 0.10680182607580483, 0.0025122368545634951,
        }
    },
    {
        41, 590, {1, 1, 0, 8, 4, 6},
        {
            0.00000000000000000, 0.00000000000000000, 0.0601009005753378758, 
            0.57735026918962576, 0.57735026918962576, 0.0018514016873890461,
            0.70404760433146996, 0.09291900596883211, 0.0018686219518306975,
            0.68084561988024238, 0.26999719217017240, 0.0018648696345606001,
            0.63723669159418917, 0.43342680786054810, 0.0018497643975168892,
            0.50447558060926046, 0.50447558060926046, 0.0018450277740822388,
            0.42175447334398773, 0.42175447334398773, 0.0018164174988262214,
            0.33201962086729379, 0.33201962086729379, 0.0017449464690023229,
            0.23917494336556047, 0.23917494336556047, 0.0016278016126848035,
            0.14024070738935403, 0.14024070738935403, 0.0015576827519901693,
            0.09161634328605240, 0.00000000000000000, 0.0012680968886048433,
            0.20326292518419433, 0.00000000000000000, 0.0011183965414769017,
            0.39364042372978295, 0.00000000000000000, 0.0017287035120530033,
            0.61262355812929648, 0.00000000000000000, 0.0018551905629473527,
            0.28114771623428322, 0.08959875911893791, 0.0014697353123693616,
            0.38175470908581117, 0.17327600238498666, 0.0016819651914742022,
            0.47452376478986998, 0.26422260656245780, 0.0017876372876796954,
            0.56127905075920534, 0.35189965873835832, 0.0018400735685528423,
            0.50324791996964975, 0.08886791018186295, 0.0018072536817113700,
            0.59768324320748616, 0.18154345643517542, 0.0018527289739424312,
        }
    },
    {
        47, 770, {1, 1, 1, 9, 4, 9},
        {
            0.00000000000000000, 0.00000000000000000, 0.0011685335608691628, 
            0.57735026918962576, 0.57735026918962576, 0.0014121215930643264,
            0.70710678118654752, 0.00000000000000000, 0.0014468645950992776,
            0.11441365123336336, 0.11441365123336336, 0.0010478418864629224,
            0.19944675708548970, 0.19944675708548970, 0.0012392547584848484,
            0.28401278368259530, 0.28401278368259530, 0.0013259295792415379,
            0.36646411416548296, 0.36646411416548296, 0.0013756097758625958,
            0.44356118052513995, 0.44356118052513995, 0.0013999348863558624,
            0.51435709575333968, 0.51435709575333968, 0.0014096221218822673,
            0.63052081196671812, 0.45264446462279973, 0.0014108746499638577,
            0.67164784337293865, 0.31269529735024947, 0.0014134887639034478,
            0.69812332010174177, 0.15889512220405632, 0.0014366946685816802,
            0.12047667931264991, 0.00000000000000000, 0.0010901543574180667,
            0.30940302315480606, 0.00000000000000000, 0.0001869137844803852,
            0.34884276430183016, 0.00000000000000000, 0.0011284267652336505,
            0.03224214285417946, 0.00000000000000000, 0.0013844558026568455,
            0.23249923409267532, 0.06616159933437003, 0.0011853923885095502,
            0.32477344409682044, 0.14568618765136356, 0.0012949021664637693,
            0.41056989039349425, 0.22832839132127622, 0.0013525857420363760,
            0.49213658085114203, 0.30714431901543855, 0.0013925025908786082,
            0.56548849812588755, 0.38271180625074657, 0.0014073257894372725,
            0.43713473693946563, 0.07970715187939190, 0.0013128954307755017,
            0.52320749473197761, 0.15892620239864833, 0.0013784632898490457,
            0.60283033994386521, 0.23667220253873893, 0.0014125450609821936,
            0.62037164721742807, 0.07982328826030880, 0.0014289835314095131,
        }
    },
    {
        53, 974, {1, 1, 0, 12, 4, 12},
        {
            0.00000000000000000, 0.00000000800000000, 0.0001438294190527431, 
            0.57735026918962576, 0.57735026918962576, 0.0011257722882870041,
            0.04292963545341347, 0.04292963545341347, 0.0004948029341949241,
            0.10514268540864042, 0.10514268540864042, 0.0007357990109125470,
            0.17500248676230874, 0.17500248676230874, 0.0008889132771304384,
            0.24776533796502568, 0.24776533796502568, 0.0009888347838921435,
            0.32065671239559574, 0.32065671239559574, 0.0010532996817094706,
            0.39165207498499835, 0.39165207498499835, 0.0010927788070145785,
            0.45908258741876237, 0.45908258741876237, 0.0011143893940632272,
            0.52145638884158605, 0.52145638884158605, 0.0011237247880515553,
            0.62531702446541989, 0.46685890569574328, 0.0011252393252438136,
            0.66379267445231699, 0.34461365423743822, 0.0011261532718159050,
            0.69104103984983007, 0.21195415185018465, 0.0011302869311238468,
            0.70529070074577603, 0.07162440144995566, 0.0011349865343639549,
            0.12366867626579899, 0.00000000000000000, 0.0006823367927109931,
            0.29407771144683870, 0.00000000000000000, 0.0009454158160447096,
            0.46977538492076491, 0.00000000000000000, 0.0010744299753856791,
            0.63345632411395669, 0.00000000000000000, 0.0011293000865691317,
            0.20291287527775228, 0.05974048614181342, 0.0008436884500901954,
            0.46026219424840539, 0.13757604084736365, 0.0010752557204488846,
            0.50306739996620357, 0.33910165263362857, 0.0011085772368644620,
            0.28176064224421343, 0.12716751914398195, 0.0009566475323783357,
            0.43315612917201574, 0.26931207404135125, 0.0010806632507173907,
            0.62561673585808142, 0.14197864526019183, 0.0011267971311962946,
            0.37983952168591567, 0.06709284600738255, 0.0010225687153580612,
            0.55175054214235205, 0.07057738183256172, 0.0011089602677131075,
            0.60296191561591869, 0.27838884778821546, 0.0011227906534357658,
            0.35896063295890958, 0.19795789389174069, 0.0010324018471174598,
            0.53486664381354765, 0.20873070611032740, 0.0011072493822838539,
            0.96749975460743735, 0.40551221378728359, 0.0011217800485199721,
        }
    },
    {
        59, 1202, {1, 1, 1, 13, 4, 16},
        {
            0.00006000000000000, 0.00000000000000000, 0.0001105189233267572, 
            0.97735026918962576, 0.57735026918962576, 0.0009133159786443561,
            0.70710678118654752, 0.00000000000000000, 0.0009205232738090741,
            0.03712636449657089, 0.03712636449657089, 0.0003690421898017899,
            0.09140060412262223, 0.09140060412262223, 0.0005603990928680660,
            0.15310778524699062, 0.15310778524699062, 0.0006865297629282609,
            0.21809288916606116, 0.21809288916606116, 0.0007720338551145630,
            0.28398745322001746, 0.28398745322001746, 0.0008301545958894795,
            0.34911776009637644, 0.34911776009637644, 0.0008686692550179628,
            0.41214314614443092, 0.41214314614443092, 0.0008927076285846890,
            0.47189936271491266, 0.47189936271491266, 0.0009060820238568219,
            0.52731454528423366, 0.52731454528423366, 0.0009119777254940867,
            0.62094753324440192, 0.47838093807695216, 0.0009128720138604181,
            0.65697227118572905, 0.36983086645942597, 0.0009130714935691735,
            0.68417883090701434, 0.25258395570071777, 0.0009152873784554116,
            0.70126043301236308, 0.12832618665972300, 0.0009187436274321654,
            0.10723822154781661, 0.00000000000000000, 0.0005176977312965694,
            0.29820689594969680, 0.00006000000000000, 0.0007331143682101417,
            0.41727529553067168, 0.00000000000000000, 0.0008463232836379928,
            0.57003669117925033, 0.00000000000000000, 0.0009031122694253992,
            0.17717740226153253, 0.05210639477011284, 0.0006485778453163257,
            0.24757164634262876, 0.11156409571564867, 0.0007435030910982369,
            0.31736152466119767, 0.17465516775786261, 0.0008101731497468018,
            0.38542911506692237, 0.23902784793817240, 0.0008556299257311812,
            0.45074225931570644, 0.30294669735289819, 0.0008850282341265444,
            0.51235184864198708, 0.36498322605976536, 0.0009022692938426915,
            0.56937024984684411, 0.42386447815223403, 0.0009105760258970126,
            0.33546162890664885, 0.05905888853235508, 0.0007998527891839054,
            0.40902684270853572, 0.12172350510959870, 0.0008483389574594331,
            0.47853206759224352, 0.18575051945473351, 0.0008811048182425720,
            0.54343035696939004, 0.24941121623622365, 0.0009010091677105086,
            0.60311616930963100, 0.31122759471496082, 0.0009107813579482705,
            0.49322211848512846, 0.06266250624154169, 0.0008803208679738260,
            0.96321230207620997, 0.12677748006842827, 0.0009021342299040653,
            0.62698055090243917, 0.19060182227792370, 0.0009131578003189435,
            0.63942796347491023, 0.06424549224220589, 0.0009158016174693465,
        }
    }
};

#include <cstdio>
int main() {

    const double x = 0.000000012345678987654321;
    const double w = 3.8282704949371616e-03;
    printf("%25.17e\n", x);
    printf("%25.17e\n", w);
    printf("%25.17e\n", table[0].data_[2]);
    printf("%25.17e\n", table[1].data_[2]);
    printf("%25.17e\n", table.back().data_[2]);

}
