#include "source_base/constants.h"
#include "source_base/tool_quit.h"
#include "read_input.h"
#include "read_input_tool.h"

#include <algorithm>
#include <array>

namespace ModuleIO
{
namespace
{

std::vector<int> parse_supersine_steps(const Input_Item& item, const int default_step)
{
    const std::vector<std::string> tokens
        = item.is_read() ? item.str_values : std::vector<std::string>(1, "default");
    std::vector<int> values;
    for (const std::string& token: tokens)
    {
        if (token == "default")
        {
            values.push_back(default_step);
        }
        else
        {
            std::vector<int> parsed;
            parse_expression(std::vector<std::string>(1, token), parsed);
            values.insert(values.end(), parsed.begin(), parsed.end());
        }
    }
    return values;
}

struct FieldParameterRule
{
    const char* name;
    int field_type;
    std::size_t size;
};

} // namespace

void check_td_efield_parameters(const Input_para& input)
{
    if (input.td_ttype.size() != input.td_vext_dire.size())
    {
        ModuleBase::WARNING_QUIT("ReadInput",
                                 "td_ttype and td_vext_dire must contain the same number of fields.");
    }
    if (input.td_stype < 0 || input.td_stype > 2)
    {
        ModuleBase::WARNING_QUIT("ReadInput", "td_stype must be 0, 1, or 2.");
    }
    for (const int direction: input.td_vext_dire)
    {
        if (direction < 1 || direction > 3)
        {
            ModuleBase::WARNING_QUIT("ReadInput", "td_vext_dire must contain values from 1 to 3.");
        }
    }

    std::array<std::size_t, 5> field_counts = {{0, 0, 0, 0, 0}};
    for (const int field_type: input.td_ttype)
    {
        if (field_type < 0 || field_type > 4)
        {
            ModuleBase::WARNING_QUIT("ReadInput", "td_ttype must contain values from 0 to 4.");
        }
        ++field_counts[field_type];
    }

    const std::array<FieldParameterRule, 24> rules = {{
        {"td_gauss_freq", 0, input.td_gauss_freq.size()},
        {"td_gauss_phase", 0, input.td_gauss_phase.size()},
        {"td_gauss_sigma", 0, input.td_gauss_sigma.size()},
        {"td_gauss_t0", 0, input.td_gauss_t0.size()},
        {"td_gauss_amp", 0, input.td_gauss_amp.size()},
        {"td_trape_freq", 1, input.td_trape_freq.size()},
        {"td_trape_phase", 1, input.td_trape_phase.size()},
        {"td_trape_t1", 1, input.td_trape_t1.size()},
        {"td_trape_t2", 1, input.td_trape_t2.size()},
        {"td_trape_t3", 1, input.td_trape_t3.size()},
        {"td_trape_amp", 1, input.td_trape_amp.size()},
        {"td_trigo_freq1", 2, input.td_trigo_freq1.size()},
        {"td_trigo_freq2", 2, input.td_trigo_freq2.size()},
        {"td_trigo_phase1", 2, input.td_trigo_phase1.size()},
        {"td_trigo_phase2", 2, input.td_trigo_phase2.size()},
        {"td_trigo_amp", 2, input.td_trigo_amp.size()},
        {"td_heavi_t0", 3, input.td_heavi_t0.size()},
        {"td_heavi_amp", 3, input.td_heavi_amp.size()},
        {"td_supsine_amp", 4, input.td_supsine_amp.size()},
        {"td_supsine_freq", 4, input.td_supsine_freq.size()},
        {"td_supsine_phase", 4, input.td_supsine_phase.size()},
        {"td_supsine_sigma", 4, input.td_supsine_sigma.size()},
        {"td_supsine_tstart", 4, input.td_supsine_tstart.size()},
        {"td_supsine_tend", 4, input.td_supsine_tend.size()},
    }};
    for (const FieldParameterRule& rule: rules)
    {
        const std::size_t field_count = field_counts[rule.field_type];
        if (field_count > 0 && rule.size != field_count)
        {
            ModuleBase::WARNING_QUIT("ReadInput",
                                     std::string(rule.name)
                                         + " must contain exactly one value for each td_ttype "
                                         + std::to_string(rule.field_type) + " occurrence.");
        }
    }

    for (std::size_t index = 0; index < field_counts[0]; ++index)
    {
        if (input.td_gauss_sigma[index] == 0.0)
        {
            ModuleBase::WARNING_QUIT("ReadInput", "td_gauss_sigma must be nonzero.");
        }
    }
    for (std::size_t index = 0; index < field_counts[1]; ++index)
    {
        if (input.td_trape_t1[index] > input.td_trape_t2[index]
            || input.td_trape_t2[index] > input.td_trape_t3[index])
        {
            ModuleBase::WARNING_QUIT(
                "ReadInput",
                "Each trapezoid field must satisfy td_trape_t1 <= td_trape_t2 <= td_trape_t3.");
        }
    }
    for (std::size_t index = 0; index < field_counts[4]; ++index)
    {
        if (input.td_supsine_freq[index] == 0.0)
        {
            ModuleBase::WARNING_QUIT("ReadInput", "td_supsine_freq must be nonzero.");
        }
        if (input.td_supsine_sigma[index] <= 0.0
            || input.td_supsine_sigma[index] >= ModuleBase::PI / 2.0)
        {
            ModuleBase::WARNING_QUIT("ReadInput",
                                     "td_supsine_sigma must be greater than 0 and less than pi/2.");
        }
        if (input.td_supsine_tstart[index] < input.td_tstart
            || input.td_supsine_tstart[index] >= input.td_supsine_tend[index]
            || input.td_supsine_tend[index] > input.td_tend)
        {
            ModuleBase::WARNING_QUIT(
                "ReadInput",
                "Each supersine pulse must satisfy td_tstart <= td_supsine_tstart < "
                "td_supsine_tend <= td_tend.");
        }
    }
}

void ReadInput::item_rt_tddft()
{
    // NOTE: The order of add_item() calls below determines the parameter order
    // in the generated documentation (docs/advanced/input_files/input-main.md).
    // Please preserve this ordering when adding new parameters.
    {
        Input_Item item("estep_per_md");
        item.annotation = "steps of force change";
        item.category = "RT-TDDFT: Real-Time Time-Dependent Density Functional Theory";
        item.type = "Integer";
        item.description = "The number of electronic propagation steps between two ionic steps.";
        item.default_value = "1";
        item.unit = "";
        read_sync_int(input.estep_per_md);
        this->add_item(item);
    }
    // real time TDDFT
    {
        Input_Item item("td_dt");
        item.annotation = "time step for evolving wavefunction";
        item.category = "RT-TDDFT: Real-Time Time-Dependent Density Functional Theory";
        item.type = "Real";
        item.description = R"(The time step used for electronic propagation. If td_dt is not specified, it is set to md_dt / estep_per_md. If td_dt is specified explicitly, md_dt is reset to td_dt * estep_per_md.)";
        item.default_value = "md_dt / estep_per_md";
        item.unit = "fs";
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.td_dt == -1.0)
            {
                GlobalV::ofs_running << "td_dt don't exist, set td_dt with md_dt" << std::endl;
                para.input.td_dt = para.input.mdp.md_dt / para.input.estep_per_md;
            }
        };
        read_sync_double(input.td_dt);
        this->add_item(item);
    }
    {
        Input_Item item("td_edm");
        item.annotation = "the method to calculate the energy density matrix";
        item.category = "RT-TDDFT: Real-Time Time-Dependent Density Functional Theory";
        item.type = "Integer";
        item.description = R"(Method used to calculate the energy-density matrix for the overlap contribution to forces in LCAO RT-TDDFT.
* 0: Use $\mathrm{EDM}_{\boldsymbol{k}}=\frac{1}{2}\left(S_{\boldsymbol{k}}^{-1}H_{\boldsymbol{k}}\rho_{\boldsymbol{k}}+\rho_{\boldsymbol{k}}H_{\boldsymbol{k}}S_{\boldsymbol{k}}^{-1}\right)$.
* 1: Use the ground-state eigenvalue-weighted expression $\mathrm{EDM}_{\mu\nu,\boldsymbol{k}}=\sum_i w_{i\boldsymbol{k}}\epsilon_{i\boldsymbol{k}}C_{\mu i,\boldsymbol{k}}C_{\nu i,\boldsymbol{k}}^*$. This expression is deprecated for RT-TDDFT and is generally not valid when the propagated wave functions are not Hamiltonian eigenstates.)";
        item.default_value = "0";
        item.unit = "";
        read_sync_int(input.td_edm);
        this->add_item(item);
    }
    {
        Input_Item item("td_print_eij");
        item.annotation = "print eij or not";
        item.category = "RT-TDDFT: Real-Time Time-Dependent Density Functional Theory";
        item.type = "Real";
        item.description = R"(Controls output of the propagated-state Hamiltonian matrix elements $E_{ij}=\Braket{\psi_i | \hat{H} | \psi_j}$ to the running log. The printed band indices $i$ and $j$ are one-based global indices. Both the threshold and the printed matrix elements are in Ry.
* $\lt 0$: Disable the output.
* $\geqslant 0$: Print an element when either $\left|\operatorname{Re}E_{ij}\right|$ or $\left|\operatorname{Im}E_{ij}\right|$ is greater than or equal to td_print_eij.)";
        item.default_value = "-1";
        item.unit = "Ry";
        read_sync_double(input.td_print_eij);
        this->add_item(item);
    }
    {
        Input_Item item("td_propagator");
        item.annotation = "method of propagator";
        item.category = "RT-TDDFT: Real-Time Time-Dependent Density Functional Theory";
        item.type = "Integer";
        item.description = R"(Method used to propagate the electronic states in a nonorthogonal LCAO basis. The formulas below use Hartree atomic units, with $S$, $H$, and $\Delta t=\mathtt{td\_dt}$ evaluated as required by each approximation.
* 0: Crank-Nicolson through an explicitly constructed evolution matrix, $U=\left[S+\mathrm{i}H\Delta t/2\right]^{-1}\left[S-\mathrm{i}H\Delta t/2\right]$.
* 1: Fourth-order Taylor approximation to the exponential. With $\mathcal{A}=-\mathrm{i}S^{-1}H\Delta t$, $U=I+\mathcal{A}+\mathcal{A}^2/2+\mathcal{A}^3/6+\mathcal{A}^4/24$.
* 2: Enforced time-reversal symmetry (ETRS), $U(t+\Delta t,t)=\exp\left[-\mathrm{i}S^{-1}H(t+\Delta t)\Delta t/2\right]\exp\left[-\mathrm{i}S^{-1}H(t)\Delta t/2\right]$. In the implementation, each exponential is replaced by the fourth-order Taylor polynomial from method 1 evaluated with a half time step.
* 3: Crank-Nicolson by directly solving $\left[S+\mathrm{i}H\Delta t/2\right]\psi(t+\Delta t)=\left[S-\mathrm{i}H\Delta t/2\right]\psi(t)$.

[NOTE] GPU execution currently supports only method 0 in both single-GPU and multi-GPU solver configurations. CPU execution supports methods 0 through 3.)";
        item.default_value = "0";
        item.unit = "";
        read_sync_int(input.propagator);
        this->add_item(item);
    }
    {
        Input_Item item("td_vext");
        item.annotation = "add extern potential or not";
        item.category = "RT-TDDFT: Real-Time Time-Dependent Density Functional Theory";
        item.type = "Boolean";
        item.description = R"(Controls whether a time-dependent external electric field is applied.
* True: Add a laser-material interaction (external electric field).
* False: No external electric field.)";
        item.default_value = "False";
        item.unit = "";
        read_sync_bool(input.td_vext);
        this->add_item(item);
    }
    {
        Input_Item item("td_vext_dire");
        item.annotation = "extern potential direction";
        item.category = "RT-TDDFT: Real-Time Time-Dependent Density Functional Theory";
        item.type = "Vector of Integer";
        item.description = R"(Specifies one absolute Cartesian direction for each external electric field when td_vext is enabled. Unlike the ground-state efield_dir parameter, these directions are not defined by lattice or reciprocal-lattice vectors. The number of values must equal that of td_ttype, and repeated directions are allowed; fields assigned to the same direction are added. For example, td_vext_dire 1 2 applies one field along Cartesian x and one along Cartesian y.
* 1: The external field direction is along the x-axis.
* 2: The external field direction is along the y-axis.
* 3: The external field direction is along the z-axis.)";
        item.default_value = "1";
        item.unit = "";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            parse_expression(item.str_values, para.input.td_vext_dire);
        };
        sync_intvec(input.td_vext_dire, para.input.td_vext_dire.size(), 0);
        this->add_item(item);
    }
    {
        Input_Item item("td_stype");
        item.annotation = "type of electric field in space domain";
        item.category = "RT-TDDFT: Real-Time Time-Dependent Density Functional Theory";
        item.type = "Integer";
        item.description = R"(Type of electric field in the space domain, i.e. the gauge of the electric field.
* 0: Length gauge.
* 1: Velocity gauge.
* 2: Hybrid gauge. See J. Chem. Theory Comput. 2025, 21, 3335-3341 for more information.)";
        item.default_value = "0";
        item.unit = "";
        read_sync_int(input.td_stype);
        this->add_item(item);
    }
    {
        Input_Item item("td_ttype");
        item.annotation = "type of electric field in time domain";
        item.category = "RT-TDDFT: Real-Time Time-Dependent Density Functional Theory";
        item.type = "Vector of Integer";
        item.description = R"(Specifies one time-domain type for each external electric field. Its number of values must equal that of td_vext_dire. Parameters belonging to each type must provide exactly one value for every occurrence of that type, in occurrence order; fields with a repeated direction are added.

The formulas below use Hartree atomic units. For every ordinary input frequency $f$, $\omega=2\pi f$; $\Delta t=\mathtt{td\_dt}$; and $E_0$ denotes the corresponding amplitude parameter. A step-valued parameter $n_q$ represents the physical time $t_q=n_q\Delta t$.
* 0: Gaussian pulse, $E(t)=E_0\cos\left[\omega(t-t_0)+\varphi\right]\mathrm{e}^{-(t-t_0)^2/(2\sigma^2)}$, where $t_0=\mathtt{td\_gauss\_t0}\Delta t$.
* 1: Trapezoid pulse, $E(t)=E_0g(t)\cos(\omega t+\varphi)$. With $t_1=\mathtt{td\_trape\_t1}\Delta t$, $t_2=\mathtt{td\_trape\_t2}\Delta t$, and $t_3=\mathtt{td\_trape\_t3}\Delta t$, the envelope is $g(t)=t/t_1$ for $0\leqslant t\lt t_1$, $g(t)=1$ for $t_1\leqslant t\lt t_2$, $g(t)=(t_3-t)/(t_3-t_2)$ for $t_2\leqslant t\lt t_3$, and $g(t)=0$ otherwise.
* 2: Trigonometric pulse, $E(t)=E_0\cos(\omega_1t+\varphi_1)\sin^2(\omega_2t+\varphi_2)$.
* 3: Heaviside pulse defined on electronic steps. With $n_0=\mathtt{td\_heavi\_t0}$, $E(n)=E_0$ for $n\lt n_0$ and $E(n)=0$ for $n\geqslant n_0$.
* 4: Finite-support supersine pulse. For $t_{\mathrm{s}}\lt t\lt t_{\mathrm{e}}$, the envelope is $f(t)=\left\{\sin\left[\pi\frac{t-t_{\mathrm{s}}}{t_{\mathrm{e}}-t_{\mathrm{s}}}\right]\right\}^{\frac{\pi}{\sigma}\left|\frac{t-t_{\mathrm{s}}}{t_{\mathrm{e}}-t_{\mathrm{s}}}-\frac{1}{2}\right|}$ and the electric field is $E(t)=E_0\left\{f(t)\cos\left[\omega\left(t-\frac{t_{\mathrm{s}}+t_{\mathrm{e}}}{2}\right)+\varphi\right]+\frac{\dot{f}(t)}{\omega}\sin\left[\omega\left(t-\frac{t_{\mathrm{s}}+t_{\mathrm{e}}}{2}\right)+\varphi\right]\right\}$. The corresponding analytic vector potential is $\boldsymbol{A}(t)=-\frac{E_0}{\omega}f(t)\sin\left[\omega\left(t-\frac{t_{\mathrm{s}}+t_{\mathrm{e}}}{2}\right)+\varphi\right]\hat{\boldsymbol{e}}$, with $\boldsymbol{E}(t)=-\partial\boldsymbol{A}(t)/\partial t$. The envelope, electric field, and vector potential are zero at the pulse boundaries and outside the interval.

In the velocity and hybrid gauges, ABACUS obtains the vector potential actually used in propagation by Simpson integration of the selected electric fields, including the supersine field, so a residual at the numerical-quadrature accuracy scale may remain.)";
        item.default_value = "0";
        item.unit = "";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            parse_expression(item.str_values, para.input.td_ttype);
        };
        sync_intvec(input.td_ttype, para.input.td_ttype.size(), 0);
        item.check_value = [](const Input_Item&, const Parameter& para) { check_td_efield_parameters(para.inp); };
        this->add_item(item);
    }
    {
        Input_Item item("td_tstart");
        item.annotation = " number of steps where electric field starts";
        item.category = "RT-TDDFT: Real-Time Time-Dependent Density Functional Theory";
        item.type = "Integer";
        item.description = R"(First electronic step at which the time-dependent electric field is active. The interval from td_tstart through td_tend includes both endpoints. On each active step $n$, the velocity and hybrid gauges integrate the field over $[n\Delta t,(n+1)\Delta t]$, where $\Delta t=\mathtt{td\_dt}$.)";
        item.default_value = "1";
        item.unit = "";
        read_sync_int(input.td_tstart);
        this->add_item(item);
    }
    {
        Input_Item item("td_tend");
        item.annotation = "number of steps where electric field ends";
        item.category = "RT-TDDFT: Real-Time Time-Dependent Density Functional Theory";
        item.type = "Integer";
        item.description = R"(Last electronic step at which the time-dependent electric field is active. The interval from td_tstart through td_tend includes both endpoints. On each active step $n$, the velocity and hybrid gauges integrate the field over $[n\Delta t,(n+1)\Delta t]$, where $\Delta t=\mathtt{td\_dt}$.)";
        item.default_value = "1000";
        item.unit = "";
        read_sync_int(input.td_tend);
        this->add_item(item);
    }
    {
        Input_Item item("td_lcut1");
        item.annotation = "cut1 of interval in length gauge";
        item.category = "RT-TDDFT: Real-Time Time-Dependent Density Functional Theory";
        item.type = "Real";
        item.description = R"(Lower fractional-coordinate cutoff for the periodic spatial modulation used in the length gauge. Let $c_1=\mathtt{td\_lcut1}$, $c_2=\mathtt{td\_lcut2}$, $D=c_2-c_1$, and $G=c_1+1-c_2$. For a fractional coordinate $x$, the field factor is $\eta(x)=1$ when $c_1\leqslant x\lt c_2$ and $\eta(x)=-D/G$ elsewhere. The reversed outer interval makes the potential periodic and continuous and gives the field zero cell average.)";
        item.default_value = "0.05";
        item.unit = "";
        read_sync_double(input.td_lcut1);
        this->add_item(item);
    }
    {
        Input_Item item("td_lcut2");
        item.annotation = "cut2 of interval in length gauge";
        item.category = "RT-TDDFT: Real-Time Time-Dependent Density Functional Theory";
        item.type = "Real";
        item.description = R"(Upper fractional-coordinate cutoff for the periodic spatial modulation used in the length gauge. Let $c_1=\mathtt{td\_lcut1}$, $c_2=\mathtt{td\_lcut2}$, $D=c_2-c_1$, and $G=c_1+1-c_2$. For a fractional coordinate $x$, the field factor is $\eta(x)=1$ when $c_1\leqslant x\lt c_2$ and $\eta(x)=-D/G$ elsewhere. The reversed outer interval makes the potential periodic and continuous and gives the field zero cell average.)";
        item.default_value = "0.95";
        item.unit = "";
        read_sync_double(input.td_lcut2);
        this->add_item(item);
    }
    {
        Input_Item item("td_gauss_freq");
        item.annotation = "frequency (freq) of Gauss type electric field";
        item.category = "RT-TDDFT: Real-Time Time-Dependent Density Functional Theory";
        item.type = "Vector of Real";
        item.description = R"(Ordinary frequency $f$ in the Gaussian-pulse formula, with $\omega=2\pi f$. Supply exactly one value for each td_ttype 0 occurrence, in occurrence order.)";
        item.default_value = "22.13";
        item.unit = "1/fs";
        item.set_availability("0 in td_ttype");
        item.read_value = [](const Input_Item& item, Parameter& para) {
            parse_expression(item.str_values, para.input.td_gauss_freq);
        };
        sync_doublevec(input.td_gauss_freq, para.input.td_gauss_freq.size(), 0.0);
        this->add_item(item);
    }
    {
        Input_Item item("td_gauss_phase");
        item.annotation = "phase of Gauss type electric field";
        item.category = "RT-TDDFT: Real-Time Time-Dependent Density Functional Theory";
        item.type = "Vector of Real";
        item.description = R"(Carrier phase $\varphi$ in the Gaussian-pulse formula. Supply exactly one value for each td_ttype 0 occurrence, in occurrence order.)";
        item.default_value = "0.0";
        item.unit = "rad";
        item.set_availability("0 in td_ttype");
        item.read_value = [](const Input_Item& item, Parameter& para) {
            parse_expression(item.str_values, para.input.td_gauss_phase);
        };
        sync_doublevec(input.td_gauss_phase, para.input.td_gauss_phase.size(), 0.0);
        this->add_item(item);
    }
    {
        Input_Item item("td_gauss_sigma");
        item.annotation = "sigma of Gauss type electric field";
        item.category = "RT-TDDFT: Real-Time Time-Dependent Density Functional Theory";
        item.type = "Vector of Real";
        item.description = R"(Nonzero standard deviation $\sigma$ of the Gaussian envelope. Supply exactly one value for each td_ttype 0 occurrence, in occurrence order.)";
        item.default_value = "30.0";
        item.unit = "fs";
        item.set_availability("0 in td_ttype");
        item.read_value = [](const Input_Item& item, Parameter& para) {
            parse_expression(item.str_values, para.input.td_gauss_sigma);
        };
        sync_doublevec(input.td_gauss_sigma, para.input.td_gauss_sigma.size(), 0.0);
        this->add_item(item);
    }
    {
        Input_Item item("td_gauss_t0");
        item.annotation = "step number of time center (t0) of Gauss type electric field";
        item.category = "RT-TDDFT: Real-Time Time-Dependent Density Functional Theory";
        item.type = "Vector of Real";
        item.description = R"(Electronic-step position of the Gaussian center, which defines $t_0=\mathtt{td\_gauss\_t0}\Delta t$. Supply exactly one value for each td_ttype 0 occurrence, in occurrence order.)";
        item.default_value = "100";
        item.unit = "";
        item.set_availability("0 in td_ttype");
        item.read_value = [](const Input_Item& item, Parameter& para) {
            parse_expression(item.str_values, para.input.td_gauss_t0);
        };
        sync_doublevec(input.td_gauss_t0, para.input.td_gauss_t0.size(), 0.0);
        this->add_item(item);
    }
    {
        Input_Item item("td_gauss_amp");
        item.annotation = "amplitude of Gauss type electric field";
        item.category = "RT-TDDFT: Real-Time Time-Dependent Density Functional Theory";
        item.type = "Vector of Real";
        item.description = R"(Electric-field scale $E_0$ in the Gaussian-pulse formula. Supply exactly one value for each td_ttype 0 occurrence, in occurrence order.)";
        item.default_value = "0.25";
        item.unit = "V/Angstrom";
        item.set_availability("0 in td_ttype");
        item.read_value = [](const Input_Item& item, Parameter& para) {
            parse_expression(item.str_values, para.input.td_gauss_amp);
        };
        sync_doublevec(input.td_gauss_amp, para.input.td_gauss_amp.size(), 0.0);
        this->add_item(item);
    }
    {
        Input_Item item("td_trape_freq");
        item.annotation = "frequency of Trapezoid type electric field";
        item.category = "RT-TDDFT: Real-Time Time-Dependent Density Functional Theory";
        item.type = "Vector of Real";
        item.description = R"(Ordinary carrier frequency $f$ in the trapezoid-pulse formula, with $\omega=2\pi f$. Supply exactly one value for each td_ttype 1 occurrence, in occurrence order.)";
        item.default_value = "1.60";
        item.unit = "1/fs";
        item.set_availability("1 in td_ttype");
        item.read_value = [](const Input_Item& item, Parameter& para) {
            parse_expression(item.str_values, para.input.td_trape_freq);
        };
        sync_doublevec(input.td_trape_freq, para.input.td_trape_freq.size(), 0.0);
        this->add_item(item);
    }
    {
        Input_Item item("td_trape_phase");
        item.annotation = "phase of Trapezoid type electric field";
        item.category = "RT-TDDFT: Real-Time Time-Dependent Density Functional Theory";
        item.type = "Vector of Real";
        item.description = R"(Carrier phase $\varphi$ in the trapezoid-pulse formula. Supply exactly one value for each td_ttype 1 occurrence, in occurrence order.)";
        item.default_value = "0.0";
        item.unit = "rad";
        item.set_availability("1 in td_ttype");
        item.read_value = [](const Input_Item& item, Parameter& para) {
            parse_expression(item.str_values, para.input.td_trape_phase);
        };
        sync_doublevec(input.td_trape_phase, para.input.td_trape_phase.size(), 0.0);
        this->add_item(item);
    }
    {
        Input_Item item("td_trape_t1");
        item.annotation = "t1 of Trapezoid type electric field";
        item.category = "RT-TDDFT: Real-Time Time-Dependent Density Functional Theory";
        item.type = "Vector of Real";
        item.description = R"(Electronic step defining the end of the linear rise, $t_1=\mathtt{td\_trape\_t1}\Delta t$. Each field must satisfy td_trape_t1 <= td_trape_t2 <= td_trape_t3. Supply exactly one value for each td_ttype 1 occurrence, in occurrence order.)";
        item.default_value = "1875";
        item.unit = "";
        item.set_availability("1 in td_ttype");
        item.read_value = [](const Input_Item& item, Parameter& para) {
            parse_expression(item.str_values, para.input.td_trape_t1);
        };
        sync_doublevec(input.td_trape_t1, para.input.td_trape_t1.size(), 0.0);
        this->add_item(item);
    }
    {
        Input_Item item("td_trape_t2");
        item.annotation = "t2 of Trapezoid type electric field";
        item.category = "RT-TDDFT: Real-Time Time-Dependent Density Functional Theory";
        item.type = "Vector of Real";
        item.description = R"(Electronic step defining the end of the plateau, $t_2=\mathtt{td\_trape\_t2}\Delta t$. Each field must satisfy td_trape_t1 <= td_trape_t2 <= td_trape_t3. Supply exactly one value for each td_ttype 1 occurrence, in occurrence order.)";
        item.default_value = "5625";
        item.unit = "";
        item.set_availability("1 in td_ttype");
        item.read_value = [](const Input_Item& item, Parameter& para) {
            parse_expression(item.str_values, para.input.td_trape_t2);
        };
        sync_doublevec(input.td_trape_t2, para.input.td_trape_t2.size(), 0.0);
        this->add_item(item);
    }
    {
        Input_Item item("td_trape_t3");
        item.annotation = "t3 of Trapezoid type electric field";
        item.category = "RT-TDDFT: Real-Time Time-Dependent Density Functional Theory";
        item.type = "Vector of Real";
        item.description = R"(Electronic step defining the end of the linear fall, $t_3=\mathtt{td\_trape\_t3}\Delta t$. Each field must satisfy td_trape_t1 <= td_trape_t2 <= td_trape_t3. Supply exactly one value for each td_ttype 1 occurrence, in occurrence order.)";
        item.default_value = "7500";
        item.unit = "";
        item.set_availability("1 in td_ttype");
        item.read_value = [](const Input_Item& item, Parameter& para) {
            parse_expression(item.str_values, para.input.td_trape_t3);
        };
        sync_doublevec(input.td_trape_t3, para.input.td_trape_t3.size(), 0.0);
        this->add_item(item);
    }
    {
        Input_Item item("td_trape_amp");
        item.annotation = "amplitude of Trapezoid type electric field";
        item.category = "RT-TDDFT: Real-Time Time-Dependent Density Functional Theory";
        item.type = "Vector of Real";
        item.description = R"(Electric-field scale $E_0$ in the trapezoid-pulse formula. Supply exactly one value for each td_ttype 1 occurrence, in occurrence order.)";
        item.default_value = "2.74";
        item.unit = "V/Angstrom";
        item.set_availability("1 in td_ttype");
        item.read_value = [](const Input_Item& item, Parameter& para) {
            parse_expression(item.str_values, para.input.td_trape_amp);
        };
        sync_doublevec(input.td_trape_amp, para.input.td_trape_amp.size(), 0.0);
        this->add_item(item);
    }
    {
        Input_Item item("td_trigo_freq1");
        item.annotation = "frequency 1 of Trigonometric type electric field";
        item.category = "RT-TDDFT: Real-Time Time-Dependent Density Functional Theory";
        item.type = "Vector of Real";
        item.description = R"(First ordinary frequency $f_1$ in the trigonometric-pulse formula, with $\omega_1=2\pi f_1$. Supply exactly one value for each td_ttype 2 occurrence, in occurrence order.)";
        item.default_value = "1.164656";
        item.unit = "1/fs";
        item.set_availability("2 in td_ttype");
        item.read_value = [](const Input_Item& item, Parameter& para) {
            parse_expression(item.str_values, para.input.td_trigo_freq1);
        };
        sync_doublevec(input.td_trigo_freq1, para.input.td_trigo_freq1.size(), 0.0);
        this->add_item(item);
    }
    {
        Input_Item item("td_trigo_freq2");
        item.annotation = "frequency 2 of Trigonometric type electric field";
        item.category = "RT-TDDFT: Real-Time Time-Dependent Density Functional Theory";
        item.type = "Vector of Real";
        item.description = R"(Second ordinary frequency $f_2$ in the trigonometric-pulse formula, with $\omega_2=2\pi f_2$. Supply exactly one value for each td_ttype 2 occurrence, in occurrence order.)";
        item.default_value = "0.029116";
        item.unit = "1/fs";
        item.set_availability("2 in td_ttype");
        item.read_value = [](const Input_Item& item, Parameter& para) {
            parse_expression(item.str_values, para.input.td_trigo_freq2);
        };
        sync_doublevec(input.td_trigo_freq2, para.input.td_trigo_freq2.size(), 0.0);
        this->add_item(item);
    }
    {
        Input_Item item("td_trigo_phase1");
        item.annotation = "phase 1 of Trigonometric type electric field";
        item.category = "RT-TDDFT: Real-Time Time-Dependent Density Functional Theory";
        item.type = "Vector of Real";
        item.description = R"(Carrier phase $\varphi_1$ in the cosine factor of the trigonometric-pulse formula. Supply exactly one value for each td_ttype 2 occurrence, in occurrence order.)";
        item.default_value = "0.0";
        item.unit = "rad";
        item.set_availability("2 in td_ttype");
        item.read_value = [](const Input_Item& item, Parameter& para) {
            parse_expression(item.str_values, para.input.td_trigo_phase1);
        };
        sync_doublevec(input.td_trigo_phase1, para.input.td_trigo_phase1.size(), 0.0);
        this->add_item(item);
    }
    {
        Input_Item item("td_trigo_phase2");
        item.annotation = "phase 2 of Trigonometric type electric field";
        item.category = "RT-TDDFT: Real-Time Time-Dependent Density Functional Theory";
        item.type = "Vector of Real";
        item.description = R"(Envelope phase $\varphi_2$ in the sine-squared factor of the trigonometric-pulse formula. Supply exactly one value for each td_ttype 2 occurrence, in occurrence order.)";
        item.default_value = "0.0";
        item.unit = "rad";
        item.set_availability("2 in td_ttype");
        item.read_value = [](const Input_Item& item, Parameter& para) {
            parse_expression(item.str_values, para.input.td_trigo_phase2);
        };
        sync_doublevec(input.td_trigo_phase2, para.input.td_trigo_phase2.size(), 0.0);
        this->add_item(item);
    }
    {
        Input_Item item("td_trigo_amp");
        item.annotation = "amplitude of Trigonometric type electric field";
        item.category = "RT-TDDFT: Real-Time Time-Dependent Density Functional Theory";
        item.type = "Vector of Real";
        item.description = R"(Electric-field scale $E_0$ in the trigonometric-pulse formula. Supply exactly one value for each td_ttype 2 occurrence, in occurrence order.)";
        item.default_value = "2.74";
        item.unit = "V/Angstrom";
        item.set_availability("2 in td_ttype");
        item.read_value = [](const Input_Item& item, Parameter& para) {
            parse_expression(item.str_values, para.input.td_trigo_amp);
        };
        sync_doublevec(input.td_trigo_amp, para.input.td_trigo_amp.size(), 0.0);
        this->add_item(item);
    }
    {
        Input_Item item("td_heavi_t0");
        item.annotation = "t0 of Heaviside type electric field";
        item.category = "RT-TDDFT: Real-Time Time-Dependent Density Functional Theory";
        item.type = "Vector of Real";
        item.description = R"(Electronic switch step $n_0$ in the Heaviside-pulse definition. The field is $E_0$ for $n\lt n_0$ and zero for $n\geqslant n_0$. Supply exactly one value for each td_ttype 3 occurrence, in occurrence order.)";
        item.default_value = "100";
        item.unit = "";
        item.set_availability("3 in td_ttype");
        item.read_value = [](const Input_Item& item, Parameter& para) {
            parse_expression(item.str_values, para.input.td_heavi_t0);
        };
        sync_doublevec(input.td_heavi_t0, para.input.td_heavi_t0.size(), 0.0);
        this->add_item(item);
    }
    {
        Input_Item item("td_heavi_amp");
        item.annotation = "amplitude of Heaviside type electric field";
        item.category = "RT-TDDFT: Real-Time Time-Dependent Density Functional Theory";
        item.type = "Vector of Real";
        item.description = R"(Electric-field scale $E_0$ in the Heaviside-pulse definition. Supply exactly one value for each td_ttype 3 occurrence, in occurrence order.)";
        item.default_value = "1.0";
        item.unit = "V/Angstrom";
        item.set_availability("3 in td_ttype");
        item.read_value = [](const Input_Item& item, Parameter& para) {
            parse_expression(item.str_values, para.input.td_heavi_amp);
        };
        sync_doublevec(input.td_heavi_amp, para.input.td_heavi_amp.size(), 0.0);
        this->add_item(item);
    }
    {
        Input_Item item("td_supsine_amp");
        item.annotation = "carrier electric-field scale of the supersine pulse";
        item.category = "RT-TDDFT: Real-Time Time-Dependent Density Functional Theory";
        item.type = "Vector of Real";
        item.description = R"(Carrier electric-field scale $E_0$ of each supersine pulse. This is not a normalization of the complete waveform maximum, because the envelope-derivative term also contributes. Supply exactly one value for each td_ttype 4 occurrence, in occurrence order.)";
        item.default_value = "0.27";
        item.unit = "V/Angstrom";
        item.set_availability("4 in td_ttype");
        item.read_value = [](const Input_Item& item, Parameter& para) {
            parse_expression(item.str_values, para.input.td_supsine_amp);
        };
        sync_doublevec(input.td_supsine_amp, para.input.td_supsine_amp.size(), 0.0);
        this->add_item(item);
    }
    {
        Input_Item item("td_supsine_freq");
        item.annotation = "carrier frequency of the supersine pulse";
        item.category = "RT-TDDFT: Real-Time Time-Dependent Density Functional Theory";
        item.type = "Vector of Real";
        item.description = R"(Nonzero ordinary carrier frequency $f$ of each supersine pulse, with $\omega=2\pi f$. Supply exactly one value for each td_ttype 4 occurrence, in occurrence order.)";
        item.default_value = "0.18737028625";
        item.unit = "1/fs";
        item.set_availability("4 in td_ttype");
        item.read_value = [](const Input_Item& item, Parameter& para) {
            parse_expression(item.str_values, para.input.td_supsine_freq);
        };
        sync_doublevec(input.td_supsine_freq, para.input.td_supsine_freq.size(), 0.0);
        this->add_item(item);
    }
    {
        Input_Item item("td_supsine_phase");
        item.annotation = "carrier phase at the center of the supersine pulse";
        item.category = "RT-TDDFT: Real-Time Time-Dependent Density Functional Theory";
        item.type = "Vector of Real";
        item.description = R"(Electric-field carrier phase $\varphi$ at the center of each supersine envelope. A value of 0 places a cosine carrier maximum at the envelope center. Supply exactly one value for each td_ttype 4 occurrence, in occurrence order.)";
        item.default_value = "0.0";
        item.unit = "rad";
        item.set_availability("4 in td_ttype");
        item.read_value = [](const Input_Item& item, Parameter& para) {
            parse_expression(item.str_values, para.input.td_supsine_phase);
        };
        sync_doublevec(input.td_supsine_phase, para.input.td_supsine_phase.size(), 0.0);
        this->add_item(item);
    }
    {
        Input_Item item("td_supsine_sigma");
        item.annotation = "shape parameter of the supersine envelope";
        item.category = "RT-TDDFT: Real-Time Time-Dependent Density Functional Theory";
        item.type = "Vector of Real";
        item.description = R"(Dimensionless shape parameter $\sigma$ of each supersine envelope. It must satisfy $0\lt\sigma\lt\pi/2$ so that the electric field approaches zero at the pulse boundaries. Supply exactly one value for each td_ttype 4 occurrence, in occurrence order.)";
        item.default_value = "0.75";
        item.unit = "";
        item.set_availability("4 in td_ttype");
        item.read_value = [](const Input_Item& item, Parameter& para) {
            parse_expression(item.str_values, para.input.td_supsine_sigma);
        };
        sync_doublevec(input.td_supsine_sigma, para.input.td_supsine_sigma.size(), 0.0);
        this->add_item(item);
    }
    {
        Input_Item item("td_supsine_tstart");
        item.annotation = "start boundary step of the supersine pulse";
        item.category = "RT-TDDFT: Real-Time Time-Dependent Density Functional Theory";
        item.type = "Vector of String";
        item.description = R"(Integer electronic step at the left, exactly zero boundary of each supersine pulse, defining $t_{\mathrm{s}}=\mathtt{td\_supsine\_tstart}\Delta t$. Supply exactly one integer or default token for each td_ttype 4 occurrence, in occurrence order; each default token inherits td_tstart. The complete pulse support must lie inside the inclusive global td_tstart to td_tend interval; hard truncation of a supersine pulse is rejected.)";
        item.default_value = "default";
        item.unit = "";
        item.set_availability("4 in td_ttype");
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            para.input.td_supsine_tstart = parse_supersine_steps(item, para.input.td_tstart);
        };
        item.get_final_value = [](Input_Item& item, const Parameter&) {
            item.final_value << (item.is_read() ? longstring(item.str_values) : "default");
        };
        add_intvec_bcast(input.td_supsine_tstart, para.input.td_supsine_tstart.size(), 0);
        this->add_item(item);
    }
    {
        Input_Item item("td_supsine_tend");
        item.annotation = "end boundary step of the supersine pulse";
        item.category = "RT-TDDFT: Real-Time Time-Dependent Density Functional Theory";
        item.type = "Vector of String";
        item.description = R"(Integer electronic step at the right, exactly zero boundary of each supersine pulse, defining $t_{\mathrm{e}}=\mathtt{td\_supsine\_tend}\Delta t$. Supply exactly one integer or default token for each td_ttype 4 occurrence, in occurrence order; each default token inherits td_tend. The complete pulse support must lie inside the inclusive global td_tstart to td_tend interval; hard truncation of a supersine pulse is rejected.)";
        item.default_value = "default";
        item.unit = "";
        item.set_availability("4 in td_ttype");
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            para.input.td_supsine_tend = parse_supersine_steps(item, para.input.td_tend);
        };
        item.get_final_value = [](Input_Item& item, const Parameter&) {
            item.final_value << (item.is_read() ? longstring(item.str_values) : "default");
        };
        add_intvec_bcast(input.td_supsine_tend, para.input.td_supsine_tend.size(), 0);
        this->add_item(item);
    }
    {
        Input_Item item("init_vecpot_file");
        item.annotation = "init vector potential through file or not";
        item.category = "RT-TDDFT: Real-Time Time-Dependent Density Functional Theory";
        item.type = "Boolean";
        item.description = R"(Selects the source of the Cartesian vector potential used by LCAO RT-TDDFT.
* True: Read vector_pot.txt from the calculation working directory. Each non-comment line must contain four columns: a conventionally one-based electronic-step label followed by $A_x$, $A_y$, and $A_z$ in atomic units. Rows are consumed sequentially; the first column is read as a label and is not used for lookup. If propagation continues beyond the available rows, the last row is reused.
* False: Obtain the vector potential by integrating the configured electric field.)";
        item.default_value = "False";
        item.unit = "";
        read_sync_bool(input.init_vecpot_file);
        this->add_item(item);
    }
    {
        Input_Item item("ocp");
        item.annotation = "change occupation or not";
        item.category = "RT-TDDFT: Real-Time Time-Dependent Density Functional Theory";
        item.type = "Boolean";
        item.description = R"(Controls fixed band occupations. In calculations other than LCAO RT-TDDFT, fixed values are applied during electronic-state setup. In LCAO RT-TDDFT, the initial ground-state SCF determines occupations normally, and fixed values from ocp_set are applied during the subsequent real-time propagation steps.
* True: Use the fixed occupations specified by ocp_set during propagation.
* False: Keep the occupations determined by the initial SCF.)";
        item.default_value = "False";
        item.unit = "";
        read_sync_bool(input.ocp);
        this->add_item(item);
    }
    {
        Input_Item item("ocp_set");
        item.annotation = "set occupation";
        item.category = "RT-TDDFT: Real-Time Time-Dependent Density Functional Theory";
        item.type = "String";
        item.description = R"(Fixed occupation weights used when ocp is true. Values are assigned band by band for each k-point, following k-point order. In LCAO RT-TDDFT, the initial ground-state SCF uses its normally determined occupations, and this array is applied only during subsequent real-time propagation steps. The repetition syntax N*x expands to N copies of x.
* Example: 1 10*1 0 1 expands to 13 values, with the 12th value equal to 0 and all other values equal to 1.
* After expansion, the array length must equal nks * nbands.
* The sum of all weights must equal nelec; otherwise the calculation terminates with an error.)";
        item.default_value = "None";
        item.unit = "";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            parse_expression(item.str_values, para.input.ocp_kb);
        };
        item.get_final_value = [](Input_Item& item, const Parameter& para) {
            if(item.is_read())
            {
                item.final_value.str(longstring(item.str_values));
            }
        };
        add_doublevec_bcast(input.ocp_kb, para.input.ocp_kb.size(), 0.0);
        this->add_item(item);
    }
}
void ReadInput::item_tdofdft()
{
    // NOTE: The order of add_item() calls below determines the parameter order
    // in the generated documentation (docs/advanced/input_files/input-main.md).
    // Please preserve this ordering when adding new parameters.
    // TD-OFDFT
    {
        Input_Item item("of_cd");
        item.annotation = "add CD Potential or not";
        item.category = "TDOFDFT: time dependent orbital free density functional theory";
        item.type = "Boolean";
        item.description = R"(Added the current dependent(CD) potential. (https://doi.org/10.1103/PhysRevB.98.144302)
* True: Added the CD potential.
* False: Not added the CD potential.)";
        item.default_value = "False";
        item.unit = "";
        item.set_availability("esolver_type==tdofdft");
        read_sync_bool(input.of_cd);
        this->add_item(item);
    }
    {
        Input_Item item("of_mcd_alpha");
        item.annotation = "parameter of modified CD Potential";
        item.category = "TDOFDFT: time dependent orbital free density functional theory";
        item.type = "Real";
        item.description = "The value of the parameter alpha in modified CD potential method. mCDPotential=alpha*CDPotential (proposed in paper PhysRevB.98.144302)";
        item.default_value = "1.0";
        item.unit = "";
        item.set_availability("esolver_type==tdofdft");
        read_sync_double(input.of_mCD_alpha);
        this->add_item(item);
    }
}
void ReadInput::item_lr_tddft()
{
    // NOTE: The order of add_item() calls below determines the parameter order
    // in the generated documentation (docs/advanced/input_files/input-main.md).
    // Please preserve this ordering when adding new parameters.
    {
        Input_Item item("xc_kernel");
        item.annotation = "exchange correlation (XC) kernel for LR-TDDFT";
        item.category = "Linear Response TDDFT (Under Development Feature)";
        item.type = "String";
        item.description = "The exchange-correlation kernel used in the calculation. Currently supported: RPA, LDA, PBE, HSE, HF.";
        item.default_value = "LDA";
        item.unit = "";
        read_sync_string(input.xc_kernel);
        this->add_item(item);
    }
    {
        Input_Item item("lr_init_xc_kernel");
        item.annotation = "The method to initalize the xc kernel";
        item.category = "Linear Response TDDFT (Under Development Feature)";
        item.type = "Vector of String (>=1 values)";
        item.description = R"(The method to initalize the xc kernel.
* "default": Calculate xc kernel from the ground-state charge density.
* "file": Read the xc kernel on grid from the provided files. The following words should be the paths of ".cube" files, where the first 1 (nspin==1) or 3 (nspin==2, namely spin-aa, spin-ab and spin-bb) will be read in. The parameter xc_kernel will be invalid. Now only LDA-type kernel is supported as the potential will be calculated by directly multiplying the transition density.
* "from_charge_file": Calculate fxc from the charge density read from the provided files. The following words should be the paths of ".cube" files, where the first nspin files will be read in.)";
        item.default_value = "\"default\"";
        item.unit = "";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            size_t count = item.get_size();
            auto& ifxc = para.input.lr_init_xc_kernel;
            for (int i = 0; i < count; i++) { ifxc.push_back(item.str_values[i]); }
            };
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.lr_init_xc_kernel.empty()) { para.input.lr_init_xc_kernel.push_back("default"); }
            };
        sync_stringvec(input.lr_init_xc_kernel, para.input.lr_init_xc_kernel.size(), "default");
        this->add_item(item);
    }
    {
        Input_Item item("lr_solver");
        item.annotation = "the eigensolver for LR-TDDFT";
        item.category = "Linear Response TDDFT (Under Development Feature)";
        item.type = "String";
        item.description = R"(The method to solve the Casida equation in LR-TDDFT under Tamm-Dancoff approximation (TDA).
* dav/dav_subspace/cg: Construct and diagonalize the Hamiltonian matrix iteratively with Davidson/Non-ortho-Davidson/CG algorithm.
* lapack: Construct the full matrix and directly diagonalize with LAPACK.
* spectrum: Calculate absorption spectrum only without solving Casida equation.)";
        item.default_value = "dav";
        item.unit = "";
        read_sync_string(input.lr_solver);
        this->add_item(item);
    }
    {
        Input_Item item("lr_thr");
        item.annotation = "convergence threshold of the LR-TDDFT eigensolver";
        item.category = "Linear Response TDDFT (Under Development Feature)";
        item.type = "Real";
        item.description = "The convergence threshold of iterative diagonalization solver for LR-TDDFT. It is a pure-math number with the same meaning as pw_diag_thr, but since the Casida equation is a one-shot eigenvalue problem, it is also the convergence threshold of LR-TDDFT.";
        item.default_value = "1e-2";
        item.unit = "";
        read_sync_double(input.lr_thr);
        this->add_item(item);
    }
    {
        Input_Item item("nocc");
        item.annotation = "the number of occupied orbitals to form the 2-particle basis ( <= nelec/2)";
        item.category = "Linear Response TDDFT (Under Development Feature)";
        item.type = "Integer";
        item.description = R"(The number of occupied orbitals (up to HOMO) used in the LR-TDDFT calculation.
* Note: If the value is illegal ( > nelec/2 or <= 0), it will be autoset to nelec/2.)";
        item.default_value = "nband";
        item.unit = "";
        read_sync_int(input.nocc);
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            const int nocc_default = std::max(static_cast<int>(para.input.nelec + 1) / 2, para.input.nbands);
            if (para.input.nocc <= 0 || para.input.nocc > nocc_default) { para.input.nocc = nocc_default; }
            };
        this->add_item(item);
    }
    {
        Input_Item item("nvirt");
        item.annotation = "the number of virtual orbitals to form the 2-particle basis (nocc + nvirt <= nbands)";
        item.category = "Linear Response TDDFT (Under Development Feature)";
        item.type = "Integer";
        item.description = "The number of virtual orbitals (starting from LUMO) used in the LR-TDDFT calculation.";
        item.default_value = "1";
        item.unit = "";
        read_sync_int(input.nvirt);
        this->add_item(item);
    }
    // Linear Responce TDDFT
    {
        Input_Item item("lr_nstates");
        item.annotation = "the number of 2-particle states to be solved";
        item.category = "Linear Response TDDFT (Under Development Feature)";
        item.type = "Integer";
        item.description = "The number of 2-particle states to be solved.";
        item.default_value = "0";
        item.unit = "";
        read_sync_int(input.lr_nstates);
        this->add_item(item);
    }
    {
        Input_Item item("lr_unrestricted");
        item.annotation = "Whether to use unrestricted construction for LR-TDDFT";
        item.category = "Linear Response TDDFT (Under Development Feature)";
        item.type = "Boolean";
        item.description = R"(Whether to use unrestricted construction for LR-TDDFT (the matrix size will be doubled).
* True: Always use unrestricted LR-TDDFT.
* False: Use unrestricted LR-TDDFT only when the system is open-shell.)";
        item.default_value = "False";
        item.unit = "";
        read_sync_bool(input.lr_unrestricted);
        this->add_item(item);
    }
    {
        Input_Item item("abs_wavelen_range");
        item.annotation = "the range of wavelength(nm) to output the absorption spectrum ";
        item.category = "Linear Response TDDFT (Under Development Feature)";
        item.type = "Real Real";
        item.description = "The range of the wavelength for the absorption spectrum calculation.";
        item.default_value = "0.0 0.0";
        item.unit = "nm";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            size_t count = item.get_size();
            for (int i = 0; i < count; i++)
            {
                para.input.abs_wavelen_range.push_back(std::stod(item.str_values[i]));
            }
            };
        sync_doublevec(input.abs_wavelen_range, 2, 0.0);
        this->add_item(item);
    }
    {
        Input_Item item("out_wfc_lr");
        item.annotation = "whether to output the eigenvectors (excitation amplitudes) in the particle-hole basis";
        item.category = "Linear Response TDDFT (Under Development Feature)";
        item.type = "Boolean";
        item.description = "Whether to output the eigenstates (excitation energy) and eigenvectors (excitation amplitude) of the LR-TDDFT calculation. The output files are OUT.{suffix}/Excitation_Amplitude_${processor_rank}.dat.";
        item.default_value = "False";
        item.unit = "";
        read_sync_bool(input.out_wfc_lr);
        this->add_item(item);
    }
    {
        Input_Item item("abs_gauge");
        item.annotation = "whether to use length or velocity gauge to calculate the absorption spectrum in LR-TDDFT";
        item.category = "Linear Response TDDFT (Under Development Feature)";
        item.type = "String";
        item.description = "Whether to use length or velocity gauge to calculate the absorption spectrum in LR-TDDFT.";
        item.default_value = "length";
        item.unit = "";
        read_sync_string(input.abs_gauge);
        this->add_item(item);
    }
    {
        Input_Item item("abs_broadening");
        item.annotation = "the broadening (eta) for LR-TDDFT absorption spectrum";
        item.category = "Linear Response TDDFT (Under Development Feature)";
        item.type = "Real";
        item.description = "The broadening factor for the absorption spectrum calculation.";
        item.default_value = "0.01";
        item.unit = "";
        read_sync_double(input.abs_broadening);
        this->add_item(item);
    }
}
}
