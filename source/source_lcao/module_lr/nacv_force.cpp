#include "nacv_force.h"
#include "cal_hs_grad.h"
#include "gradient_inputs.h"
#include "nacv_timer.h"
#include "nacv_workflow.h"
namespace LR
{
NacvForce::NacvForce(const NacvInputs<double> &data, const ForceSettings &settings)
    : LR_Force<double>(
          data.gradient.ucell,
          data.gradient.kv.kvec_d,
          data.gradient.pmat,
          data.density_basis,
          data.smooth_basis,
          data.local_potential,
          data.structure_factor,
          data.gradient.gd,
          data.two_center_bundle,
          settings,
          data.gradient.ofs
#ifdef __EXX
          ,
          data.gradient.exx_lri,
          data.gradient.hybrid_alpha
#endif
          ),
      data_(data)
{
}
const ForceGeometry &NacvForce::geometry() { return *geometry_terms(); }
const ForceGeometry *NacvForce::geometry_terms()
{
    if (!geometry_)
    {
        prepare_geometry();
    }
    return geometry_.get();
}
void NacvForce::prepare_geometry()
{
    const NacvTimerScope timer("NacvForce", "prepare_geometry");
    const auto &g = data_.gradient;
    overlap_ = cal_hs_grad(
        'S',
        g.ucell,
        g.pmat,
        g.gd,
        data_.two_center_bundle);
    kinetic_ = cal_hs_grad(
        'T',
        g.ucell,
        g.pmat,
        g.gd,
        data_.two_center_bundle);
    local_ = local_potential();
    ground_ = dm_to_hxc_potential(data_.ground_density);
    geometry_.reset(new ForceGeometry{overlap_,
        kinetic_,
        *local_,
        *ground_});
}
} // namespace LR
