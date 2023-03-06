#ifndef LMC_LMC_MC_INCLUDE_KINETICMCFIRSTMPI_H_
#define LMC_LMC_MC_INCLUDE_KINETICMCFIRSTMPI_H_
#include <mpi.h>
#include "KineticMcAbstract.h"
namespace mc {
class KineticMcFirstMpi : public KineticMcFirstAbstract {
  public:
    KineticMcFirstMpi(cfg::Config config,
                      unsigned long long int log_dump_steps,
                      unsigned long long int config_dump_steps,
                      unsigned long long int maximum_steps,
                      unsigned long long int thermodynamic_averaging_steps,
                      unsigned long long int restart_steps,
                      double restart_energy,
                      double restart_time,
                      double temperature,
                      const std::set<Element> &element_set,
                      const std::string &json_coefficients_filename,
                      const std::string &time_temperature_filename,
                      bool is_rate_corrector);
    ~KineticMcFirstMpi() override;
  protected:
    void BuildEventList() override;
    double CalculateTime() override;
};
} // mc

#endif //LMC_LMC_MC_INCLUDE_KINETICMCFIRSTMPI_H_
