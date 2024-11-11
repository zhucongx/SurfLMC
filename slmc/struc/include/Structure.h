#ifndef SLMC_SLMC_STRUC_INCLUDE_STRUCTURE_H_
#define SLMC_SLMC_STRUC_INCLUDE_STRUCTURE_H_

#include "Config.h"
namespace cfg {

std::pair<std::vector<Lattice>, std::vector<Atom> > CreateOneLayerA(size_t n_x, size_t n_y, size_t layer_index);
std::pair<std::vector<Lattice>, std::vector<Atom> > CreateOneLayerB(size_t n_x, size_t n_y, size_t layer_index);
std::pair<std::vector<Lattice>, std::vector<Atom> > CreateOneLayerC(size_t n_x, size_t n_y, size_t layer_index);

std::pair<std::vector<Lattice>, std::vector<Atom> > CreateOneOrthLayerA(size_t n_x, size_t n_y, size_t layer_index);
std::pair<std::vector<Lattice>, std::vector<Atom> > CreateOneOrthLayerB(size_t n_x, size_t n_y, size_t layer_index);
std::pair<std::vector<Lattice>, std::vector<Atom> > CreateOneOrthLayerC(size_t n_x, size_t n_y, size_t layer_index);
cfg::Config CreateLayers(size_t n_x, size_t ny, const std::string &layers_type);
cfg::Config CreateOrthLayers(size_t n_x, size_t ny, const std::string &layers_type);
cfg::Config CreateOrthLayersSlab(size_t n_x, size_t ny, const std::string &layers_type);

} // cfg

#endif //SLMC_SLMC_STRUC_INCLUDE_STRUCTURE_H_
