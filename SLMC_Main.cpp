#include "Atom.hpp"
#include "Element.hpp"
#include "Lattice.hpp"
#include "Config.h"
#include "Structure.h"
int main(int argc, char *argv[]) {
  cfg::Config config = cfg::CreateLayers(1, 1, "ABABABCABC");
  config.WriteConfig("out.cfg");
  config.WritePoscar("POSCAR_out");

  return 0;
}
