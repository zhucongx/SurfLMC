#include "Atom.hpp"
#include "Element.hpp"
#include "Lattice.hpp"
#include "Config.h"
#include "Structure.h"
int main(int argc, char *argv[]) {
  cfg::Config config;

  config = cfg::CreateOrthLayers(7, 7, "CBABABABABABABAB");
  config.WriteConfig("GaN.cfg");
  config.WritePoscar("POSCAR");


  // config = cfg::CreateOrthLayers(1, 1, "ABCABC");
  // config.WriteConfig("zincblende.cfg");
  // config.WritePoscar("POSCAR_zincblende");
  //
  // config = cfg::CreateOrthLayers(1, 1, "ABABAB");
  // config.WriteConfig("wurtzite.cfg");
  // config.WritePoscar("POSCAR_wurtzite");

  return 0;
}
