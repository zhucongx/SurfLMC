#include "Atom.hpp"
#include "Element.hpp"
#include "Lattice.hpp"
#include "Config.h"
int main(int argc, char *argv[]) {
  auto config = cfg::Config::ReadConfig("zincblende.cfg");

  config.WriteConfig("zincblende_out.cfg");
  return 0;
}
