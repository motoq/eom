#ifndef EOM_EPHEM_PRINTER_H
#define EOM_EPHEM_PRINTER_H

#include <memory>
#include <string>
#include <vector>
#include <deque>

#include <eom_command.h>

//#include <cal_duration.h>
#include <cal_julian_date.h>
#include <astro_ephemeris.h>

namespace eom_app {

class EomEphemPrinter : public EomCommand {
public:
  // Consumes tokesn
  EomEphemPrinter(std::deque<std::string>& tokens,
      const eom::JulianDate& jdEphStart, const eom::JulianDate& jdEphStop,
      const std::shared_ptr<std::unordered_map<std::string, int>>& orbit_ndxs,
      const std::shared_ptr<std::vector<std::shared_ptr<eom::Ephemeris>>>&
                                                                  orbit_ephems);

  void execute() const override;

private:
  int vndx;
  eom::EphemFrame frame;
  std::string file_name;
  eom::JulianDate jdStart;
  eom::JulianDate jdStop;
  std::shared_ptr<std::vector<std::shared_ptr<eom::Ephemeris>>> orbits;
};


}

#endif

