#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <memory>
#include <vector>
#include <deque>
#include <unordered_map>
#include <stdexcept>


#include <Eigen/Dense>
#include <phy_const.h>
#include <utl_units.h>

#include <eom_config.h>
#include <eom_parse.h>
#include <astro_orbit_def.h>
#include <astro_ephemeris.h>
#include <astro_ecfeci_sys.h>
#include <astro_build.h>
#include <astro_print.h>

/**
 * Parses an input file and passes each line to the parser.
 */
int main(int argc, char* argv[])
{

    // Check for filename
  if (argc != 2) {
    std::cerr << "\nProper use is:  " << argv[0] << " <input_file_name>\n";
    return 0;
  }
    // Try to open for input
  std::ifstream ifs(argv[1]);
  if (!ifs.is_open()) {
    std::cerr << "\nError opening " << argv[1] << "\n";
    return 0;
  }
  std::cout << "\nOpened " << argv[1] << '\n';

    // Parse input file and generate the simulation configuration
    // parameters along with modeling component definitions (that will
    // be used to create the actual modeling components)
  eom_app::EomConfig cfg;
  std::vector<eom::OrbitDef> orbit_defs;
    // Read each line and pass to parser while tracking line number
  int line_number {0};
  std::string input_line;
  bool parse_tokens = false;
  bool input_error = false;
  std::deque<std::string> tokens;
  while (std::getline(ifs,input_line)) {
    line_number++;
    std::istringstream iss(input_line);
    std::string token;
    while (iss >> token  &&  !input_error) {
      if (token.front() == '#') {
        break;
      } else {
        if (token.back() == ';') {
          parse_tokens = true;
          if (token.size() > 1) {
            token.pop_back();
            tokens.push_back(token);
          }
        } else {
          tokens.push_back(token);
        }
        if (parse_tokens) {
          input_error = true;
          parse_tokens = false;
          if (tokens.size() > 0) {
            auto make = tokens[0];
            tokens.pop_front();
            // Start Input Types
            if (make == "SimStart") {
              cfg.setStartTime(tokens);
              input_error = !cfg.isValid();
            } else if (make == "SimDuration") {
              cfg.setDuration(tokens);
              input_error = !cfg.isValid();
            } else if (make == "LeapSeconds") {
              cfg.setLeapSeconds(tokens);
              input_error = !cfg.isValid();
            } else if (make == "EcfEciRate") {
              cfg.setEcfEciRate(tokens);
              input_error = !cfg.isValid();
            } else if (make == "end") {
              input_error = false;
              ifs.seekg(0, std::ios::end);
            } else if (make == "ToKilometers") {
              cfg.setToKilometers(tokens);
              input_error = !cfg.isValid();
            } else if (make == "ToSeconds") {
              cfg.setToSeconds(tokens);
              input_error = !cfg.isValid();
            } else if (make == "Orbit") {
              orbit_defs.push_back(eom_app::parse_orbit_def(tokens, cfg));
              input_error = false;
            } else if (make == "Command") {
              //commands.push_back(eom_app::parse_command(tokens);
              input_error = false;
            }
            // End Input Types
          }
          tokens.clear();
        }
      }
    }
    if (input_error) {
      std::cout << "\nError on line: " << line_number;
      std::cout << '\n' << cfg.getError();
      std::cout << '\n';
      break;
    }
  }
  ifs.close();
  if (input_error) {
    return 0;
  }
  cfg.print(std::cout);

    // Determine time span that must be supported by the simulation
    // based on the input scenario time and orbit epoch times.
  eom::JulianDate minJd = cfg.getStartTime();
  eom::JulianDate maxJd = cfg.getStopTime();
  for (const auto& orbit : orbit_defs) {
    if (orbit.getEpoch() < minJd) {
      minJd = orbit.getEpoch();
    }
    if (maxJd < orbit.getEpoch()) {
      maxJd = orbit.getEpoch();
    }
  }

  //
  // Create resources and modeling components
  //
  // Integrity of the application relies on the lists containing
  // models within the environment to remain unchanged during the
  // duration of the simulation once initialization has completed.
  //

    // Ecf to Eci transformation service, pass as
    // const std::shared_ptr<const EcfEciSys>&
    // for maximum safety
  auto f2iSys = std::make_shared<eom::EcfEciSys>(minJd, maxJd,
                                                 cfg.getEcfEciRate());
    // Internal numeric orbit ID is the location of the orbit in the
    // ephemeris vector.  Generate orbits, and note Name/ID association.
  std::unordered_map<std::string, int> orbit_ids;
  std::vector<std::shared_ptr<eom::Ephemeris>> orbits;
  int ii {0};
  for (const auto& orbit : orbit_defs) {
    orbit_ids[orbit.getOrbitName()] = ii;
    orbits.emplace_back(eom::build_orbit(orbit, f2iSys));
    ii++;
  }

  using namespace utl_units;
  std::cout << "\nSize of orbits: " << orbits.size();
  if (orbits.size() > 0) {
    eom::print_ephemeris("test_i.e", cfg.getStartTime(), cfg.getStopTime(),
                         {60.0, phy_const::tu_per_sec},
                         eom::EphemFrame::eci, orbits[0]);
  }

  f2iSys->print(std::cout);



  std::cout << "\n\n";

}

  //if (surface_table.at("SPH")->getType() == mth::SurfaceType::SPHERE) {
  //  std::cout << "\nYes, a sphere";
  //}

              //try {
              //  std::string name = tokens[0];
              //  mth::Sphere* sptr = 
              //    dynamic_cast<mth::Sphere*>(surface_table.at(name).get());
              //  std::cout << " A sphere of size " << sptr->getRadius();
              //  plt_sphere(sptr);
              //  input_error = false;
              //} catch (std::out_of_range const& exc) {
              //  std::cout << exc.what() << '\n';
              //  input_error = false;
              //}
