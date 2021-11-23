#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <memory>
#include <vector>
#include <deque>
//#include <map>
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
 * Table of surfaces
 */
//std::map<std::string, std::unique_ptr<mth::Surface>> surface_table;

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

  eom_app::EomConfig cfg;
  std::vector<eom::OrbitDef> orbitDefs;

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
              /*
              std::shared_ptr<IShoot> sp_m14 =
                                      std::make_shared<M14>("make_shared");
              shooters_shared.emplace_back(sp_m14);

              void shoot_sp(std::vector<std::shared_ptr<IShoot>>& shooter_lst)
              {
                unsigned int n = static_cast<unsigned int>(shooter_lst.size());
                for (unsigned int ii=0;  ii<n; ++ii) {
                  shooter_lst[ii]->fire();
                }
              }
              */
              //std::shared_ptr<eom::Ephemeris> orbit =
              //              eom_app::parse_orbit(tokens, cfg, f2iSys);
              //orbits.push_back(orbit);
              orbitDefs.push_back(eom_app::parse_orbit_def(tokens, cfg));
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

  eom::JulianDate minJd = cfg.getStartTime();
  eom::JulianDate maxJd = cfg.getStopTime();
  for (auto& orbit : orbitDefs) {
    if (orbit.getEpoch() < minJd) {
      minJd = orbit.getEpoch();
    }
    if (maxJd < orbit.getEpoch()) {
      maxJd = orbit.getEpoch();
    }
  }
  auto f2iSys = std::make_shared<eom::EcfEciSys>(minJd, maxJd,
                                                 cfg.getEcfEciRate());
  std::cout << "\nMinJd " << minJd.to_str();
  std::cout << "\nMaxJd " << maxJd.to_str();

  std::vector<std::shared_ptr<eom::Ephemeris>> orbits;
  for (auto& orbit : orbitDefs) {
    orbits.emplace_back(eom::build_orbit(orbit, f2iSys));
  }

  using namespace utl_units;
  std::cout << "\nSize of orbits: " << orbits.size();
  if (orbits.size() > 0) {
    eom::print_ephemeris("test_i.e", cfg.getStartTime(), cfg.getStopTime(),
                         {60.0, phy_const::tu_per_sec},
                         eom::EphemFrame::eci, orbits[0]);
  }

  f2iSys->print(std::cout);


/*
  auto jdTmp = cfg.getStartTime() + 0.5*cfg.getEcfEciRate().getDays();
  auto f2i = f2iSys->getEcfEciData(jdTmp);
  std::cout << "\nJD2000: " << f2i.mjd2000;
  std::cout << '\n' << f2i.bpn.w();
  std::cout << '\n' << f2i.bpn.vec();
  std::cout << '\n' << f2i.pm.w();
  std::cout << '\n' << f2i.pm.vec();

  Eigen::Matrix<double, 6, 1> eci0;
  eci0(0,0) = -5552.0_km;
  eci0(1,0) = -2563.0_km;
  eci0(2,0) =  3258.0_km;
  eci0(3,0) =     2.149_kms;
  eci0(4,0) =    -7.539_kms;
  eci0(5,0) =    -2.186_kms;
  Eigen::Matrix<double, 6, 1> ecf = f2iSys->eci2ecf(cfg.getStartTime(),
                                                    eci0.block<3,1>(0,0),
                                                    eci0.block<3,1>(3,0));
  Eigen::Matrix<double, 6, 1> eci = f2iSys->ecf2eci(cfg.getStartTime(),
                                                    ecf.block<3,1>(0,0),
                                                    ecf.block<3,1>(3,0));

*/

  //eom::Kepler orbit(cfg.getStartTime(), eci, f2iSys);


/*
  std::cout << '\n';
  std::cout.precision(17);
  std::cout << "\nJD Start: " << cfg.getStartTime().getJdHigh() +
                                 cfg.getStartTime().getJdLow();
  std::cout << '\n' << phy_const::km_per_du*eci0.block<3,1>(0,0);
  std::cout << '\n' << phy_const::km_per_du*phy_const::tu_per_sec*
                       eci0.block<3,1>(3,0);
  std::cout << '\n';
  std::cout << '\n' << phy_const::km_per_du*ecf.block<3,1>(0,0);
  std::cout << '\n' << phy_const::km_per_du*phy_const::tu_per_sec*
                       ecf.block<3,1>(3,0);
  std::cout << '\n';
  std::cout << '\n' << phy_const::km_per_du*eci.block<3,1>(0,0);
  std::cout << '\n' << phy_const::km_per_du*phy_const::tu_per_sec*
                       eci.block<3,1>(3,0);
*/

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
