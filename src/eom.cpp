#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <memory>
//#include <vector>
#include <deque>
//#include <map>
#include <stdexcept>

#include <eom_config.h>
#include <astro_ecfeci_sys.h>


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

  auto f2iSys = std::make_shared<eom::EcfEciSys>(cfg.getStartTime(),
                                                 cfg.getStopTime(),
                                                 cfg.getEcfEciRate());
  auto jdTmp = cfg.getStartTime() + 0.5*cfg.getEcfEciRate().getDays();
  auto f2i = f2iSys->getEcfEciData(jdTmp);
  std::cout << "\nJD2000: " << f2i.jd2000;

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
