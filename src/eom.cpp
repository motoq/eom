#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <memory>
//#include <vector>
#include <deque>
//#include <map>
#include <stdexcept>



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
    // Read each line and pass to parser while tracking line number
  int line_number {0};
  std::string input_line;
  bool parse_tokens = false;
  bool no_error = true;
  std::deque<std::string> tokens;
  while (std::getline(ifs,input_line)) {
    line_number++;
    std::istringstream iss(input_line);
    std::string token;
    while (iss >> token  &&  no_error) {
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
          no_error = false;
          parse_tokens = false;
          if (tokens.size() > 0) {
            std::string make = tokens[0];
            tokens.pop_front();
            if (make == "SimStart"  &&  tokens.size() > 1) {
              //std::string name = tokens[1];
              //tokens.erase(tokens.begin()+1);
              //surface_table.emplace(std::make_pair(name,
              //                      cog_make_surface(tokens)));
              std::cout << "\nSimStart model: " << tokens[0];
              no_error = true;
            } else if (make == "SimDuration"  &&  tokens.size() > 1) {
              /*
              try {
                std::string name = tokens[0];
                mth::Sphere* sptr = 
                  dynamic_cast<mth::Sphere*>(surface_table.at(name).get());
                std::cout << " A sphere of size " << sptr->getRadius();
                plt_sphere(sptr);
                no_error = true;
              } catch (std::out_of_range const& exc) {
                std::cout << exc.what() << '\n';
                no_error = false;
              }
              */
              std::cout << "\nSimDuration model: " << tokens[0];
              no_error = true;
            }
          }
          tokens.clear();
        }
      }
    }
    if (!no_error) {
      std::cout << "\nError on line: " << line_number;
      break;
    }
  }
  ifs.close();

  //if (surface_table.at("SPH")->getType() == mth::SurfaceType::SPHERE) {
  //  std::cout << "\nYes, a sphere";
  //}
  

  std::cout << "\n\n";

}


