/*
 * Copyright 2022 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <iostream>
#include <fstream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>

#include <astro_math.h>

/**
 * Structure of spherical harmonic terms.
 */
struct egm_rec {
  int degree {0};            ///< 'n' subscriipt
  int order {0};             ///< 'm' subscript
  double cnm {0.0};          ///< Cosine terms
  double snm {0.0};          ///< Sine terms

  /**
   * Initialize
   *
   * @param  n  Degree
   * @param  m  Order (<= Degree to be valid)
   * @param  c  Cosine terms
   * @param  s  Sine terms
   */
  egm_rec(int n, int m, double c, double s) :
           degree {n}, order {m}, cnm {c}, snm {s}
  {
  }

    // Sorting operator
  bool operator<(const egm_rec& rec) const
  {
    if (order ==  rec.order) {
      return degree < rec.degree;
    } else {
      return order < rec.order  ; 
    }
  }
};


/**
 * Reads an EGM96, EGM2008, or similar formatted file where each line
 * begins with space separated degree, order, cosine, and sine
 * coefficients (trailing uncertainty values are ignored).  The
 * coefficients are assumed to be normalized.
 *
 * For output, the spherical harmonic coefficients are sorted such that
 * zonal terms are first (order == 0).  Remaining terms continue to sort
 * such that all coefficients of the same order are grouped together.
 * They are unormalized and printed to either cout by rows, or written
 * to a file as constexpr std::array<> types.
 *
 * For both output options, the 'D' exponent character is converted
 * to a 'e'.
 *
 * @author  Kurt Motekew
 * @date    2022/07/10  Initial
 */
int main(int argc, char* argv[])
{
  if (argc < 4  ||  argc > 5) {
    std::cerr << "\nProper use is:  " << argv[0] <<
                 " <egm_filename> <deg> <order> <out_filename>\n";
    return 0;
  }

    // Check for valid degree and order
  int max_degree  {std::stoi(argv[2])};
  int max_order {std::stoi(argv[3])};
  if (max_order > max_degree) {
    std::cerr << "Order must be <= deg: Degree " << max_degree <<
                                   "   Order = " << max_order << '\n';
    return 0;
  }

    // Open/validate input filename
  std::string ifname = argv[1];
  std::ifstream fin(ifname);
  if (!fin.is_open()) {
    std::cerr << "EopSys::EopSys Can't open " << ifname << '\n';
    return 0;
  }

  int degree;
  int order;
  double cnm;
  double snm;
  std::string input_line;
  std::vector<egm_rec> egm_data;
    // Read required number of values and unormalize.
  while (std::getline(fin, input_line)) {
    std::replace(input_line.begin(), input_line.end(), 'D', 'e');
    std::stringstream ss(input_line);
    if (!(ss >> degree >> order >> cnm >> snm)) {
      std::cerr << "Error parsing: " << input_line;
      return 0;
    }
      // Use degree to mark end of file parsing
    if (degree > max_degree) {
      break;
    }
      // Only include if order is also within desired range
    if (order <= max_order) {
      auto norm_fact = astro_math::kaula_norm(static_cast<double>(degree),
                                              static_cast<double>(order));
      cnm /= norm_fact;
      snm /= norm_fact;
      egm_data.emplace_back(degree, order, cnm, snm);
    }
  }
  fin.close();

  std::sort(egm_data.begin(), egm_data.end());

    // Output to C++ compatible file or print to stdout
  if (argc == 5) {
    std::string ofname = argv[4];
    std::ofstream fout(argv[4]);
    if (!fout.is_open()) {
      std::cerr << "Error opening : " << ofname << " for output";
      return 0;
    }
      // Degree index
    fout << "constexpr int degree {" << max_degree << "};";
    fout << "\nconstexpr int order {" << max_order << "};";
    fout << "\nconstexpr int nc {" << egm_data.size() << "};";
    fout << "\nconstexpr std::array<int, nc> xn = \n  {";
    for (const auto& rec : egm_data) {
      fout << '\n';
      fout.width(5);
      fout << rec.degree;
      fout << ",";
    }
    fout.seekp(-1, std::ios_base::cur);
      // Order index
    fout << "\n  };";
    fout << "\nconstexpr std::array<int, nc> xm = \n  {";
    for (const auto& rec : egm_data) {
      fout << '\n';
      fout.width(5);
      fout << rec.order;
      fout << ",";
    }
    fout.seekp(-1, std::ios_base::cur);
      // Cosine terms
    fout << "\n  };";
    fout.precision(15);
    fout << std::scientific;
    fout << "\nconstexpr std::array<double, nc> cnm = \n  {";
    for (const auto& rec : egm_data) {
      fout << '\n';
      fout.width(25);
      fout << rec.cnm;
      fout << ",";
    }
    fout.seekp(-1, std::ios_base::cur);
      // Sine terms
    fout << "\n  };";
    fout << "\nconstexpr std::array<double, nc> snm = \n  {";
    for (const auto& rec : egm_data) {
      fout << '\n';
      fout.width(25);
      fout << rec.snm;
      fout << ",";
    }
    fout.seekp(-1, std::ios_base::cur);
    fout << "\n  };";
    fout << '\n';
    fout.close();
  } else {
    std::cout << '\n';
    for (const auto& rec : egm_data) {
      std::cout << '\n';
      std::cout.width(5);
      std::cout << rec.degree;
      std::cout.width(5);
      std::cout << rec.order;
      std::cout.precision(15);
      std:: cout << std::scientific;
      std::cout.width(25);
      std::cout << rec.cnm;
      std::cout.width(25);
      std::cout << rec.snm;
    }
  }

  std::cout << '\n';
}
