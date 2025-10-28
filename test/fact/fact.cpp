#include <iostream>
#include <limits>
#include <string>

#include <utl_numeric.h>
#include <mth_util.h>

template <typename T>
void print_fact(T n, T d)
{

  using namespace mth_util;
  std::cout << '\n' << "factorial(" << n << ")    = " << factorial(n);
  std::cout << '\n';
  std::cout << '\n' << "factorial(" << n << "/"  << d <<
                                            ")  = " << factorial(n)/
                                                         factorial(d);
  std::cout << '\n' << "factorial(" << n << ", " << d <<
                                            ") = " << factorial(n, d);
}

/*
 * Quick test program for the factorial functions
 *
 * $ fact n d
 *
 * n!/d! is computed using int, long, and double
 */
int main(int argc, char* argv[])
{

  std::cout << "\nmachin<short>:        " << utl_numeric::machin<short>();
  std::cout << "\nmachin<int>:          " << utl_numeric::machin<int>();
  std::cout << "\nmachin<long>:         " << utl_numeric::machin<long>();
  std::cout << '\n';
  std::cout << "\nmachin<float>:        " << utl_numeric::machin<float>();
  std::cout << "\nmachin<double>:       " << utl_numeric::machin<double>();
  std::cout << "\nmachin<long double>:  " << utl_numeric::machin<double>();
  std::cout << '\n';
  std::cout << "\nlimits<float>:        " <<
      std::numeric_limits<float>::min() + std::numeric_limits<float>::min();
  std::cout << "\nlimits<double>:        " <<
      std::numeric_limits<double>::min() + std::numeric_limits<double>::min();
  std::cout << "\nlimits<long double>:  " <<
      std::numeric_limits<long double>::min() +
      std::numeric_limits<long double>::min();
  std::cout << '\n';

  if (argc != 3) {
    std::cerr << "\nProper use is:  " << argv[0] << " n d\n";
    return 0;
  }

  {
    std::cout << "\n--- int";
    int n {std::stoi(argv[1])};
    int d {std::stoi(argv[2])};
    print_fact(n, d);
  }
  std::cout << '\n';
  {
    std::cout << "\n--- long";
    long n {std::stol(argv[1])};
    long d {std::stol(argv[2])};
    print_fact(n, d);
  }
  std::cout << '\n';
  {
    std::cout << "\n--- double";
    double n {std::stod(argv[1])};
    double d {std::stod(argv[2])};
    print_fact(n, d);
  }

  std::cout << '\n';
}
