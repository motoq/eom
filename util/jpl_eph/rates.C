// ROOT script prints angular velocities of the moon around the earth
// and of the planets around the sun
//
{
  double leo_rev_per_day {14.0};
  auto leo_deg_per_min = 360.0*leo_rev_per_day/1440.0;

  cout << '\n' << "LEO Satellite" <<
          "\n  " << leo_rev_per_day << " rev/day = " <<
                    leo_deg_per_min << " deg/min";

  double moon_rev_per_day {1.0/27.32};
  auto moon_deg_per_min = 360.0*moon_rev_per_day/1440.0;
  auto moon_min_per_deg = 1.0/moon_deg_per_min;
  auto moon_day_per_deg = moon_min_per_deg/1440.0;

  const int nplanets {9};
  array<string, nplanets> planet_names = { "Mercury",
                                           "Venus",
                                           "Earth",
                                           "Mars",
                                           "Jupiter",
                                           "Saturn",
                                           "Uranus",
                                           "Neptune",
                                           "Pluto" };
    // Sidereal Orbital Periods
    // <https://ssd.jpl.nasa.gov/planets/phys_par.html>
  array<double, nplanets> years_per_rev = { 0.2408467, 
                                            0.61519726,
                                            1.0000174,
                                            1.8808476,
                                            11.862615,
                                            29.447498,
                                            84.016846,
                                            164.79132,
                                            247.92065 };

  cout << '\n';
  cout << '\n' << "Moon" << "\t\t" << moon_day_per_deg << " day/deg";
  for (int ii=0; ii<nplanets; ++ii) {
    cout << '\n' << planet_names[ii] << "\t\t" <<
                    365.25*years_per_rev[ii]/360.0 << " day/deg";
  }

  cout << "\n\n";
}
