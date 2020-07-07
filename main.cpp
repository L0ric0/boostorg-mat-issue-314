#include <boost/math/interpolators/barycentric_rational.hpp>
#include <boost/range/adaptors.hpp>
#include <boost/math/tools/roots.hpp>
#include <boost/units/systems/si.hpp>
#include <boost/units/io.hpp>

#include <map>
#include <iostream>

using namespace boost::units;

int main()
{
    // The lithium potential is given in Kohn's paper, Table I.
    // (We could equally easily use an unordered_map, a list of tuples or pairs, or a 2-dimensional array).
    std::map<quantity<si::energy, double>, quantity<si::electric_charge, double>> r;

    r[0.02 *si::joule] = 5.727 * si::coulomb;
    r[0.04 *si::joule] = 5.544 * si::coulomb;
    r[0.06 *si::joule] = 5.450 * si::coulomb;
    r[0.08 *si::joule] = 5.351 * si::coulomb;
    r[0.10 *si::joule] = 5.253 * si::coulomb;
    r[0.12 *si::joule] = 5.157 * si::coulomb;
    r[0.14 *si::joule] = 5.058 * si::coulomb;
    r[0.16 *si::joule] = 4.960 * si::coulomb;
    r[0.18 *si::joule] = 4.862 * si::coulomb;
    r[0.20 *si::joule] = 4.762 * si::coulomb;
    r[0.24 *si::joule] = 4.563 * si::coulomb;
    r[0.28 *si::joule] = 4.360 * si::coulomb;
    r[0.32 *si::joule] = 4.1584 * si::coulomb;
    r[0.36 *si::joule] = 3.9463 * si::coulomb;
    r[0.40 *si::joule] = 3.7360 * si::coulomb;
    r[0.44 *si::joule] = 3.5429 * si::coulomb;
    r[0.48 *si::joule] = 3.3797 * si::coulomb;
    r[0.52 *si::joule] = 3.2417 * si::coulomb;
    r[0.56 *si::joule] = 3.1209 * si::coulomb;
    r[0.60 *si::joule] = 3.0138 * si::coulomb;
    r[0.68 *si::joule] = 2.8342 * si::coulomb;
    r[0.76 *si::joule] = 2.6881 * si::coulomb;
    r[0.84 *si::joule] = 2.5662 * si::coulomb;
    r[0.92 *si::joule] = 2.4242 * si::coulomb;
    r[1.00 *si::joule] = 2.3766 * si::coulomb;
    r[1.08 *si::joule] = 2.3058 * si::coulomb;
    r[1.16 *si::joule] = 2.2458 * si::coulomb;
    r[1.24 *si::joule] = 2.2035 * si::coulomb;
    r[1.32 *si::joule] = 2.1661 * si::coulomb;
    r[1.40 *si::joule] = 2.1350 * si::coulomb;
    r[1.48 *si::joule] = 2.1090 * si::coulomb;
    r[1.64 *si::joule] = 2.0697 * si::coulomb;
    r[1.80 *si::joule] = 2.0466 * si::coulomb;
    r[1.96 *si::joule] = 2.0325 * si::coulomb;
    r[2.12 *si::joule] = 2.0288 * si::coulomb;
    r[2.28 *si::joule] = 2.0292 * si::coulomb;
    r[2.44 *si::joule] = 2.0228 * si::coulomb;
    r[2.60 *si::joule] = 2.0124 * si::coulomb;
    r[2.76 *si::joule] = 2.0065 * si::coulomb;
    r[2.92 *si::joule] = 2.0031 * si::coulomb;
    r[3.08 *si::joule] = 2.0015 * si::coulomb;
    r[3.24 *si::joule] = 2.0008 * si::coulomb;
    r[3.40 *si::joule] = 2.0004 * si::coulomb;
    r[3.56 *si::joule] = 2.0002 * si::coulomb;
    r[3.72 *si::joule] = 2.0001 * si::coulomb;

    // Let's discover the absissa that will generate a potential of exactly 3.0,
    // start by creating 2 ranges for the x and y values:
    auto x_range = boost::adaptors::keys(r);
    auto y_range = boost::adaptors::values(r);
    boost::math::barycentric_rational<quantity<si::energy, double>, quantity<si::electric_charge, double>> b(x_range.begin(), x_range.end(), y_range.begin());
    //
    // We'll use a lambda expression to provide the functor to our root finder, since we want
    // the abscissa value that yields 3, not zero.  We pass the functor b by value to the
    // lambda expression since barycentric_rational is trivial to copy.
    // Here we're using simple bisection to find the root:
    boost::uintmax_t iterations = (std::numeric_limits<boost::uintmax_t>::max)();
    double abscissa_3 = boost::math::tools::bisect([=](quantity<si::energy, double> x) { return b(x) - 3.* si::joule; }, 0.44 *si::joule, 1.24 *si::joule, boost::math::tools::eps_tolerance<double>(), iterations).first;
    std::cout << "Abscissa value that yields a potential of 3 = " << abscissa_3 << std::endl;
    std::cout << "Root was found in " << iterations << " iterations." << std::endl;
    //
    // However, we have a more efficient root finding algorithm than simple bisection:
    iterations = (std::numeric_limits<boost::uintmax_t>::max)();
    abscissa_3 = boost::math::tools::bracket_and_solve_root([=](quantity<si::energy, double> x) { return b(x) - 3. *si::joule; }, 0.6 *si::joule, 1.2 *si::joule, false, boost::math::tools::eps_tolerance<double>(), iterations).first;
    std::cout << "Abscissa value that yields a potential of 3 = " << abscissa_3 << std::endl;
    std::cout << "Root was found in " << iterations << " iterations." << std::endl;
}
