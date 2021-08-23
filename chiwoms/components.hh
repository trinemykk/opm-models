#ifndef COMPONENTS_HH
#define COMPONENTS_HH

#include "chiwoms.h"

#include <opm/material/IdealGas.hpp>
#include <opm/material/components/Component.hpp>
#include <opm/material/components/SimpleCO2.hpp>
#include <opm/material/components/CO2.hpp>
#include <opm/material/components/H2O.hpp>

namespace Opm {
/*!
 * \ingroup Components
 *
 * \brief A simple representation of linear octane
 *
 * \tparam Scalar The type used for scalar values
 */
template <class Scalar>
class Octane : public Opm::Component<Scalar, Octane<Scalar> >
{
public:
        /// Chemical name
        static const char* name() { return "C8"; }

        /// Molar mass in \f$\mathrm{[kg/mol]}\f$
        static Scalar molarMass() { return 0.11423; }

        /// Critical temperature in \f$\mathrm[K]}\f$
        static Scalar criticalTemperature() { return 568.7; }

        /// Critical pressure in \f$\mathrm[Pa]}\f$
        static Scalar criticalPressure() { return 2.49e6; }

        /// Acentric factor
        static Scalar acentricFactor() { return 0.398; }
};

template <class Scalar>
class NDekane : public Opm::Component<Scalar, NDekane<Scalar> >
{
public:
        /// Chemical name
        static const char* name() { return "C10"; }

        /// Molar mass in \f$\mathrm{[kg/mol]}\f$
        static Scalar molarMass() { return 0.1423; }

        /// Critical temperature in \f$\mathrm[K]}\f$
        static Scalar criticalTemperature() { return 617.7; }

        /// Critical pressure in \f$\mathrm[Pa]}\f$
        static Scalar criticalPressure() { return 2.103e6; }

        /// Acentric factor
        static Scalar acentricFactor() { return 0.4884; }

        // Critical volume
        static Scalar criticalVolume() {return 6.0976e-4; }
};

template <class Scalar>
class Methane : public Opm::Component<Scalar, Methane<Scalar> >
{
public:
        /// Chemical name
        static const char* name() { return "CH4"; }

        /// Molar mass in \f$\mathrm{[kg/mol]}\f$
        static Scalar molarMass() { return 0.0160; }

        /// Critical temperature in \f$\mathrm[K]}\f$
        static Scalar criticalTemperature() { return 190.5640; }

        /// Critical pressure in \f$\mathrm[Pa]}\f$
        static Scalar criticalPressure() { return 4.599e6; }

        /// Acentric factor
        static Scalar acentricFactor() { return 0.0114; }

        // Critical volume
        static Scalar criticalVolume() {return 9.8628e-5; }
};


template <class Scalar>
class ChiwomsCO2 : public Opm::SimpleCO2<Scalar>
{
public:
         /// Chemical name
        static const char* name() { return "CO2"; }

        /// Molar mass in \f$\mathrm{[kg/mol]}\f$
        static Scalar molarMass() { return 0.0440095; }

        /// Critical temperature in \f$\mathrm[K]}\f$
        static Scalar criticalTemperature() { return 304.1; }

        /// Critical pressure in \f$\mathrm[Pa]}\f$
        static Scalar criticalPressure() { return 7.38e6; }

        /// Acentric factor
        static Scalar acentricFactor() { return 0.225; }

        // Critical volume
        static Scalar criticalVolume() {return 9.4118e-5; }
};

template <class Scalar>
class ChiwomsBrine : public Opm::H2O<Scalar>
{
public:
        /// Molar mass in \f$\mathrm{[kg/mol]}\f$
        static Scalar molarMass() { return 0.0180158; }

        /// Acentric factor
        static Scalar acentricFactor() { return 0.344; }
};

struct EOS
{
        // template<typename LhsEval>
        // static LhsEval oleic_density(LhsEval T, LhsEval p, LhsEval x) {
        //         assert(T == (TEMPERATURE + 273.15));
        //         return 650;
        // }

        // template<typename LhsEval>
        // static LhsEval aqueous_density(LhsEval T, LhsEval p, LhsEval x) {
        //         assert(T == (TEMPERATURE + 273.15));
        //         return 1000;
        // }

        template<typename LhsEval>
        static LhsEval oleic_viscosity(LhsEval T, LhsEval p, LhsEval x) {
                assert(T == (TEMPERATURE + 273.15));
                return 5e-3;
        }

        template<typename LhsEval>
        static LhsEval aqueous_viscosity(LhsEval T, LhsEval p, LhsEval x) {
                assert(T == (TEMPERATURE + 273.15));
                return 1e-3;
        }

        template<typename LhsEval>
        static LhsEval oleic_enthalpy(LhsEval T, LhsEval p, LhsEval x) {
                assert(T == (TEMPERATURE + 273.15));
                return 0;
        }

        template<typename LhsEval>
        static LhsEval aqueous_enthalpy(LhsEval T, LhsEval p, LhsEval x) {
                assert(T == (TEMPERATURE + 273.15));
                return 0;
        }
};

} // namespace opm

#endif // COMPONENTS_HH
