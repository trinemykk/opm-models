#ifndef __CHIWOMS_H__
#define __CHIWOMS_H__

// values are from the paper "Field-scale implications of density-driven
// convection in CO2-EOR reservoirs", to be presented at the Fifth CO2
// Geological Storage Workshop, at 21–23 November 2018, in Utrecht,
// The Netherlands.

constexpr double TEMPERATURE = 80;   /* degree Celsius */
const double MIN_PRES = 180;     /* bars */
const double MAX_PRES = 220;     /* bars */
constexpr double SIM_TIME = 1;      /* days */
constexpr double X_SIZE = 100;        /* centimeter */
constexpr double Y_SIZE = 100;        /* centimeter */
const unsigned NX = 5;         /* number of cells horizontally */
const unsigned NY = 5;         /* number of cells vertically */
const double POROSITY = 0.1;     /* non-dimensional */
const double PERMEABILITY = 100; /* milli-Darcy */
const double DIFFUSIVITY = 1e-9; /* square meter per second */
const double PERTUBATION = 1e-3; /* of concentration in boundary layer */
const unsigned ZONE = 1;         /* blocks on each side in pertubation zone */
constexpr double WAVE_LENGTH = 0.1;  /* critical wavelength fraction of height */

/* "random" fields will be equal as long as this is set the same */
const double SEED = 5163166242092481088;

#endif /* __CHIWOMS_H__ */
