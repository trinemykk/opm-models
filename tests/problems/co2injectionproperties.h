#ifndef __CO2INJECTIONPARAMETERS_H__
#define __CO2INJECTIONPARAMETERS_H__

// values are from the paper "Field-scale implications of density-driven
// convection in CO2-EOR reservoirs", to be presented at the Fifth CO2
// Geological Storage Workshop, at 21–23 November 2018, in Utrecht,
// The Netherlands.

constexpr double TEMPERATURE = 80;   /* degree Celsius */
constexpr double GRAVITYFACTOR = 1;   /* fraction og gravity */
constexpr double MIN_PRES = 10;     /* bars */
const double MAX_PRES = 220;     /* bars */
constexpr double SIM_TIME = 1;      /* days */
constexpr double Y_SIZE = 1.0;        /* meter */
constexpr double X_SIZE = 10.0;        /* meter */
constexpr double Z_SIZE = 5.0;        /* meter */
const unsigned NX = 10;         /* number of cells x-dir */
const unsigned NY = 1;         /* number of cells y-dir */
const unsigned NZ = 5;         /* number of cells z-dir */
const double POROSITY = 0.2;     /* non-dimensional */
const double PERMEABILITY = 100; /* milli-Darcy */
const double DIFFUSIVITY = 1e-9; /* square meter per second */
const double MFCOMP0 = 0.9;
const double MFCOMP1 = 0.1;
const double MFCOMP2 = 0.0;
constexpr double INFLOW_RATE = -1e-4; /* unit kg/s ? */

/* "random" fields will be equal as long as this is set the same */
const double SEED = 5163166242092481088;

#endif /* __CO2INJECTIONPARAMETERS_H__ */