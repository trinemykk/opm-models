#include "padua.h"
#include <opm/material/densead/Evaluation.hpp>

using scalar_t = double;
constexpr int numPhases = 2;
constexpr int numComponents = 3;
constexpr int numDerivs = numPhases + numComponents;

template<typename Scalar>
void
padua_eval(const struct poly2d_t * const poly,
           const int num_pts,
           const Scalar * const x_vals,
           const Scalar * const y_vals,
           Scalar * const f_out) {

	// Padua interpolation is a two-level system: at the outer level, there is
	// a Chebyshev polynomial which can be evaluated together with the y-values
	// to get an approximation to the function value. each column in the primary
	// coefficient matrix represents a degree in the outer polynomial, and by
	// inputting this column as coefficients in a Chebyshev polynomial together
	// with the x-values, we get the coefficient for this degree in the outer
	// polynomial.

	// we use nested Clenshaw summation to evaluate the polynomials.

	// since we only use each coefficient once, when evaluating the term for a
	// particular degree in the outer polynomial, then we can generate the
	// coefficients on the fly as we need them, by nesting another recurrence
	// within the outer recurrence.

	// both x and y must be in the domain [-1, 1]. these are the loop-invariant
	// factors of the translation: (1-(-1))*(x-lo_x)/(up_x-lo_x)+(-1)
	const double span_x = 2. / (poly->up_x - poly->lo_x);
	const double span_y = 2. / (poly->up_y - poly->lo_y);

	// local alias
	const int n = poly->degree;

	for(int i = 0; i < num_pts; ++i) {
		// change of variables; xx is the alpha value of the inner recurrence,
		// and yy is the alpha value of the outer recurrence
		const Scalar xx = 2. * ((x_vals[i] - poly->lo_x) * span_x - 1.);
		const Scalar yy = 2. * ((y_vals[i] - poly->lo_y) * span_y - 1.);

		// the coefficients will be traversed from the right side of the matrix
		// up from the diagonal to the top. each coefficient matrix element is
		// only used once, so assuming that the coefficients are stored in this
		// order, we can traverse the matrix in fire-hose fashion with a single
		// pointer. this pointer has to be reinitialized for every point, here.
		const double * coeff = &poly->coeffs[0];

		// initialize the upper items of the recurrence
		Scalar b_j = 0.;
		Scalar b_jp1 = 0.;  // jp1 = j + 1

		// compute recurrence steps from the top, moving old values up
		for(int j = n; ; --j) {
			// --- begin inner recurrence for coefficient from x-value ---

			// g is the b series, but for the inner recurrence
			Scalar g_k = 0.;
			Scalar g_kp1 = 0.;  // kp1 = k + 1

			// often we see variants of this doing a double step to avoid
			// moving the temporary up (one variable then holds the even betas
			// and one variable holds the odd betas), but if we are not doing a
			// vectorized version, it is probably easier for the processor to do
			// a move than to try to accomodate the if statement afterwards that
			// handles the last odd recurrence.
			for(int k = n-j; k > 0; --k) {
				const Scalar prev_gk = g_k;  // temp. swap var.
				g_k = *(coeff++) + xx * g_k - g_kp1;
				g_kp1 = prev_gk;
			}

			// since coefficients are fetched directly out of a memory array,
			// we can do the last step explicitly here
			const Scalar c_j = *(coeff++) + 0.5 * xx * g_k - g_kp1;

			// --- end inner recurrence for coefficient from x-value ---

			// to avoid having to code the inner recurrence twice, once inside
			// the loop that does the outer recurrence step, and once for the
			// last step, we always do the calculation of the coefficient, but
			// move the end condition of the loop from the while statement at
			// the top into the middle of the loop.

			if(j != 0) {
				const Scalar prev_bj = b_j;  // temp. swap var.
				b_j = c_j + yy * b_j - b_jp1;
				b_jp1 = prev_bj;
			}
			else {
				// first coefficient only counts half
				f_out[i] = c_j + 0.5 * yy * b_j - b_jp1;

				// outer recurrence is finished; progress to next point
				break;
			}
		}
	}
}

// force the linker to instantiate a version for double-precision floats
template
void
padua_eval(const struct poly2d_t * const poly,
           const int num_pts,
           const scalar_t * const x_vals,
           const scalar_t * const y_vals,
           scalar_t * const f_out);

// force the linker to instantiate a version for automatic differentiation
template
void
padua_eval(const struct poly2d_t * const poly,
           const int num_pts,
           const Opm::DenseAd::Evaluation<scalar_t, numDerivs> * const x_vals,
           const Opm::DenseAd::Evaluation<scalar_t, numDerivs> * const y_vals,
           Opm::DenseAd::Evaluation<scalar_t, numDerivs> * const f_out);
