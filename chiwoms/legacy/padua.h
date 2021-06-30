#ifndef __PADUA_H__
#define __PADUA_H__

// the coefficients are just a data dump of calculated values. there is no point
// in initializing Scalar values for them; instead we rely on the Scalar type to
// be able to do multiplication and comparison with double precision floating
// point constants.

#ifdef __cplusplus
extern "C" {
#endif // __cplusplus

/**
 * Description of a two-dimensional Padua interpolation polynomial.
 */
struct poly2d_t {
	/// Degree of the highest term in the polynomial.
	int degree;

	/// Range in parameter space in which the polynomial is valid.
	double lo_x;
	double up_x;
	double lo_y;
	double up_y;

	/// Coefficient matrix for the bivariate polynomial. This array is
	/// organized in a peculiar way for optimal performance: only the upper
	/// left triangle of the matrix is stored, in a column-major fashion,
	/// and it is furhermore stored from the right, with the highest indexed
	/// column, from the diagonal and up, i.e. both the column and the row
	/// index decreases when traversing. There are (n+1)*(n+2)/2 items in
	/// this matrix, where n is the highest degree in the polynomial.
	double coeffs[];
};

#ifdef __cplusplus
} // extern "C"

// this function is now only available to C++ code. we need the header to work
// for plain old C as well, since the coefficient listing is printed in that
// language.

/**
 * Evaluate the interpolation polynomial to find function value approximations
 * in a set of points.
 *
 * @param[in] poly
 * Description of the polynomial to use for interpolation.
 *
 * @param[in] num_pts
 * Number of points to evaluate.
 *
 * @param[in] x_vals
 * First axis coordinate for evaluation points. This must be between the
 * lower_x and upper_x fields in the polynomial description. This argument
 * must point to an array with num_pts items.
 *
 * @param[in] y_vals
 * Second axis coordinate for evaluation points. This must be between the
 * lower_y and upper_y fields in the polynomial description. This argument
 * must point to an array with num_pts items.
 *
 * @param[out] f_out
 * Approximated function value for each of the points specified by x_vals
 * and y_vals. This argument must point to a vector that can receive num_pts
 * values. The previous content of the array will be overwritten.
 */
template<typename Scalar>
void
padua_eval(const struct poly2d_t * const poly,
           const int num_pts,
           const Scalar * const x_vals,
           const Scalar * const y_vals,
           Scalar * const f_out);
#endif // __cpluscplus

#endif
