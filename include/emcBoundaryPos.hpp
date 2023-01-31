#ifndef EMC_BOUNDARY_POS_HPP
#define EMC_BOUNDARY_POS_HPP

#include <emcUtil.hpp>

/// enum that links all different positions for boundaries to
/// the index of grid belonging to that position
enum struct emcBoundaryPos : SizeType {
  XMIN = 0, /**< index for boundary at x = 0 */
  XMAX,     /**< index for boundary at x = deviceMaxPos[x] */
  YMIN,     /**< index for boundary at y = 0 */
  YMAX,     /**< index for boundary at y = deviceMaxPos[y] */
  ZMIN,     /**< index for boundary at z = 0 */
  ZMAX,     /**< index for boundary at z = deviceMaxPos[z] */
  INVALID
};

#endif // EMC_BOUNDARY_POS_HPP