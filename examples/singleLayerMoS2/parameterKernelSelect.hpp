#ifndef MOS2_PARAMETER_KERNEL_SELECT_HPP
#define MOS2_PARAMETER_KERNEL_SELECT_HPP
// Compile-time choice of the monolayer-MoS2 transport kernel used by the validation
// drivers: default = single-valley K kernel (parameterKaasbjerg.hpp, first-principles
// couplings); -DMOS2_TWO_VALLEY = two-valley K+Q kernel v3 (parameterTwoValley.hpp);
// -DMOS2_TWO_VALLEY_V2 = the frozen v2 kernel of cluster array 196933 (parameterTwoValleyV2.hpp).
// Both expose the same add*() functions under the alias MoS2Kernel.
#if defined(MOS2_TWO_VALLEY_V2)
#include "parameterTwoValleyV2.hpp"
namespace MoS2Kernel = MoS2TwoValleyV2;
#elif defined(MOS2_TWO_VALLEY)
#include "parameterTwoValley.hpp"
namespace MoS2Kernel = MoS2TwoValley;
#else
#include "parameterKaasbjerg.hpp"
namespace MoS2Kernel = MoS2Kaasbjerg;
#endif
#endif
