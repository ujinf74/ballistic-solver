#pragma once

// Umbrella header: modular bs/ components in dependency order.
#include "bs/vec3.hpp"
#include "bs/math_utils.hpp"
#include "bs/params.hpp"
#include "bs/integrator.hpp"
#include "bs/target.hpp"
#include "bs/closest_approach.hpp"
#include "bs/vacuum_lead.hpp"
#include "bs/residual.hpp"
#include "bs/lm.hpp"
#include "bs/solve.hpp"
#include "bs/predictor.hpp"
// Note: bs/kalman.hpp (TargetTracker, a position-only estimator) is intentionally
// NOT included here. It is an optional add-on, not part of the core solver; include
// it explicitly with `#include "bs/kalman.hpp"` when you want the tracker.
