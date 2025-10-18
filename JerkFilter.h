#pragma once

#include <vector>

#include "Trajectory.h"

/**
 * @brief Container for the diagnostics emitted by the jerk smoothing pipeline.
 *
 * The filter keeps track of the reconstructed phase trajectory samples as well as a
 * handful of scalar metrics that make it easy to validate the jerk reduction
 * claimed by the guide in TOTG_jerk_smoothing_guide.md.
 */
struct JerkFilterResult
{
        /// Resampled phase trajectory (positions/velocities/accelerations over time)
        /// after jerk smoothing has been applied.
        std::vector<Trajectory::PhaseSample> filteredSamples;

        /// Largest jerk encountered in the original trajectory discretization.
        double maxJerkOriginal = 0.0;

        /// Largest jerk encountered after filtering.
        double maxJerkFiltered = 0.0;

        /// Difference between filtered and original terminal path velocity.
        double endpointVelocityError = 0.0;

        /// Difference between filtered and original terminal path position.
        double endpointPositionError = 0.0;
};

/**
 * @brief Apply the jerk smoothing recipe documented in TOTG_jerk_smoothing_guide.md.
 *
 * The procedure mirrors the written guide: mirror-padded binomial filtering reduces
 * high-frequency acceleration content, a bi-directional jerk clamp enforces the
 * requested jerk limit, and a low-order moment correction restores the position and
 * velocity boundary conditions of the original TOTG solution.  The final pass
 * re-integrates the corrected acceleration profile so that downstream code can reuse
 * the trajectory samples without modification.
 *
 * @param samples    Uniform time samples of the original TOTG phase trajectory.
 * @param jerkLimit  Maximum |jerk| (path-acceleration change per second) to enforce.
 *
 * @return Diagnostic information and the smoothed trajectory samples.
 */
JerkFilterResult smoothTrajectoryWithJerkLimit(const std::vector<Trajectory::PhaseSample> &samples, double jerkLimit);

