#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <list>
#include <vector>

#include "FirSmoother.h"
#include "Path.h"

int main() {
        // 构造一个一维路径，使用 Kunz TOTG 生成时间律后再进行平滑。
        std::list<Eigen::VectorXd> waypoints;
        Eigen::VectorXd waypoint(1);
        waypoint << 0.0;
        waypoints.push_back(waypoint);
        waypoint << 1.0;
        waypoints.push_back(waypoint);
        waypoint << 2.0;
        waypoints.push_back(waypoint);

        Eigen::VectorXd maxVel(1);
        maxVel << 1.0;
        Eigen::VectorXd maxAcc(1);
        maxAcc << 2.0;

        Trajectory trajectory(Path(waypoints, 0.0), maxVel, maxAcc, 0.001);
        assert(trajectory.isValid());

        FirSmoother smoother;
        FirSmoothingConfig config;
        config.samplePeriod = 0.01;
        config.primaryKernelLength = 5;
        config.jerkLimit = 50.0;
        config.allowTimeScale = true;

        const std::vector<Trajectory::Sample> original = trajectory.sampleTimeLaw(config.samplePeriod);
        SmoothedTimeLaw smoothed = smoother.smooth(trajectory, config);

        const double boundaryTolerance = 1e-3;
        assert(std::fabs(smoothed.position.front() - original.front().position) < boundaryTolerance);
        assert(std::fabs(smoothed.position.back() - original.back().position) < boundaryTolerance);
        assert(std::fabs(smoothed.velocity.front() - original.front().velocity) < boundaryTolerance);
        assert(std::fabs(smoothed.velocity.back() - original.back().velocity) < boundaryTolerance);

        double originalPeakJerk = 0.0;
        double smoothedPeakJerk = 0.0;
        for(std::size_t i = 1; i < original.size(); ++i) {
                const double jerk = std::fabs((original[i].acceleration - original[i - 1].acceleration) / config.samplePeriod);
                originalPeakJerk = std::max(originalPeakJerk, jerk);
        }
        for(std::size_t i = 1; i < smoothed.acceleration.size(); ++i) {
                const double jerk = std::fabs((smoothed.acceleration[i] - smoothed.acceleration[i - 1]) / config.samplePeriod);
                smoothedPeakJerk = std::max(smoothedPeakJerk, jerk);
        }
        assert(smoothedPeakJerk <= originalPeakJerk + 1e-6);

        std::cout << "Jerk before: " << originalPeakJerk << " after: " << smoothedPeakJerk << std::endl;
        std::cout << "All smoother checks passed.\n";
        return 0;
}

