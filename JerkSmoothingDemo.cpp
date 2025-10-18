#include <cmath>
#include <iomanip>
#include <iostream>
#include <list>
#include <vector>

#include <Eigen/Core>

#include "JerkFilter.h"
#include "Path.h"
#include "Trajectory.h"

// This executable demonstrates the jerk-limited smoothing procedure described in
// TOTG_jerk_smoothing_guide.md by generating a time-optimal trajectory, sampling the
// phase profile, filtering it, and printing a brief before/after summary.

namespace
{

std::list<Eigen::VectorXd> makeDemoWaypoints()
{
        std::list<Eigen::VectorXd> waypoints;
        Eigen::VectorXd waypoint(3);
        waypoint << 0.0, 0.0, 0.0;
        waypoints.push_back(waypoint);
        waypoint << 0.0, 0.2, 1.0;
        waypoints.push_back(waypoint);
        waypoint << 0.0, 3.0, 0.5;
        waypoints.push_back(waypoint);
        waypoint << 1.1, 2.0, 0.0;
        waypoints.push_back(waypoint);
        waypoint << 1.0, 0.0, 0.0;
        waypoints.push_back(waypoint);
        waypoint << 0.0, 1.0, 0.0;
        waypoints.push_back(waypoint);
        waypoint << 0.0, 0.0, 1.0;
        waypoints.push_back(waypoint);

        return waypoints;
}

void printSample(const Trajectory::PhaseSample &sample)
{
        std::cout << std::fixed << std::setprecision(4)
                  << sample.time << "  "
                  << std::setw(8) << sample.pathPos << "  "
                  << std::setw(8) << sample.pathVel << "  "
                  << std::setw(8) << sample.pathAcc << std::endl;
}

} // namespace

int main()
{
        const std::list<Eigen::VectorXd> waypoints = makeDemoWaypoints();

        Eigen::VectorXd maxAcceleration(3);
        maxAcceleration << 1.0, 1.0, 1.0;
        Eigen::VectorXd maxVelocity(3);
        maxVelocity << 1.0, 1.0, 1.0;

        Trajectory trajectory(Path(waypoints, 0.1), maxVelocity, maxAcceleration);
        if(!trajectory.isValid())
        {
                std::cerr << "Time-optimal trajectory generation failed." << std::endl;
                return 1;
        }

        const double sampleTimeStep = 0.005;
        std::vector<Trajectory::PhaseSample> samples = trajectory.samplePhaseTrajectory(sampleTimeStep);
        if(samples.size() < 3)
        {
                std::cerr << "Sampling produced an insufficient number of points." << std::endl;
                return 1;
        }

        std::cout << "Sample count: " << samples.size() << ", last time: " << samples.back().time << std::endl;
        std::cout << "Last original samples (time, s, s_dot, s_ddot):" << std::endl;
        for(std::size_t offset = 5; offset > 0; --offset)
        {
                const std::size_t index = samples.size() >= offset ? samples.size() - offset : 0;
                if(index < samples.size())
                {
                        printSample(samples[index]);
                }
        }

        std::vector<double> reconstructedVel(samples.size(), samples.front().pathVel);
        std::vector<double> reconstructedPos(samples.size(), samples.front().pathPos);
        for(std::size_t i = 1; i < samples.size(); ++i)
        {
                const double dt = samples[i].time - samples[i - 1].time;
                reconstructedVel[i] = reconstructedVel[i - 1] + 0.5 * (samples[i - 1].pathAcc + samples[i].pathAcc) * dt;
                reconstructedPos[i] = reconstructedPos[i - 1] + 0.5 * (samples[i - 1].pathVel + reconstructedVel[i]) * dt;
        }
        std::cout << "Reconstructed end state from original acceleration: " << reconstructedPos.back() << ", " << reconstructedVel.back() << std::endl;

        const double jerkLimit = 5.0; // rad/s^3 equivalent in path coordinates
        const JerkFilterResult result = smoothTrajectoryWithJerkLimit(samples, jerkLimit);

        const Trajectory::PhaseSample &originalEnd = samples.back();
        std::cout << "Original end state (s, s_dot): " << originalEnd.pathPos << ", " << originalEnd.pathVel << std::endl;

        std::cout << "Original max jerk : " << result.maxJerkOriginal << std::endl;
        std::cout << "Filtered max jerk : " << result.maxJerkFiltered << std::endl;
        std::cout << "Endpoint velocity error : " << result.endpointVelocityError << std::endl;
        std::cout << "Endpoint position error : " << result.endpointPositionError << std::endl;

        std::cout << "\nFirst five samples (time, s, s_dot, s_ddot) after smoothing:" << std::endl;
        for(std::size_t i = 0; i < std::min<std::size_t>(5, result.filteredSamples.size()); ++i)
        {
                printSample(result.filteredSamples[i]);
        }

        std::cout << "\nLast five samples after smoothing:" << std::endl;
        const std::size_t n = result.filteredSamples.size();
        for(std::size_t offset = 5; offset > 0; --offset)
        {
                const std::size_t index = n >= offset ? n - offset : 0;
                if(index < n)
                {
                        printSample(result.filteredSamples[index]);
                }
        }

        const bool jerkSatisfied = result.maxJerkFiltered <= jerkLimit + 1e-6;
        const bool endpointsAligned = std::abs(result.endpointVelocityError) < 1e-6 && std::abs(result.endpointPositionError) < 1e-6;
        std::cout << "\nJerk constraint satisfied: " << (jerkSatisfied ? "yes" : "no") << std::endl;
        std::cout << "Endpoints preserved: " << (endpointsAligned ? "yes" : "no") << std::endl;

        return jerkSatisfied && endpointsAligned ? 0 : 2;
}

