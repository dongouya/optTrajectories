#include "JerkFilter.h"

#include <algorithm>
#include <cmath>
#include <numeric>

// The implementation in this file follows the derivation sketched in
// TOTG_jerk_smoothing_guide.md.  In broad strokes we:
//   1. Mirror-pad and convolve the acceleration profile with a binomial kernel to
//      attenuate high-frequency components while preserving the trajectory endpoints.
//   2. Clamp the finite differences to the requested jerk bound in both forward and
//      backward passes (a discrete total-variation-style projection).
//   3. Correct the first two moments of the acceleration error so that integration
//      recovers the original position/velocity boundary conditions.
//   4. Re-integrate the corrected acceleration profile to yield smoothed samples that
//      can replace the original TOTG discretization without further adjustments.
// Each helper below performs one of these conceptual steps.

namespace
{

/// Build a (normalized) binomial filter kernel of the requested length.  The kernel
/// acts as a simple low-pass filter when convolved with the acceleration profile.
std::vector<double> buildBinomialKernel(std::size_t length)
{
        std::vector<double> kernel(length, 0.0);
        if(length == 0)
        {
                return kernel;
        }

        const std::size_t order = length - 1;
        double coefficient = 1.0;
        double sum = 0.0;
        for(std::size_t i = 0; i < length; ++i)
        {
                kernel[i] = coefficient;
                sum += coefficient;
                if(i < order)
                {
                        coefficient *= static_cast<double>(order - i) / static_cast<double>(i + 1);
                }
        }

        if(sum > 0.0)
        {
                for(double &value : kernel)
                {
                        value /= sum;
                }
        }

        return kernel;
}

/// Mirror-pad the input signal so the convolution in convolveSymmetric behaves as if
/// the acceleration continued smoothly past the boundaries.  This dramatically reduces
/// ringing that would otherwise spoil the end-point constraints.
std::vector<double> mirrorPad(const std::vector<double> &data, std::size_t pad)
{
        if(pad == 0 || data.empty())
        {
                return data;
        }

        const std::size_t n = data.size();
        std::size_t effectivePad = std::min(pad, n - 1);
        std::vector<double> padded(2 * effectivePad + n);

        // mirror at the beginning
        for(std::size_t i = 0; i < effectivePad; ++i)
        {
                padded[effectivePad - 1 - i] = data[std::min(i, n - 1)];
        }

        // copy original data
        std::copy(data.begin(), data.end(), padded.begin() + effectivePad);

        // mirror at the end
        for(std::size_t i = 0; i < effectivePad; ++i)
        {
                padded[effectivePad + n + i] = data[n - 1 - std::min(i, n - 1)];
        }

        return padded;
}

/// Symmetric convolution that applies the kernel produced by buildBinomialKernel to the
/// acceleration samples.  The combination of mirrorPad + convolution corresponds to the
/// filtering stage from the guide.
std::vector<double> convolveSymmetric(const std::vector<double> &data, const std::vector<double> &kernel)
{
        if(kernel.empty() || data.empty())
        {
                return data;
        }

        const std::size_t pad = kernel.size() / 2;
        std::vector<double> padded = mirrorPad(data, pad);
        std::vector<double> result(data.size(), 0.0);

        for(std::size_t index = 0; index < data.size(); ++index)
        {
                double value = 0.0;
                for(std::size_t k = 0; k < kernel.size(); ++k)
                {
                        value += kernel[k] * padded[index + k];
                }
                result[index] = value;
        }

        return result;
}

/// Enforce |jerk| <= jerkLimit by clamping finite differences in forward and backward
/// passes.  This discrete projection keeps the filtered signal within the requested
/// jerk budget while remaining inexpensive to evaluate.
void limitJerk(std::vector<double> &acceleration, const std::vector<double> &timeSteps, double jerkLimit)
{
        if(acceleration.empty())
        {
                return;
        }

        for(std::size_t i = 1; i < acceleration.size(); ++i)
        {
                const double limit = jerkLimit * timeSteps[i - 1];
                const double delta = acceleration[i] - acceleration[i - 1];
                if(delta > limit)
                {
                        acceleration[i] = acceleration[i - 1] + limit;
                }
                else if(delta < -limit)
                {
                        acceleration[i] = acceleration[i - 1] - limit;
                }
        }

        for(std::size_t i = acceleration.size(); i-- > 1;)
        {
                const double limit = jerkLimit * timeSteps[i - 1];
                const double delta = acceleration[i - 1] - acceleration[i];
                if(delta > limit)
                {
                        acceleration[i - 1] = acceleration[i] + limit;
                }
                else if(delta < -limit)
                {
                        acceleration[i - 1] = acceleration[i] - limit;
                }
        }
}

/// Utility function used for the diagnostic report in JerkFilterResult.
double computeMaxJerk(const std::vector<double> &acceleration, const std::vector<double> &timeSteps)
{
        double maxJerk = 0.0;
        if(acceleration.size() < 2 || timeSteps.empty())
        {
                return maxJerk;
        }

        for(std::size_t i = 1; i < acceleration.size(); ++i)
        {
                const double dt = timeSteps[i - 1];
                if(dt <= 0.0)
                {
                        continue;
                }
                const double jerk = std::abs(acceleration[i] - acceleration[i - 1]) / dt;
                maxJerk = std::max(maxJerk, jerk);
        }
        return maxJerk;
}

} // namespace

// The public entry point simply orchestrates the steps above and performs the moment
// corrections/integration required by the guide.
JerkFilterResult smoothTrajectoryWithJerkLimit(const std::vector<Trajectory::PhaseSample> &samples, double jerkLimit)
{
        JerkFilterResult result;
        result.filteredSamples = samples;

        if(samples.size() < 3 || jerkLimit <= 0.0)
        {
                if(samples.size() >= 2)
                {
                        std::vector<double> originalAcc(samples.size());
                        for(std::size_t i = 0; i < samples.size(); ++i)
                        {
                                originalAcc[i] = samples[i].pathAcc;
                        }
                        std::vector<double> timeSteps(samples.size() - 1, 1.0);
                        for(std::size_t i = 0; i + 1 < samples.size(); ++i)
                        {
                                double dt = samples[i + 1].time - samples[i].time;
                                if(dt > 0.0)
                                {
                                        timeSteps[i] = dt;
                                }
                        }
                        result.maxJerkOriginal = computeMaxJerk(originalAcc, timeSteps);
                        result.maxJerkFiltered = result.maxJerkOriginal;
                }
                return result;
        }

        const double startTime = samples.front().time;
        const double endTime = samples.back().time;
        double duration = endTime - startTime;
        if(duration <= 0.0)
        {
                duration = static_cast<double>(samples.size() - 1);
        }
        double averageTimeStep = duration / static_cast<double>(samples.size() - 1);
        if(averageTimeStep <= 0.0)
        {
                averageTimeStep = 1.0;
        }

        std::vector<double> timeSteps(samples.size() - 1, averageTimeStep);
        for(std::size_t i = 0; i < timeSteps.size(); ++i)
        {
                double dt = samples[i + 1].time - samples[i].time;
                if(dt > 0.0)
                {
                        timeSteps[i] = dt;
                }
        }

        std::vector<double> acceleration(samples.size());
        std::transform(samples.begin(), samples.end(), acceleration.begin(), [](const Trajectory::PhaseSample &sample) {
                return sample.pathAcc;
        });

        double maxAbsAcceleration = 0.0;
        for(const double value : acceleration)
        {
                maxAbsAcceleration = std::max(maxAbsAcceleration, std::abs(value));
        }
        if(maxAbsAcceleration == 0.0)
        {
                maxAbsAcceleration = jerkLimit * averageTimeStep;
        }

        const double smoothingTime = maxAbsAcceleration / jerkLimit;
        std::size_t kernelSize = static_cast<std::size_t>(std::ceil(smoothingTime / averageTimeStep));
        if(kernelSize % 2 == 0)
        {
                kernelSize += 1;
        }
        if(kernelSize < 3)
        {
                kernelSize = 3;
        }
        if(kernelSize > samples.size())
        {
                kernelSize = samples.size() - (samples.size() % 2 == 0 ? 1 : 0);
        }
        if(kernelSize < 3)
        {
                kernelSize = 3;
        }

        const std::vector<double> kernel = buildBinomialKernel(kernelSize);
        std::vector<double> smoothedAcceleration = convolveSymmetric(acceleration, kernel);

        limitJerk(smoothedAcceleration, timeSteps, jerkLimit);

        auto computeBasis = [&](double xi) {
                const double xi2 = xi * xi;
                const double xi3 = xi2 * xi;
                const double xi4 = xi3 * xi;
                const double xi5 = xi4 * xi;
                const double b0 = 10.0 * xi3 - 15.0 * xi4 + 6.0 * xi5;
                const double b1 = xi * (1.0 - xi) * (1.0 - 2.0 * xi);
                return std::make_pair(b0, b1);
        };

        std::vector<double> basis0(samples.size(), 0.0);
        std::vector<double> basis1(samples.size(), 0.0);
        for(std::size_t i = 0; i < samples.size(); ++i)
        {
                const double timeFromStart = samples[i].time - startTime;
                const double xi = (duration <= 0.0) ? 0.0 : timeFromStart / duration;
                const auto basis = computeBasis(xi);
                basis0[i] = basis.first;
                basis1[i] = basis.second;
        }

        double phi00 = 0.0, phi01 = 0.0, phi10 = 0.0, phi11 = 0.0;
        for(std::size_t i = 1; i < samples.size(); ++i)
        {
                const double prevTime = samples[i - 1].time - startTime;
                const double currTime = samples[i].time - startTime;
                const double dt = timeSteps[i - 1];
                phi00 += 0.5 * (basis0[i - 1] + basis0[i]) * dt;
                phi01 += 0.5 * ((duration - prevTime) * basis0[i - 1] + (duration - currTime) * basis0[i]) * dt;
                phi10 += 0.5 * (basis1[i - 1] + basis1[i]) * dt;
                phi11 += 0.5 * ((duration - prevTime) * basis1[i - 1] + (duration - currTime) * basis1[i]) * dt;
        }

        std::vector<double> correctedAcceleration = smoothedAcceleration;
        const double determinant = phi00 * phi11 - phi10 * phi01;
        bool converged = false;

        for(int iteration = 0; iteration < 6; ++iteration)
        {
                limitJerk(correctedAcceleration, timeSteps, jerkLimit);

                std::vector<double> difference(correctedAcceleration.size());
                for(std::size_t i = 0; i < correctedAcceleration.size(); ++i)
                {
                        difference[i] = correctedAcceleration[i] - acceleration[i];
                }

                double deltaVelocity = 0.0;
                double deltaPosition = 0.0;
                for(std::size_t i = 1; i < difference.size(); ++i)
                {
                        const double prevDifference = difference[i - 1];
                        const double currDifference = difference[i];
                        const double dt = timeSteps[i - 1];
                        deltaVelocity += 0.5 * (prevDifference + currDifference) * dt;

                        const double prevTime = samples[i - 1].time - startTime;
                        const double currTime = samples[i].time - startTime;
                        const double prevIntegrand = (duration - prevTime) * prevDifference;
                        const double currIntegrand = (duration - currTime) * currDifference;
                        deltaPosition += 0.5 * (prevIntegrand + currIntegrand) * dt;
                }

                if(std::abs(deltaVelocity) < 1e-9 && std::abs(deltaPosition) < 1e-9)
                {
                        converged = true;
                        break;
                }

                double alpha = 0.0;
                double beta = 0.0;
                if(std::abs(determinant) > 1e-9)
                {
                        alpha = (-deltaVelocity * phi11 + deltaPosition * phi10) / determinant;
                        beta = (-deltaPosition * phi00 + deltaVelocity * phi01) / determinant;
                }

                for(std::size_t i = 0; i < correctedAcceleration.size(); ++i)
                {
                        correctedAcceleration[i] += alpha * basis0[i] + beta * basis1[i];
                }
        }
        if(!converged)
        {
                limitJerk(correctedAcceleration, timeSteps, jerkLimit);
        }

        // Integrate acceleration to recover velocity and position
        auto integrateAcceleration = [&]() {
                result.filteredSamples.front().pathPos = samples.front().pathPos;
                result.filteredSamples.front().pathVel = samples.front().pathVel;
                result.filteredSamples.front().pathAcc = correctedAcceleration.front();

                for(std::size_t i = 1; i < correctedAcceleration.size(); ++i)
                {
                        const double dt = timeSteps[i - 1];
                        const double velocity = result.filteredSamples[i - 1].pathVel + 0.5 * (correctedAcceleration[i - 1] + correctedAcceleration[i]) * dt;
                        const double position = result.filteredSamples[i - 1].pathPos + 0.5 * (result.filteredSamples[i - 1].pathVel + velocity) * dt;

                        result.filteredSamples[i].time = samples[i].time;
                        result.filteredSamples[i].pathVel = velocity;
                        result.filteredSamples[i].pathPos = position;
                        result.filteredSamples[i].pathAcc = correctedAcceleration[i];
                }
        };

        integrateAcceleration();

        for(int iteration = 0; iteration < 3; ++iteration)
        {
                const double velocityError = result.filteredSamples.back().pathVel - samples.back().pathVel;
                const double positionError = result.filteredSamples.back().pathPos - samples.back().pathPos;
                if(std::abs(velocityError) < 1e-6 && std::abs(positionError) < 1e-6)
                {
                        break;
                }

                double alpha = 0.0;
                double beta = 0.0;
                if(std::abs(determinant) > 1e-9)
                {
                        alpha = (-velocityError * phi11 + positionError * phi10) / determinant;
                        beta = (-positionError * phi00 + velocityError * phi01) / determinant;
                }

                for(std::size_t i = 0; i < correctedAcceleration.size(); ++i)
                {
                        correctedAcceleration[i] += alpha * basis0[i] + beta * basis1[i];
                }
                limitJerk(correctedAcceleration, timeSteps, jerkLimit);
                integrateAcceleration();
        }

        result.maxJerkOriginal = computeMaxJerk(acceleration, timeSteps);
        result.maxJerkFiltered = computeMaxJerk(correctedAcceleration, timeSteps);
        result.endpointVelocityError = result.filteredSamples.back().pathVel - samples.back().pathVel;
        result.endpointPositionError = result.filteredSamples.back().pathPos - samples.back().pathPos;

        return result;
}

