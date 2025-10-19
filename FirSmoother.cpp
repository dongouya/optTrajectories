#include "FirSmoother.h"

#include <algorithm>
#include <cmath>
#include <numeric>

using std::size_t;

namespace {

/// @brief 对称镜像延拓索引，保证 FIR 卷积访问合法。
static size_t reflectIndex(long index, size_t size) {
        if(size == 0) {
                return 0;
        }
        const long maxIndex = static_cast<long>(size) - 1;
        while(index < 0 || index > maxIndex) {
                if(index < 0) {
                        index = -index - 1;
                }
                else if(index > maxIndex) {
                        index = 2 * maxIndex - index + 1;
                }
        }
        return static_cast<size_t>(index);
}

/// @brief 根据边界误差估计全局时间拉伸因子。
static double computeTimeScale(double originalDuration, double allowedScale, bool allowScaling,
                double startVelocity, double endVelocity, double targetPosition, double achievedPosition)
{
        if(!allowScaling || originalDuration <= 0.0) {
                return 1.0;
        }
        const double velTolerance = 1e-12;
        if(std::abs(startVelocity) > velTolerance || std::abs(endVelocity) > velTolerance) {
                return 1.0;
        }
        if(achievedPosition <= 0.0) {
                return 1.0;
        }
        double scale = targetPosition / achievedPosition;
        if(scale < 1.0) {
                scale = 1.0;
        }
        if(scale > allowedScale) {
                scale = allowedScale;
        }
        return scale;
}

} // namespace

SmoothedTimeLaw FirSmoother::smooth(const Trajectory &trajectory, const FirSmoothingConfig &config) const {
        // 先从轨迹对象中按统一采样周期取样，再复用另一个重载进行平滑。
        std::vector<Trajectory::Sample> samples = trajectory.sampleTimeLaw(config.samplePeriod);
        double dt = 0.0;
        if(samples.size() >= 2) {
                dt = samples[1].time - samples[0].time;
        }
        return smooth(samples, dt, config);
}

SmoothedTimeLaw FirSmoother::smooth(const std::vector<Trajectory::Sample> &samples, double samplePeriod,
                const FirSmoothingConfig &config) const
{
        SmoothedTimeLaw result;
        if(samples.empty()) {
                // 无采样数据时直接返回空结构，避免后续计算。
                return result;
        }

        const size_t n = samples.size();
        if(samplePeriod <= 0.0 && n >= 2) {
                samplePeriod = samples[1].time - samples[0].time;
        }
        if(samplePeriod <= 0.0) {
                samplePeriod = config.samplePeriod;
        }
        double currentSamplePeriod = samplePeriod;

        result.time.resize(n);
        result.position.resize(n);
        result.velocity.resize(n);
        result.acceleration.resize(n);

        for(size_t i = 0; i < n; ++i) {
                result.time[i] = samples[i].time;
        }

        const double totalTime = result.time.back();
        result.diagnostics.originalDuration = totalTime;

        // 使用三次基底重建首末一致的基准轨迹，确保卷积只作用于残差部分。
        std::vector<double> basePosition(n);
        std::vector<double> baseVelocity(n);
        std::vector<double> baseAcceleration(n);
        const double s0 = samples.front().position;
        const double sT = samples.back().position;
        const double v0 = samples.front().velocity;
        const double vT = samples.back().velocity;

        if(totalTime > 0.0) {
                const double a2 = (3.0 * (sT - s0) - (2.0 * v0 + vT) * totalTime) / (totalTime * totalTime);
                const double a3 = (2.0 * (s0 - sT) + (v0 + vT) * totalTime) / (totalTime * totalTime * totalTime);
                for(size_t i = 0; i < n; ++i) {
                        const double t = result.time[i];
                        const double t2 = t * t;
                        basePosition[i] = s0 + v0 * t + a2 * t2 + a3 * t2 * t;
                        baseVelocity[i] = v0 + 2.0 * a2 * t + 3.0 * a3 * t2;
                        baseAcceleration[i] = 2.0 * a2 + 6.0 * a3 * t;
                }
        }
        else {
                for(size_t i = 0; i < n; ++i) {
                        basePosition[i] = s0;
                        baseVelocity[i] = v0;
                        baseAcceleration[i] = 0.0;
                }
        }

        std::vector<double> residualAcceleration(n);
        std::vector<double> rawAcceleration(n);
        for(size_t i = 0; i < n; ++i) {
                rawAcceleration[i] = samples[i].acceleration;
                residualAcceleration[i] = rawAcceleration[i] - baseAcceleration[i];
        }

        std::size_t kernelLength1 = config.primaryKernelLength;
        if(std::isfinite(config.jerkLimit) && samplePeriod > 0.0) {
                double maxDelta = 0.0;
                for(size_t i = 1; i < n; ++i) {
                        maxDelta = std::max(maxDelta, std::abs(rawAcceleration[i] - rawAcceleration[i - 1]));
                }
                if(maxDelta > 0.0 && config.jerkLimit > 0.0) {
                        const double requiredLength = std::ceil(maxDelta / (config.jerkLimit * samplePeriod));
                        if(requiredLength > static_cast<double>(kernelLength1)) {
                                kernelLength1 = static_cast<size_t>(requiredLength);
                        }
                }
        }
        ensureOddLength(kernelLength1);
        std::vector<double> kernel1 = makeMovingAverageKernel(kernelLength1);
        // 第一阶段移动平均滤波，降低 jerk 台阶。
        std::vector<double> filteredResidual = applySymmetricFir(residualAcceleration, kernel1);

        std::size_t kernelLength2 = 0;
        if(config.enableSecondStage) {
                kernelLength2 = config.secondaryKernelLength == 0 ? kernelLength1 : config.secondaryKernelLength;
                ensureOddLength(kernelLength2);
                std::vector<double> kernel2 = makeMovingAverageKernel(kernelLength2);
                // 若开启二次滤波则继续平滑残差。
                filteredResidual = applySymmetricFir(filteredResidual, kernel2);
        }

        std::vector<double> filteredAcceleration(n);
        for(size_t i = 0; i < n; ++i) {
                filteredAcceleration[i] = baseAcceleration[i] + filteredResidual[i];
        }

        integrateKinematics(filteredAcceleration, currentSamplePeriod, v0, s0, result.velocity, result.position);

        const double achievedPosition = result.position.back();
        double timeScale = computeTimeScale(totalTime, config.maxTimeScale, config.allowTimeScale,
                        v0, vT, sT, achievedPosition);
        if(timeScale > 1.0 + 1e-12) {
                // 若末端位置偏差较大且允许拉伸，则等比扩展时间，确保覆盖目标位移。
                currentSamplePeriod *= timeScale;
                for(size_t i = 0; i < n; ++i) {
                        result.time[i] *= timeScale;
                        filteredAcceleration[i] /= (timeScale * timeScale);
                }
                integrateKinematics(filteredAcceleration, currentSamplePeriod, v0, s0, result.velocity, result.position);
                result.diagnostics.timeScaled = true;
                result.diagnostics.timeScale = timeScale;
        }
        else {
                result.diagnostics.timeScaled = false;
                result.diagnostics.timeScale = 1.0;
        }

        const double activeTotalTime = result.time.back();
        const double finalVelocityError = vT - result.velocity.back();
        if(std::abs(finalVelocityError) > config.tolerance && activeTotalTime > 0.0) {
                const double correction = finalVelocityError / activeTotalTime;
                for(size_t i = 0; i < n; ++i) {
                        filteredAcceleration[i] += correction;
                }
                integrateKinematics(filteredAcceleration, currentSamplePeriod, v0, s0, result.velocity, result.position);
        }
        result.diagnostics.velocityBoundaryError = vT - result.velocity.back();

        const double positionError = sT - result.position.back();
        if(std::abs(positionError) > config.tolerance && activeTotalTime > 0.0) {
                // 基于三次多项式的校正项以保持首末位置完全一致。
                for(size_t i = 0; i < n; ++i) {
                        const double t = result.time[i];
                        filteredAcceleration[i] += 6.0 * positionError / (activeTotalTime * activeTotalTime)
                                - 12.0 * positionError * t / (activeTotalTime * activeTotalTime * activeTotalTime);
                }
                integrateKinematics(filteredAcceleration, currentSamplePeriod, v0, s0, result.velocity, result.position);
        }
        result.diagnostics.positionBoundaryError = sT - result.position.back();

        result.acceleration = filteredAcceleration;
        result.diagnostics.smoothedDuration = result.time.back();
        result.diagnostics.kernelLengthStage1 = kernelLength1;
        result.diagnostics.kernelLengthStage2 = config.enableSecondStage ? kernelLength2 : 0;
        evaluateJerkMetrics(filteredAcceleration, currentSamplePeriod, result.diagnostics.jerkPeak, result.diagnostics.jerkRms);

        return result;
}

std::vector<double> FirSmoother::makeMovingAverageKernel(std::size_t length) {
        ensureOddLength(length);
        if(length == 0) {
                length = 1;
        }
        const double gain = 1.0 / static_cast<double>(length);
        std::vector<double> kernel(length, gain);
        return kernel;
}

void FirSmoother::ensureOddLength(std::size_t &length) {
        if(length == 0) {
                length = 1;
        }
        if(length % 2 == 0) {
                ++length;
        }
}

std::vector<double> FirSmoother::applySymmetricFir(const std::vector<double> &signal, const std::vector<double> &kernel) {
        const size_t n = signal.size();
        const size_t k = kernel.size();
        std::vector<double> output(n, 0.0);
        if(n == 0 || k == 0) {
                        return output;
        }
        const size_t half = k / 2;
        for(size_t i = 0; i < n; ++i) {
                double acc = 0.0;
                for(size_t j = 0; j < k; ++j) {
                        const long offset = static_cast<long>(j) - static_cast<long>(half);
                        const size_t index = reflectIndex(static_cast<long>(i) + offset, n);
                        acc += kernel[j] * signal[index];
                }
                output[i] = acc;
        }
        return output;
}

void FirSmoother::integrateKinematics(const std::vector<double> &acceleration, double samplePeriod, double initialVelocity,
                double initialPosition, std::vector<double> &velocityOut, std::vector<double> &positionOut)
{
        const size_t n = acceleration.size();
        velocityOut.resize(n);
        positionOut.resize(n);
        if(n == 0) {
                return;
        }
        velocityOut[0] = initialVelocity;
        positionOut[0] = initialPosition;
        for(size_t i = 1; i < n; ++i) {
                const double dv = 0.5 * (acceleration[i - 1] + acceleration[i]) * samplePeriod;
                velocityOut[i] = velocityOut[i - 1] + dv;
                const double avgVel = 0.5 * (velocityOut[i - 1] + velocityOut[i]);
                positionOut[i] = positionOut[i - 1] + avgVel * samplePeriod;
        }
}

void FirSmoother::evaluateJerkMetrics(const std::vector<double> &acceleration, double samplePeriod,
                double &jerkPeak, double &jerkRms)
{
        jerkPeak = 0.0;
        jerkRms = 0.0;
        const size_t n = acceleration.size();
        if(n < 2 || samplePeriod <= 0.0) {
                return;
        }
        std::vector<double> jerk(n - 1);
        for(size_t i = 1; i < n; ++i) {
                jerk[i - 1] = (acceleration[i] - acceleration[i - 1]) / samplePeriod;
        }
        double sumSquares = 0.0;
        for(double j : jerk) {
                jerkPeak = std::max(jerkPeak, std::abs(j));
                sumSquares += j * j;
        }
        if(!jerk.empty()) {
                jerkRms = std::sqrt(sumSquares / static_cast<double>(jerk.size()));
        }
}

