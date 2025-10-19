#pragma once

#include <cstddef>
#include <vector>
#include <limits>

#include "Trajectory.h"

/**
 * @brief FIR 平滑器配置项。
 * 这些参数控制采样周期、卷积核长度以及允许的时间伸缩范围，用以满足 jerk 限制并保持边界条件。
 */
struct FirSmoothingConfig {
        double samplePeriod = 0.002;                 ///< 期望的离散采样周期（秒）。
        std::size_t primaryKernelLength = 41;        ///< 主平滑窗口长度，必须为奇数。
        std::size_t secondaryKernelLength = 0;       ///< 第二阶段平滑窗口长度，0 表示与第一阶段一致。
        bool enableSecondStage = false;              ///< 若为 true，则执行双阶段移动平均平滑。
        double jerkLimit = std::numeric_limits<double>::infinity(); ///< 路径域 jerk 约束，用于自适应估计核长。
        double tolerance = 1e-6;                     ///< 数值修正阈值，保持首末速度/位置一致性。
        bool allowTimeScale = true;                  ///< 允许整体拉伸时间轴以进一步降低 jerk。
        double maxTimeScale = 1.3;                   ///< 时间拉伸因子上限，防止持续扩张。
};

/**
 * @brief 平滑诊断信息，帮助分析滤波效果与约束满足情况。
 */
struct FirSmoothingDiagnostics {
        double originalDuration = 0.0;       ///< 原始时间律总时长。
        double smoothedDuration = 0.0;       ///< 平滑后的时间律总时长。
        double timeScale = 1.0;              ///< 实际使用的时间拉伸因子。
        bool timeScaled = false;             ///< 是否应用了时间拉伸。
        double jerkPeak = 0.0;               ///< 平滑后 jerk 峰值估计。
        double jerkRms = 0.0;                ///< 平滑后 jerk 均方根估计。
        std::size_t kernelLengthStage1 = 0;  ///< 第一阶段卷积核长度。
        std::size_t kernelLengthStage2 = 0;  ///< 第二阶段卷积核长度。
        double velocityBoundaryError = 0.0;  ///< 终端速度误差，用于验证边界保持。
        double positionBoundaryError = 0.0;  ///< 终端位置误差，用于验证边界保持。
};

/**
 * @brief 离散化后的平滑时间律结果。
 */
struct SmoothedTimeLaw {
        std::vector<double> time;            ///< 等间隔采样时间戳。
        std::vector<double> position;        ///< 路径参数位置序列 s(t)。
        std::vector<double> velocity;        ///< 路径参数速度序列 \dot{s}(t)。
        std::vector<double> acceleration;    ///< 路径参数加速度序列 \ddot{s}(t)。
        FirSmoothingDiagnostics diagnostics; ///< 关联的诊断信息。
};

class FirSmoother {
public:
        /// @brief 直接从 Trajectory 对象采样并执行平滑。
        SmoothedTimeLaw smooth(const Trajectory &trajectory, const FirSmoothingConfig &config) const;
        /// @brief 对外部提供的离散采样执行平滑处理。
        SmoothedTimeLaw smooth(const std::vector<Trajectory::Sample> &samples, double samplePeriod, const FirSmoothingConfig &config) const;

private:
        /// @brief 构造单位增益的移动平均核。
        static std::vector<double> makeMovingAverageKernel(std::size_t length);
        /// @brief 将窗口长度调整为合法的奇数。
        static void ensureOddLength(std::size_t &length);
        /// @brief 采用镜像延拓的对称 FIR 卷积。
        static std::vector<double> applySymmetricFir(const std::vector<double> &signal, const std::vector<double> &kernel);
        /// @brief 对加速度序列积分，获得速度与位置。
        static void integrateKinematics(const std::vector<double> &acceleration, double samplePeriod, double initialVelocity,
                        double initialPosition, std::vector<double> &velocityOut, std::vector<double> &positionOut);
        /// @brief 估计 jerk 的峰值与均方根指标。
        static void evaluateJerkMetrics(const std::vector<double> &acceleration, double samplePeriod,
                        double &jerkPeak, double &jerkRms);
};

