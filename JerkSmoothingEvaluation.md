# Jerk-Limited Post-Processing for Kunz–Stilman TOTG

## 概述

为了验证 `TOTG_jerk_smoothing_guide.md` 中提出的滤波方案是否可行，本次实现遵循其主要步骤：

1. **统一采样并补齐端点**：为原始 TOTG 轨迹提供均匀的时间采样接口（`Trajectory::samplePhaseTrajectory`）。
2. **对称核平滑**：使用归一化二项式核对路径加速度作零相位卷积。
3. **jerk 斜率限幅**：在离散采样上执行双向扫描，确保 \\(|\Delta a| \le J_{\max} \Delta t\\)。
4. **两矩投影纠偏**：利用与指南一致的基函数（五次多项式）迭代矫正速度增量和位置增量。
5. **端点精调**：为应对 TOTG 离散加速度与其几何端点轻微不一致的问题，再次使用同一基函数对末端速度、末端位置残差做小幅修正，同时保持 jerk 限制。
6. **重新积分与验证**：积分得到平滑后的 \\(s,\dot s,\ddot s\\)，并输出 jerk、端点误差等指标。

## 方案评估

- **对称核与限幅**：二项式核配合镜像延拓能够有效抑制跳变，之后的 jerk 限幅把原始 TOTG 的 570 rad/s^3 量级 jerk 压缩到给定上界 5 rad/s^3。【F:jerk_demo_output.txt†L10-L12】
- **端点保持**：两阶段的“两个矩”纠偏后，末速度、末位置与原始 TOTG 完全一致，输出误差均为 0。【F:jerk_demo_output.txt†L12-L13】【F:jerk_demo_output.txt†L20-L27】
- **离散积分偏差**：实测 TOTG 原始的离散加速度在最后一个不完整采样步上存在微小误差（原始加速度积分与轨迹文件差了约 0.058 rad）。新增的端点精调循环能可靠地抵消这一残差。【F:jerk_demo_output.txt†L9-L11】
- **jerk 约束验证**：`JerkFilterResult` 在求解后提供最大 jerk、端点误差等指标，可直接用于回归测试。

综上，`TOTG_jerk_smoothing_guide.md` 所描述的滤波方案在 Kunz–Stilman TOTG 上可行，且经过上述实现能够生成 jerk 限制明确、端点与原规划一致的平滑轨迹。

## 关键实现文件

- `Trajectory::samplePhaseTrajectory`：输出统一采样的相位轨迹，为滤波模块提供输入。【F:Trajectory.cpp†L466-L534】
- `JerkFilter.{h,cpp}`：指南的核心步骤（平滑、jerk 限幅、两矩纠偏、端点修正与重新积分）。【F:JerkFilter.cpp†L1-L363】
- `JerkSmoothingDemo.cpp`：演示与验证入口，生成统计数据并对照原始 TOTG 结果。【F:JerkSmoothingDemo.cpp†L1-L115】

## 使用方法

```bash
g++ -std=c++17 JerkSmoothingDemo.cpp Trajectory.cpp Path.cpp JerkFilter.cpp -I. -o jerk_demo
./jerk_demo
```

程序会输出 jerk 约束前后的最大值、端点误差以及首尾若干采样点，可用于检查滤波效果。

