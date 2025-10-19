# FIR-Based Jerk-Constrained Smoothing Guide  
*for `tobiaskunz/trajectories` outputs*  
**File:** `markdown.md`  
**Purpose:** Instruct Codex to implement an offline, zero-phase FIR smoothing pipeline that (i) takes the **velocity/acceleration-limited** time law from `tobiaskunz/trajectories`, (ii) **suppresses jerk** via FIR shaping, (iii) **preserves** start/end **position & velocity**, and (iv) **allows reasonable duration growth** to respect the jerk bound.

---

## 0) Goals & Non-Goals

**Goals**
- 输入：已知几何路径（或已离散的关节轨迹） + Kunz TOTG 输出的**时间参数化** `s(t)`（满足 \( \dot s, \ddot s \) 上限，但不含 jerk 约束）。
- 输出：平滑后的时间律 `s̃(t)`（或离散序列），满足：
  - 起止位置/速度 **完全一致**（不破坏边界条件）；
  - **Jerk 显著抑制/有界**（可设置软/硬阈值）；
  - **总时长允许增加**，且**尽量最小化**时长增长；
  - 原速度/加速度上限仍被满足（如需，进行投影/裁剪）。

**非目标**
- 不改变**几何路径**（仅对时间律整形）。
- 不要求在线因果实现；采用**离线零相位 FIR**（双向）以避免相位延迟。

---

## 1) 名词与符号

- 路径参数：\( s \in [0, s_{\mathrm{end}}] \)，时间律：\( s(t) \)；采样周期 \( \Delta t \)。
- 一阶到三阶时间导数：\( \dot s, \ddot s, \dddot s \)。
- 关节/笛卡尔量：\( \mathbf q(s) \)；链式法则：
  \[
  \dot{\mathbf q}=\mathbf q'(s)\,\dot s,\quad
  \ddot{\mathbf q}=\mathbf q''(s)\,\dot s^2+\mathbf q'(s)\,\ddot s,\quad
  \dddot{\mathbf q}=\mathbf q'''(s)\,\dot s^3+3\mathbf q''(s)\,\dot s\,\ddot s+\mathbf q'(s)\,\dddot s.
  \]
- 速度/加速度/jerk 上界：\( \dot s \le \dot s_{\max},\ \ddot s \le \ddot s_{\max},\ |\dddot q_i|\le J_{i,\max} \)（或在路径域上给定）。
- FIR 核：对称系数 \( h[m],\ m=-(M\!-\!1)/2,\dots,(M\!-\!1)/2 \)，零相位且 \( \sum h[m]=1 \)。

---

## 2) 系统集成概览

**Pipeline（离线）**
1. **从 Kunz 输出取样**：从其分段时间律（恒加速度段）生成等间隔序列 \( \{ s_k,\ \dot s_k,\ \ddot s_k \} \)，\( k=0..K \)。
2. **端点基底**（保持边界）：构造首末三次多项式基底 \( b(t) \)，使 \( b \) 与原 \( s \) 在 \( t=0,T_0 \) 处的 \( (s,\dot s,\ddot s) \) 一致；把信号写为
   \[
   s(t) = b(t) + r(t),\quad \ddot s(t) = \ddot b(t) + \ddot r(t),
   \]
   只对残差 \( \ddot r \) 做滤波（避免边界被卷积破坏）。
3. **FIR 设计**：选对称 FIR（单级/双级 MA 或优化 FIR），给出**核长**以满足目标 jerk 上界。
4. **零相位卷积**：对 \( \ddot r_k \) 做前向-后向（或直接对称 FIR）卷积，得 \( \ddot r_k^{(f)} \)。
5. **积分重建**：\(\ddot s^{(f)}=\ddot b+\ddot r^{(f)}\)\ \(\rightarrow\)\ 数值积分得 \( \dot s^{(f)},s^{(f)} \)。
6. **约束投影**：若 \( \dot s^{(f)}, \ddot s^{(f)} \) 超界，进行幅值投影/缩放；保持**单调性** \( \dot s^{(f)}\ge 0 \)。
7. **时标延展**：允许 \( T \uparrow \)。通过“**速度面积守恒 + 尾段延展**”或“**全局时间缩放**”使 \( s^{(f)}(T_{\mathrm{new}})=s_{\mathrm{end}} \) 且终端 \( \dot s^{(f)}(T_{\mathrm{new}})=\dot s(T_0) \)（常为 0）。
8. **关节 jerk 校核**：由链式法则评估 \( \dddot{\mathbf q} \) 峰值与 RMS；若不达标，增大核长或放宽 \( T \)。
9. **输出**：平滑后的离散时间律或新的分段时间律（支持回写为分段多项式）。

---

## 3) 输入/输出与模块划分

### 3.1 输入接口（建议）
- `PathSamples`: 采样路径几何（或直接提供 \( \mathbf q(s), \mathbf q'(s), \mathbf q''(s) \) 的评估器）。
- `TotgSchedule`: Kunz 输出的分段时间律或等间隔采样 `s_k, sdot_k, sddot_k`。
- `Constraints`: `{s_dot_max, s_ddot_max, JerkMaxPerJoint[]}`。
- `Config`: `{dt, fir_type, fir_len, two_stage, endpoint_mode, allow_time_growth=true, growth_cap_ratio, projection_iters, tolerance}`。

### 3.2 输出接口（建议）
- `SmoothedSchedule`: `t_new[], s_new[], sdot_new[], sddot_new[]`（均等间隔或非均匀时间网格）。
- `Diagnostics`: 峰值 jerk、RMS、总时长变化、是否投影/缩放、是否触发尾段延展等。
- （可选）回写为**分段多项式**，便于控制器插补。

---

## 4) 采样与端点基底

### 4.1 等间隔采样
- 选择 \( \Delta t \) 使每个恒加速度段至少有 10–20 个采样点（典型 \( \Delta t=1\text{–}5\ \mathrm{ms} \)）。
- 得到 `s[k], sdot[k], sddot[k]`。若只有解析段，可直接采样；若只有 `sdot`，用差分近似 \( \ddot s \)。

### 4.2 端点基底（严守边界的关键）
- 用三次多项式 \( b(t)=\alpha_0+\alpha_1 t+\alpha_2 t^2+\alpha_3 t^3 \)，满足：
  \[
  b(0)=s(0),\ b(T)=s(T),\ \dot b(0)=\dot s(0),\ \dot b(T)=\dot s(T).
  \]
- 令 \( r(t)=s(t)-b(t) \)，对 \( \ddot r \) 做 FIR。这样滤波不会改变首末 \( s,\dot s \)。  
- 实做：离散上先拟合 \( b_k \)（或解析求系数），计算 \( \ddot r_k=\ddot s_k-\ddot b_k \)。

---

## 5) FIR 设计与核长定尺

### 5.1 推荐 FIR 族
- **单级等权移动平均（MA）**：长度 \( N \)；对加速度**阶跃**输出为线性爬坡，离散 jerk 上界近似
  \[
  |\dddot s|_{\max}\ \approx\ \frac{|\Delta \ddot s|}{N\,\Delta t}.
  \]
- **双级 MA（MA\*MA）**：等价三角 FIR，jerk 峰值更低，过渡更圆滑；保守起见用**同一公式**定尺，实际更安全。
- **优化对称 FIR**：以“对阶跃的最大斜率 ≤ 目标”与“单位 DC 增益”为约束，最小化通带波动；后续可迭代升级。

### 5.2 核长初值（保守选型）
1. 扫描 Kunz 序列中所有加速度**台阶**，取最大 \( |\Delta \ddot s|_{\max} \)。
2. 目标 jerk 上界 \( J_{\mathrm{path}} \)（由关节 \( J_{i,\max} \) 保守映射或经验设定）。
3. 令
   \[
   N\ \ge\ \left\lceil \frac{|\Delta \ddot s|_{\max}}{J_{\mathrm{path}}\,\Delta t} \right\rceil,\quad N\text{ 取奇数}.
   \]
4. 双级 MA 时，两级可取 \( N_1=N_2=N \)；若需缩短总时窗，可取 \( (N_1,N_2)=(N,\ \lfloor 0.6N \rfloor) \) 做折中。

> **提示**：FIR 需对称、系数和为 1（保面积），避免改变速度/位移的 DC 成分。

---

## 6) 零相位卷积与数值重建

### 6.1 零相位实现
- **方式 A**：对称系数直接非因果卷积（需要端点外延：镜像或样条外延 ≥ 半窗）。
- **方式 B（更稳）**：**Forward–Backward**：先正向用因果 FIR 过滤，再反向（时间倒序）重复一次，得到零相位等效。

### 6.2 积分重建与单调性
1. \( \ddot s_k^{(f)}=\ddot b_k+(\ddot r*h)_k \)  
2. \( \dot s_{k+1}^{(f)}=\max\{0,\ \min(\dot s_{\max},\ \dot s_k^{(f)}+\ddot s_k^{(f)}\Delta t)\} \)  
3. \( s_{k+1}^{(f)}=s_k^{(f)}+\dot s_{k+1}^{(f)}\Delta t \)

> **注意**：第 2 步同时确保**单调**（\(\dot s\ge0\)）与**限速**投影。

---

## 7) 约束投影与时长延展

### 7.1 限速/限加投影（如需）
- 若 \( \dot s^{(f)} \) 或 \( \ddot s^{(f)} \) 超界：  
  - **幅值剪裁**：点态裁剪后再做一次轻度平滑（短核）以去除尖角；  
  - 或 **缩放核长**（增大 \( N \)）并重跑。
- 迭代上限例如 3 次；每次记录越界率与峰值。

### 7.2 终端对齐与时长增长
- 经过平滑后，通常 \( s^{(f)}(T_0) < s_{\mathrm{end}} \)（速度平台降低导致**面积不足**）。  
- 方案：
  1) **尾段延展**：保持末端 \(\dot s\) 收敛到目标（常为 0），追加一段最短 S 段（jerk 合法），直到 \( s = s_{\mathrm{end}} \)。  
  2) **全局时间缩放**：把时间轴按 \( \kappa>1 \) 等比拉伸（\(\dot s,\ddot s\) 等比缩小），选最小 \(\kappa\) 使 \( s(T_{\mathrm{new}})=s_{\mathrm{end}} \) 且各上界满足。  

> **产物**：新的时长 \( T_{\mathrm{new}}\ge T_0 \)。记录 \( \Delta T=T_{\mathrm{new}}-T_0 \)。

---

## 8) 多轴 jerk 校核（关键验证）

- 若可获得 \( \mathbf q'(s),\mathbf q''(s),\mathbf q'''(s) \)，按链式法则计算 \( \dddot{\mathbf q}(t_k) \)。  
- 若仅有离散 \( \mathbf q(t) \)，用三阶向后差分近似：
  \[
  \dddot q_i[k] \approx \frac{ -q_i[k-3]+3q_i[k-2]-3q_i[k-1]+q_i[k] }{\Delta t^3 }.
  \]
- 指标：
  - 峰值：\( \max_i \max_k |\dddot q_i[k]| \le J_{i,\max} \)（硬阈值）；
  - RMS：\( \mathrm{RMS}_J=\sqrt{\sum_i \sum_k \dddot q_i[k]^2\,\Delta t / \sum_i T_{\mathrm{new}} } \)（越小越好）。

---

## 9) 主要边界与异常场景处理

- **短段 < FIR 支撑**：对切换点附近使用**自适应更短核**或直接插入**解析 S 段**替换（保证 jerk 硬上界）。
- **密集拐点**：优先采用**双级 MA**（更圆滑），必要时提高 \( \Delta T \) 或略增核长。
- **端点振铃**：启用**基底法 + 外延**；外延可用镜像/样条。
- **负速度伪影**：投影 \( \dot s\ge0 \)，并在该处**延时**（时间延展）而不是强行裁剪面积。

---

## 10) 代码结构建议（C++，供 Codex 生成）

> 仅给出**骨架与关键函数签名**，便于 Codex 填充实现。

```
/project
  /include
    FirSmoother.hpp
    ScheduleIO.hpp
    EndpointBasis.hpp
    ConstraintProjector.hpp
    JerkEvaluator.hpp
  /src
    FirSmoother.cpp
    ScheduleIO.cpp
    EndpointBasis.cpp
    ConstraintProjector.cpp
    JerkEvaluator.cpp
  /apps
    smooth_totg.cpp   // CLI: read Kunz schedule -> smooth -> write
  /tests
    test_sanity.cpp
    test_endpoints.cpp
    test_constraints.cpp
    test_short_segments.cpp
```

**核心 API**
```cpp
struct TotgSchedule { std::vector<double> t, s, sdot, sddot; double T0, send; };
struct SmoothedSchedule { std::vector<double> t, s, sdot, sddot; double Tnew; };

struct Constraints {
  double sdot_max, sddot_max;
  std::vector<double> Jmax_per_joint; // optional
};

struct FirConfig {
  double dt;
  enum Type { MovingAverage, TwoStageMA, Optimized } type;
  int N1, N2; // for TwoStageMA
  bool use_forward_backward = true;
  bool allow_time_growth = true;
  double max_growth_ratio = 1.5;
  int projection_max_iters = 3;
  double tolerance = 1e-9;
};

SmoothedSchedule smoothWithFIR(
    const TotgSchedule& in,
    const Constraints& cons,
    const FirConfig& cfg,
    std::function<void(const SmoothedSchedule&)> progress_cb = nullptr);

struct JerkReport {
  double peak_abs;
  double rms;
  std::vector<double> per_joint_peak; // if q'(s),q''(s),q'''(s) available
};

JerkReport evaluateJerk(
    const SmoothedSchedule& sch,
    std::function<void(double s, double& q1, double& q2, double& q3)> pathDerivsOpt // optional
);
```

**实现要点（供 Codex）**
- `EndpointBasis`：构造三次多项式 \( b(t) \)；生成 `b, bdot, bddot` 离散序列。
- `FirSmoother`：
  - 端点外延函数（镜像/样条）；
  - 前向滤波 + 反向滤波（零相位）；
  - 积分重建（梯形/更高阶积分），确保数值稳定；
  - 单调性与限速投影；
  - 若 `allow_time_growth`：实现**尾段延展**（以 \( \ddot s,\dddot s \) 限制生成最小 S 段）。
- `ConstraintProjector`：幅值裁剪 + 去尖角二次平滑的“小回路”，限迭代。
- `JerkEvaluator`：若提供 \( \mathbf q'(s),\mathbf q''(s),\mathbf q'''(s) \) 回调，用链式法则；否则用离散三阶差分。
- `ScheduleIO`：读写 CSV/JSON；可视化导出（gnuplot/matplotlib 数据）。

**CLI 示例**
```
./smooth_totg \
  --in kunz_schedule.csv \
  --out smooth.csv \
  --dt 0.002 \
  --fir two_stage --N1 41 --N2 31 \
  --sdot_max 1.2 --sddot_max 3.5 \
  --jerk_joint_max 200,180,180,160,160,120 \
  --allow_time_growth 1 --max_growth_ratio 1.3
```

---

## 11) 参数调优流程（实践指南）

1. **设定 \( \Delta t \)**：根据控制回路/目标带宽选 1–5 ms。  
2. **估计 \( |\Delta \ddot s|_{\max} \)**：扫描原序列台阶。  
3. **目标 jerk**（路径域）\( J_{\mathrm{path}} \)**初值**：由最严格关节 \( J_{i,\max} \) 反推保守值，或直接取经验值（如 10–30×路径加速度单位/秒）。  
4. **核长初选**：用定尺公式得到 \( N \)；双级 MA 设置 \( (N_1,N_2)=(N,\lfloor0.8N\rfloor) \)。  
5. **运行 + 校核**：看 jerk 峰值/RMS、是否越界、\( \Delta T \)。  
6. **若超界**：优先增大 \( N \)；若 \( \Delta T \) 过大则用**优化 FIR**或混合“局部 S 段”。  
7. **稳定后**：微调 \( N_1,N_2 \) 降低 \( \Delta T \) 与中频失真。

---

## 12) 验证与检验方法（通过/不通过判据）

### 12.1 边界与几何
- **必过**：\( s^{(f)}(0)=0,\ \dot s^{(f)}(0)=\dot s(0) \)；
- **必过**：\( s^{(f)}(T_{\mathrm{new}})=s_{\mathrm{end}},\ \dot s^{(f)}(T_{\mathrm{new}})=\dot s(T_0) \)（常为 0）；
- **必过**：\( \dot s^{(f)}(t)\ge 0 \)（单调）；
- **应过**：几何路径未变（仅改时间律）。

### 12.2 约束与 jerk
- **硬指标**：  
  - \( \max_t \dot s^{(f)} \le \dot s_{\max} \)，\( \max_t \ddot s^{(f)} \le \ddot s_{\max} \)；  
  - \( \max_{i,t} |\dddot q_i(t)| \le J_{i,\max} \)（若只在路径域评估，则以 \( J_{\mathrm{path}} \) 为硬阈值）。  
- **软指标**：jerk RMS 与峰值降低比例 ≥ 指定阈值（如 ≥ 50%）。

### 12.3 时长与频域
- 记录 \( \Delta T/T_0 \)；若 > 预设上限（如 30%），需调参/更优 FIR。  
- 频域检查（可选）：FIR 幅频在机构共振频带应提供 ≥ 20 dB 衰减。

### 12.4 可视化与报表（建议自动生成）
- 曲线：\( \dot s,\ \ddot s,\ \dddot s \)（或关节 jerk 热力图）。  
- 统计：峰值、RMS、超界样本率、\(\Delta T\)。  
- CSV 导出：原 vs 平滑的对照数据。

---

## 13) 单元测试与集成测试用例

1. **Sanity**：单轴、单个加速度台阶输入 → jerk 变成平台且受核长控制。  
2. **Endpoints**：随机起止速度（含非零）→ 滤波后仍精确一致。  
3. **Short Segments**：交替短/长恒加速度段 → 自适应核/局部 S 段不超界。  
4. **Multi-Joint Check**：给定 \( \mathbf q'(s),\mathbf q''(s),\mathbf q'''(s) \) 的合成路径 → 关节 jerk 峰值 ≤ 阈值。  
5. **Time Growth**：限制 `max_growth_ratio` → 若不能满足 jerk，则返回“需提高增长上限”的明确状态码。

---

## 14) 失败模式与应对

- **仍有 jerk 超界尖角**：增大核长；在尖角处插入**最小 S 段**（jerk-bounded analytical blend）。  
- **边界被拉歪**：确保使用**基底扣除**或**充足外延**；检查前后向步骤。  
- **速度出现负值**：在积分时投影 \( \dot s\ge 0 \) 并允许时间增大；禁止“硬裁剪”导致几何回退。  
- **\(\Delta T\) 过大**：换“优化 FIR”或“非均匀核长”（只在拐点周围用长核）。  
- **中频失真**：改用双级 MA 或优化 FIR，降低通带起伏。

---

## 15) 结果期望（交付标准）

- 相对 Kunz 原轨迹：
  - **jerk 峰值**显著降低（例如 ↓50–90%），RMS 同步下降；
  - 起止 \( (s,\dot s) \)**严格一致**；
  - **总时长合理增长**（可设 0–30% 目标区间）；
  - \( \dot s,\ddot s \) 无越界；如必须裁剪，越界比率 = 0。

---

## 16) 附：数学与实现细节

### 16.1 三次端点基底系数（示意）
给定 \( b(0)=s_0,\ \dot b(0)=v_0,\ b(T)=s_T,\ \dot b(T)=v_T \)，
\[
b(t)=s_0+v_0 t + a t^2 + c t^3,\quad
a=\frac{3(s_T-s_0)- (2v_0+v_T)T}{T^2},\ 
c=\frac{2(s_0-s_T)+(v_0+v_T)T}{T^3}.
\]

### 16.2 FIR 对阶跃的 jerk 上界（MA 近似）
单级 MA 长度 \( N \) 对加速度台阶 \( \Delta \ddot s \) 的输出斜率为 \( \Delta \ddot s/(N\Delta t) \)，作为离散 jerk 的保守上界。双级 MA 更缓和，取同一公式定尺更保守。

---

## 17) 交付清单（给 Codex）

- 完整 C++ 工程骨架（见 §10）。
- `smooth_totg` CLI（输入 Kunz CSV/JSON，输出平滑 CSV/JSON + 报表）。
- 单元测试与基准脚本（jerk 峰值/RMS/ΔT 对照）。
- 可视化脚本（gnuplot 指令或 matplotlib 数据导出）。

---

## 18) 实施 Checklist

- [ ] 能正确解析 Kunz 输出并等间隔采样  
- [ ] 端点基底构造 & 验证首末 \(s,\dot s\) 不变  
- [ ] FIR 零相位实现（前后向）  
- [ ] 积分重建 + 单调性保证  
- [ ] 限速/限加投影闭环  
- [ ] 时长延展策略（尾段 S 段或全局缩放）  
- [ ] 关节 jerk 评估（链式或差分）  
- [ ] 自动报表（峰值/RMS/ΔT/越界率）  
- [ ] 异常与短段处理  
- [ ] 全部测试用例通过

---

**备注**：以上流程在**不改变路径几何**的前提下，通过**对时间律的 FIR 整形**把加速度阶跃摊平成有限斜率的爬坡段，因此 jerk 被有效限幅；由于过渡被展宽，**总时长合理增加**是**必要现象**。通过“端点基底 + 零相位滤波 + 约束投影 + 时标延展”，可以**严格**满足**始末位置与速度**不变的要求。
