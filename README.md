# 古树倒伏风险模型 / Wind–root anchorage safety factor

一个最小可复现的风荷载—根系锚固力学模型。用来回答两个问题：

1. 扛过数百年台风的古树，为什么会在一场暴雨里连根倒？
2. 有限的人力与经费，现场应该先量哪几个参数？

**不是**倒伏概率预测，**不是**单株安全性判定工具。

---

## 背景

2026 年 9 月 3 日，台风"沙德尔"带来的持续强降雨使福建南靖云水谣河水暴涨，
一株树龄 624 年的古榕（一级保护古树，与相邻 674 年古榕合称"夫妻树"）在树基被冲刷、
根系失去抓地力后整株连根倒伏，9 月 9 日晚由两台 300 吨吊车扶正。

它扛过了 624 年的风，却没能扛过被水掏空的树基。这提示：决定古树倒伏的不是风速，
而是根系与土体状态。本项目据此建模，并检验这个判断在数量上是否成立。

## 模型

```
风压        q      = 0.5 · ρ · v²                      ρ = 1.225 kg/m³
迎风面积    A      = 0.65 · (0.75H) · W                0.65 为冠形折减
倾覆弯矩    M_wind = Cd · q · A · (0.6H) · (1 + 偏冠率)
抗力弯矩    M_root = 0.3 · su · D · h_eff²             h_eff = max(h − d, 0)
安全系数    SF     = M_root / M_wind
```

- `H` 树高，`W` 冠幅，`Cd` 冠层阻力系数
- `su` 土体不排水抗剪强度，`D` 根盘直径，`h` 根盘埋深，`d` 树基冲刷深度
- 抗力式是根盘侧向抗力的简化表达，系数 `0.3` 取文献量级（0.25–0.5）中值

## 参数

| 参数 | 取值 | 来源 |
|---|---|---|
| 基本风压 `w0` | 0.70 kN/m²（福州，50 年重现期） | GB 50009 量级 |
| 极值风速 `v` | Gumbel，50 年重现期 33.8 m/s，CV 0.15，截断 12–60 | 由 `w0` 反算 |
| 冠层阻力系数 `Cd` | 0.15–0.30 均匀 | 文献量级 |
| 根盘直径 `D` | 8–16 m 均匀 | 量级假设 |
| 根盘埋深 `h` | 1.5–3.5 m 均匀 | 量级假设 |
| 土体抗剪强度 `su` | 25–80 kPa 均匀 | 一般黏性土 |
| 树基冲刷深度 `d` | 80% 在 0–0.5 m，20% 在 0.5–3.0 m | 现场情形假设 |
| 基准株 | H = 22 m，冠幅 17 m，偏冠率 0.15 | 参照"总重 100 多吨"的古榕 |

**全部参数来自公开规范或文献量级，无一项为实测标定值。**

## 结果

| 指标 | 值 |
|---|---|
| 健康株、设计风速 33.8 m/s、无冲刷 | **SF = 2.71** |
| 各参数取中位 | SF ≈ 6.0 |
| 临界冲刷深度（SF 跌至 1.0） | **1.62 m** |
| 20000 组工况 | SF 中位 4.5，5% 分位 0.002，95% 分位 19.3 |

Sobol 全局敏感性（Saltelli 估计量，16384 次模型调用）：

| 参数 | S1 | ST |
|---|---|---|
| 根盘埋深 | 0.233 | **0.385** |
| 极值风速 | 0.153 | **0.284** |
| 树基冲刷深度 | 0.159 | **0.236** |
| 土体抗剪强度 | 0.107 | 0.162 |
| 冠层阻力系数 | 0.043 | 0.074 |
| 根盘直径 | 0.038 | 0.073 |

根系与土体四项的一阶效应合计 **0.54**，约为风速（0.15）的 3.5 倍。

**结论**：在这个模型里，决定安不安全的主要是根和土，不是风。所以现场调查应优先测
根盘埋深与树基冲刷深度，而不是树高冠幅。

## 图

| 文件 | 内容 |
|---|---|
| `figures/fig1_sf_cdf.png` | 20000 组工况的安全系数累积分布 |
| `figures/fig2_scour_vs_sf.png` | 安全系数随树基冲刷深度的变化，标出临界深度 |
| `figures/fig3_sobol.png` | Sobol 一阶与总效应指数 |

`figures/metrics.json` 是代码输出的全部数字，README 与图中每个数都能在这里对上。

## 运行

```bash
pip install -r requirements.txt
python src/mech_model.py
```

约 3 秒。固定随机种子（Monte Carlo 42，Sobol 7），结果逐位可复现。

## 边界（请连同结果一起引用）

1. 参数未做现场标定，模型未做倾覆试验验证；
2. 阻力式与抗力式均为简化式，只用于排序与敏感性分析，不用于单株安全性判定；
3. "13.7% 工况低于 1.0" 描述的是给定参数范围下的抽样分布，**不是真实古树的倒伏比例**；
4. 正式对外使用前，须以根盘探测、土工试验与（有条件时）倾覆试验标定。

## 后续

- 用无人机与地基激光雷达点云提取 `H`、`W`、偏冠率、树冠投影（当前为估读值）
- 以根盘探测与土工试验标定 `h`、`D`、`su`
- 把单株安全系数接入分级清单与巡检优先级排序

## 许可

MIT。若用于研究请注明出处。

---

## English abstract

A minimal model of wind-induced overturning of ancient trees against root-plate
anchorage. Wind pressure is derived from the 50-year basic wind pressure of
Fuzhou; extreme wind speed follows a Gumbel distribution; root resistance is a
simplified lateral-capacity expression. Monte Carlo propagation and Saltelli
Sobol indices show that root-plate embedment depth and scour depth dominate the
variance, while wind speed ranks second — consistent with the 624-year-old
banyan that fell after its root plate was scoured away by a flooded river in
September 2026. All parameters are order-of-magnitude values from codes and the
literature; **none is field-calibrated**. The code is intended for ranking and
sensitivity analysis only.
