# PSIS-MPIP

**PSIS-MPIP** 是基于 *乘法本原不可约多项式的渐进式秘密图像共享（Progressive Secret Image Sharing using Multiplicative Primitive Irreducible Polynomials）* 的实验实现。  
该项目实现了渐进式秘密图像共享与渐进恢复算法，便于研究人员复现论文结果，同时进行方案的二次开发和算法优化。

---

## 功能特点

- 基于乘法本原不可约多项式的渐进式秘密图像共享（PSIS）方案  
- 多级渐进恢复：可根据分享份数量逐步揭示秘密图像  
- 模块化设计，便于实验和算法扩展  

---

## 目录结构

```text
PSIS-MPIP/
├── PSIS_utils/           # 工具
│   ├── PSISutils.py
│   ├── cal_utils.py
│   └── format_utils.py
├── config.py             # 配置文件
├── config_poly_group.py  # 多项式组配置
├── PSISALL.py            # 核心 PSIS 实现
└── README.md             # 项目说明
```

---

## 项目声明

- 本项目的作者及单位：
- The author and affiliation of this project:
- 
```text
- 项目名称（Project Name）：PSIS-MPIP
- 项目作者（Author）：Xinyue Wang, Xiaotian Wu
- 作者单位（Affiliation）：暨南大学信息科学技术学院（School of Information Science and Technology, Jinan University）
```

