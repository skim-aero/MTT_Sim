# MTT_Sim: Simple Multi-target Tracking Simulation
## 1. Software Introduction
This is a simple simulation code for a multi-target tracking (MTT) example with clutter and the comparison of MTT algorithms including global nearest neighbout (GNN), kernel GNN, join probability data association (JPDA) filter and kernel JPDA filter. Only two targets in the noisy measurement (with clutter) are considered at the moment.

The kernel GNN and kernel JPDA map the state into the kernel space and the data association is done in the kernel space. Currently kernel JPDA doesn't work properly, need update.

For the detail of the JPDA, please refer to the paper follows:
[AFJPDA]](https://arc.aiaa.org/doi/full/10.2514/1.I011301)

**Coded by Sukkeun Kim, Cranfield University, UK**
* Email: <s.kim.aero@gmail.com>

## 2. Version Information

* 12th Jun 2023 Beta 0.0: First commit
* 11th Nov 2024 Release 1.0: First public
* 12th Nov 2024 Release 1.1: Minor modifications
*  4th Dec 2024 Release 1.2: Equation update
