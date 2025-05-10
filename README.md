﻿# RCWA 求解器简介

严格耦合波分析 （RCWA） 求解器可用于分析入射到多层结构上的平面波的光学响应。
RCWA 求解器可用于在层几何结构中具有周期性变化的结构，例如光子晶体和衍射光栅。
RCWA 求解器的仿真时间通常比 FDTD 短得多，因此是分析这些类型的周期性结构的理想工具。

该项目基于C++编写，使用数学库Armadillo实现矩阵运算，使用FFTW3计算二维卷积矩阵。
器件结构和仿真所需参数的输入文件在Input目录，仿真结果的输出文件在Output目录。

## 原理
RCWA 方法是一种半解析技术，用于求解多层结构中的麦克斯韦方程组。在这种方法中，结构沿传播方向被分成一系列均匀的层。
沿传播方向具有逐渐变化的横截面的结构可以用一系列均匀的层来近似。例如，在下面显示的几何结构中，梯形形状 近似为一系列六层 ：

![](examples1/images/RCWA_discrete_layers.png)

 增加截面上的层数可以提高模拟的精度，但会增加模拟时间。

 将结构划分为多个层后，麦克斯韦方程组将在傅里叶域的每一层中解析求解。傅里叶模式的波矢称为 k 矢量。由于结构的周期性，只允许使用离散的 k 向量。增加 k 向量的数量可以提高精度，但代价是仿真时间增加。

 然后，每个部分的解双向传播，以计算整个器件的 S 矩阵。计算出 S 矩阵后，来自入射平面波的光就可以传播到结构中。由于几何形状的周期性，入射平面波被衍射成一组有限的平面波，称为“光栅级”。
 计算 S 矩阵后，可以计算出传输和反射的入射功率分数、每个光栅阶数的功率以及结构内部的电场和磁场等结果。
## 建模

![](examples1/images/RCWA_model.png)

Note: 平面波的传播轴为z轴，沿着z轴正向传播。



关于Input文件夹的数据介绍

| 名称 | 含义|
|---|---|
| ku | x方向最大谐波  | 
| kv | y方向最大谐波  |
 | lambda | 波长  |
 | LayerPos| 每一层中心z坐标|
 | n_lower  | 入射区域的折射率|
|   n_upper  |反射区域折射率|
| phi  | 入射平面波旋转角|
| theta | 入射角| 
| x  |横向x网格坐标|
| y |横向y网格坐标|
| z | z向分层坐标|
|Index_real_z_i_j | 第i层第j个波长的实折射率|
|Index_imag_z_i_j| 第i层第j个波长的虚折射率|


 

## 使用

使用的软件为Visual Studio 2022 ， 库环境安装如下。

```
tree /f > filename 
git clone https://github.com/microsoft/vcpkg.git
.\bootstrap-vcpkg.bat 
vcpkg install armadillo
vcpkg install fftw3
vcpkg integrate install
```

## 案例

仿真案例采用如图的器件结构 ，

![](examples1/images/struct.png)

RCWA的仿真模型为,

![](examples1/images/simulation.png)

仿真参数设置

| 名称 | 值 |
| --- | --- |
| ku | 4 |
| kv | 4 |
| theta | 45 |
| phi | 45 |
| lambda | 1/0.55 ~ 1/0.5 um |

计算的结果如下图，其中R，T代表反射率和透射率，s,p分别表示平面波垂直极化和平行极化。
实线的结果是Lumerical RCWA，虚线的结果是本项目的RCWA求解器的结果，表明代码的有效性。
![](examples1/images/RT.png)

## 总结

觉得有意思的话点个赞🥰 O(∩_∩)O 