# FuturePark

run MGT.m

towerC 柱子已完成。

梁及斜撑 ELEMENT 已完成。
板 ELEMENT 已完成。

边界 已完成。

板面荷载

自重、风

地震

荷载组合


stair 程序有问题。

2018年4月22日
两个角塔中心定位待定。
各个筒的角度未加入。

2018年4月24日
电梯筒的楼梯部分，未加载楼梯荷载(考虑加梁上荷载)。电梯荷载未加载。

2018年4月26日
幕墙定位目前是根据剖面的点进行定位，如层高修改或层数修改则需要同时修改幕墙定位点坐标。
Misc. CAD里样条曲线的函数如何得到

2018年4月27日
CAD里样条曲线，是NURBS(默认三阶)算法。检索NURBS的矩阵表达式。

2018年5月4日
Daniel的曲线是所有点权重均为1的，3阶(degree)，spline。

2018年5月6日
nurbsfunctions里节点向量设为t = [0 0 0 0 0.5 1 1 1 1]时与Daniel的曲线一致。

2018年5月9日
model.m里层高数据更新，另外的m文件未更新，facade数据未与CAD同步。

2018年5月20日
coorLxC_sym为通过符号运算求解直线与圆交点，耗时稍长。
可以从coorLxCp入手，升级，避免符号运算。
角塔未完成。

2018年5月21日
角塔已建好。

2018年5月28日
楼梯荷载，改恒3活2。(因加在梯板上，钢楼梯，且两跑楼梯，故每层6/4，)
所有楼面荷载考虑改为压力荷载。

2018年6月19日
3个停车塔楼，全部不转角度。(考虑停车AGV方便)

2018年8月1日
圆心及角点坐标updated
MGT_facade_S2.m 未改 <=

2018年8月7日
用coorLxCp 替换了 coorLxC_sym
外边线已建好
精细化建模 可 以直代曲

2018年8月15日
以直代曲 完成