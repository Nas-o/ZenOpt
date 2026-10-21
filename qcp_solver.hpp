
/**
 * @file    qcp_solver.hpp
 * @brief   二维「线性约束 + 单个凸二次约束」规划问题解析求解器
 * @details 问题形式:
 *
 *              maximize    a^T·x
 *              subject to  x ∈ P                        (凸多边形, 由线性约束裁剪得到)
 *                          0.5·x^T·Qx + p^T·x + r ≤ 0       (Q 半正定, 凸二次约束)
 *
 *          求解步骤 :
 *            1) 拉格朗日乘子法求二次约束边界上使线性目标最大的切点 ;
 *            2) 切点落在多边形内 -> 该点即全局最优解, 直接返回;
 *            3) 否则用二次约束裁剪多边形, 枚举"多边形顶点(在约束内)"与
 *               "多边形边与二次曲线的交点", 取目标函数最大者。
 *
 *          凸规划问题的局部最优解即全局最优, 而线性目标在凸集上的最优解
 *          一定取在极点上, 故步骤 3 的枚举是完备的。
 *
 *          实现约束: 无动态内存、无异常、定点数量的浮点运算, 可用于 1kHz 控制循环。
 */
#ifndef __QCP_SOLVER_H__
#define __QCP_SOLVER_H__

#include <cstdint>

namespace robot::algorithm
{

/* ============================ 基础类型 ============================ */

/** @brief 二维点/向量 */
struct Vec2f
{
    float x = 0.0f;
    float y = 0.0f;
};

/** @brief 凸二次约束: 0.5·x^T·Qx + p^T·x + r ≤ 0, Q 为对称阵 */
struct Conic2f
{
    float q00 = 0.0f;   /**< Q(0,0) */
    float q01 = 0.0f;   /**< Q(0,1) = Q(1,0) */
    float q11 = 0.0f;   /**< Q(1,1) */
    float p0  = 0.0f;   /**< p(0) */
    float p1  = 0.0f;   /**< p(1) */
    float r   = 0.0f;   /**< 常数项 */
};

/** @brief 凸多边形顶点容量: 初始 4 个顶点, 每裁剪一个半平面最多增加 1 个 */
constexpr int kQcpMaxPoints = 16;

/** @brief 凸多边形 (顶点按同一绕向排列) */
struct Polygon2f
{
    Vec2f pts[kQcpMaxPoints]{};
    int   count = 0;

    void clear();
    void add(const Vec2f &pt);
};

/* ============================ 求解器 ============================ */

class QcpSolver2D
{
public:
    /* ---------------- 基础工具 ---------------- */

    /** @brief 二次约束左端项 g(x) = 0.5·x^T·Qx + p^T·x + r */
    static float evalConic(const Conic2f &conic, const Vec2f &x);

    /**
     * @brief  用半平面 nx·x + ny·y ≤ c 裁剪凸多边形 (Sutherland-Hodgman, 伪代码见文档 5.3.6.2.1)
     * @param  poly 输入输出, 顶点按序排列
     */
    static void clipHalfPlane(Polygon2f &poly, float nx, float ny, float c);

    /** @brief 判断点是否在凸多边形内 (含边界), 不要求知道绕向 */
    static bool pointInPolygon(const Polygon2f &poly, const Vec2f &pt, float tol = 1e-5f);

    /* ---------------- 主求解入口 ---------------- */

    /**
     * @brief  在凸多边形内最大化 a^T·x, 同时满足凸二次约束
     * @param  poly     线性约束构成的凸多边形
     * @param  conic    凸二次约束
     * @param  a        目标函数方向
     * @param  feasible 输出: 返回点是否满足二次约束
     * @return 最优解; 当二次约束不可满足时返回可行域内使 g(x) 最小的点(尽力而为)
     */
    static Vec2f solve(const Polygon2f &poly, const Conic2f &conic, const Vec2f &a, bool &feasible);
};

} // namespace robot::algorithm

#endif /* __QCP_SOLVER_H__ */