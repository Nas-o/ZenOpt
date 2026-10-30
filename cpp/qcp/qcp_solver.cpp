/**
 * @file    qcp_solver.cpp
 * @brief   二维「线性约束 + 单个凸二次约束」规划问题解析求解器实现
 */
#include <cmath>
#include "qcp_solver.hpp"

namespace robot::algorithm
{

void Polygon2f::clear()
{
    count = 0;
}

void Polygon2f::add(const Vec2f &pt)
{
    if (count < kQcpMaxPoints)
        pts[count++] = pt;
}

float QcpSolver2D::evalConic(const Conic2f &conic, const Vec2f &x)
{
    return 0.5f * (conic.q00 * x.x * x.x + 2.0f * conic.q01 * x.x * x.y + conic.q11 * x.y * x.y)
         + conic.p0 * x.x + conic.p1 * x.y + conic.r;
}

void QcpSolver2D::clipHalfPlane(Polygon2f &poly, float nx, float ny, float c)
{
    if (poly.count <= 0)
        return;

    Polygon2f out;
    out.clear();

    if (poly.count == 1)
    {
        const Vec2f &p = poly.pts[0];
        if (nx * p.x + ny * p.y <= c)
            out.add(p);
        poly = out;
        return;
    }

    for (int i = 0; i < poly.count; i++)
    {
        const Vec2f &cur  = poly.pts[i];
        const Vec2f &prev = poly.pts[(i + poly.count - 1) % poly.count];

        const float d_cur  = nx * cur.x  + ny * cur.y  - c;
        const float d_prev = nx * prev.x + ny * prev.y - c;

        if (d_cur <= 0.0f)
        {
            if (d_prev > 0.0f)
            {
                const float t = d_prev / (d_prev - d_cur);
                out.add(Vec2f{prev.x + t * (cur.x - prev.x), prev.y + t * (cur.y - prev.y)});
            }
            out.add(cur);
        }
        else if (d_prev <= 0.0f)
        {
            const float t = d_prev / (d_prev - d_cur);
            out.add(Vec2f{prev.x + t * (cur.x - prev.x), prev.y + t * (cur.y - prev.y)});
        }
    }

    poly = out;
}

bool QcpSolver2D::pointInPolygon(const Polygon2f &poly, const Vec2f &pt, float tol)
{
    if (poly.count < 3)
        return false;

    int sign = 0;
    for (int i = 0; i < poly.count; i++)
    {
        const Vec2f &p1 = poly.pts[i];
        const Vec2f &p2 = poly.pts[(i + 1) % poly.count];

        const float cross = (p2.x - p1.x) * (pt.y - p1.y) - (p2.y - p1.y) * (pt.x - p1.x);
        if (fabsf(cross) <= tol)
            continue;                     

        const int s = (cross > 0.0f) ? 1 : -1;
        if (sign == 0)
            sign = s;
        else if (s != sign)
            return false;
    }
    return true;
}

Vec2f QcpSolver2D::solve(const Polygon2f &poly, const Conic2f &conic, const Vec2f &a, bool &feasible)
{
    feasible = false;

    if (poly.count <= 0)
        return Vec2f{};

    Vec2f best{};
    float best_g = 1e30f;
    for (int i = 0; i < poly.count; i++)
    {
        const float g = evalConic(conic, poly.pts[i]);
        if (g < best_g)
        {
            best_g = g;
            best   = poly.pts[i];
        }
    }

    /* 不可行时返回边界上 g(x) 最小的点，线段退化时同样有效 */
    if (poly.count >= 2)
    {
        for (int i = 0; i < poly.count; i++)
        {
            const Vec2f &x1 = poly.pts[i];
            const Vec2f &x2 = poly.pts[(i + 1) % poly.count];
            const float dx = x2.x - x1.x;
            const float dy = x2.y - x1.y;
            const float qa = 0.5f * (conic.q00 * dx * dx + 2.0f * conic.q01 * dx * dy + conic.q11 * dy * dy);
            const float qb = conic.q00 * x1.x * dx + conic.q01 * (x1.x * dy + x1.y * dx) + conic.q11 * x1.y * dy
                           + conic.p0 * dx + conic.p1 * dy;

            if (qa > 1e-9f)
            {
                const float edge_t = -qb / (2.0f * qa);
                if (edge_t > 0.0f && edge_t < 1.0f)
                {
                    const Vec2f edge_min{x1.x + edge_t * dx, x1.y + edge_t * dy};
                    const float g = evalConic(conic, edge_min);
                    if (g < best_g)
                    {
                        best_g = g;
                        best = edge_min;
                    }
                }
            }
        }
    }

    Conic2f analytic_conic = conic;
    float det = analytic_conic.q00 * analytic_conic.q11 -
                analytic_conic.q01 * analytic_conic.q01;
    const float det_scale = fabsf(analytic_conic.q00 * analytic_conic.q11) +
                            analytic_conic.q01 * analytic_conic.q01;
    const bool needs_regularization = analytic_conic.q00 <= 0.0f ||
                                      analytic_conic.q11 <= 0.0f ||
                                      det <= 1e-6f * fmaxf(det_scale, 1e-12f);
    if (needs_regularization)
    {
        const float q_scale = fmaxf(fmaxf(fabsf(analytic_conic.q00),
                                          fabsf(analytic_conic.q11)),
                                    fabsf(analytic_conic.q01));
        const float epsilon = fmaxf(1e-6f * q_scale, 1e-9f);
        analytic_conic.q00 += epsilon;
        analytic_conic.q11 += epsilon;
        det = analytic_conic.q00 * analytic_conic.q11 -
              analytic_conic.q01 * analytic_conic.q01;
    }

    const bool invertible = analytic_conic.q00 > 0.0f &&
                            analytic_conic.q11 > 0.0f && det > 0.0f;
    float area_twice = 0.0f;
    for (int i = 0; i < poly.count; i++)
    {
        const Vec2f &p1 = poly.pts[i];
        const Vec2f &p2 = poly.pts[(i + 1) % poly.count];
        area_twice += p1.x * p2.y - p2.x * p1.y;
    }

    if (invertible && fabsf(area_twice) > 1e-8f)
    {
        const float inv00 =  analytic_conic.q11 / det;
        const float inv01 = -analytic_conic.q01 / det;
        const float inv11 =  analytic_conic.q00 / det;

        const float s = a.x * (inv00 * a.x + inv01 * a.y)
                      + a.y * (inv01 * a.x + inv11 * a.y);
        const float t = analytic_conic.p0 * (inv00 * analytic_conic.p0 + inv01 * analytic_conic.p1)
                      + analytic_conic.p1 * (inv01 * analytic_conic.p0 + inv11 * analytic_conic.p1)
                      - 2.0f * analytic_conic.r;

        if (s > 1e-9f && t > 1e-9f)
        {
            const float lambda = sqrtf(s / t);

            Vec2f tangent;
            tangent.x = inv00 * (a.x / lambda - analytic_conic.p0) +
                        inv01 * (a.y / lambda - analytic_conic.p1);
            tangent.y = inv01 * (a.x / lambda - analytic_conic.p0) +
                        inv11 * (a.y / lambda - analytic_conic.p1);

            /* 正则化切点相对原二次约束略偏向内部，使用原约束检查可行性 */
            if (evalConic(conic, tangent) <= 1e-4f && pointInPolygon(poly, tangent))
            {
                feasible = true;
                return tangent;
            }
        }
    }

    /* ---- 2. 二次约束裁剪多边形: 枚举顶点与交点 ---- */
    bool  found    = false;
    float best_obj = 0.0f;

    for (int i = 0; i < poly.count; i++)
    {
        if (evalConic(conic, poly.pts[i]) > 1e-4f)
            continue;

        const float obj = a.x * poly.pts[i].x + a.y * poly.pts[i].y;
        if (!found || obj > best_obj)
        {
            best_obj = obj;
            best     = poly.pts[i];
            found    = true;
        }
    }

    for (int i = 0; i < poly.count; i++)
    {
        const Vec2f &x1 = poly.pts[i];
        const Vec2f &x2 = poly.pts[(i + 1) % poly.count];

        const float dx = x2.x - x1.x;
        const float dy = x2.y - x1.y;

        const float qa = 0.5f * (conic.q00 * dx * dx + 2.0f * conic.q01 * dx * dy + conic.q11 * dy * dy);
        const float qb = conic.q00 * x1.x * dx + conic.q01 * (x1.x * dy + x1.y * dx) + conic.q11 * x1.y * dy
                       + conic.p0 * dx + conic.p1 * dy;
        const float qc = evalConic(conic, x1);

        float roots[2];
        int   n_roots = 0;

        if (fabsf(qa) < 1e-9f)
        {
            if (fabsf(qb) > 1e-9f)
                roots[n_roots++] = -qc / qb;
        }
        else
        {
            const float delta = qb * qb - 4.0f * qa * qc;
            if (delta >= -1e-6f)
            {
                const float sq = sqrtf(fmaxf(delta, 0.0f));
                roots[n_roots++] = (-qb + sq) / (2.0f * qa);
                roots[n_roots++] = (-qb - sq) / (2.0f * qa);
            }
        }

        for (int k = 0; k < n_roots; k++)
        {
            if (roots[k] < -1e-5f || roots[k] > 1.0f + 1e-5f)
                continue;

            const float root = fminf(fmaxf(roots[k], 0.0f), 1.0f);
            const Vec2f pt{x1.x + root * dx, x1.y + root * dy};
            const float obj = a.x * pt.x + a.y * pt.y;
            if (!found || obj > best_obj)
            {
                best_obj = obj;
                best     = pt;
                found    = true;
            }
        }
    }

    feasible = found;
    return best;
}

} // namespace robot::algorithm