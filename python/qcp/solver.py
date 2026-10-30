"""二维线性目标 + 凸多边形 + 单个凸二次约束的解析求解器。

移植自 qcp_solver.hpp / qcp_solver.cpp，保留原算法及数值阈值。
Q 应为半正定，顶点须按一致绕向排列。Python float 为双精度。
"""
from dataclasses import dataclass, field, replace
from math import sqrt
from typing import List, Tuple

kQcpMaxPoints = 16


@dataclass
class Vec2f:
    x: float = 0.0
    y: float = 0.0


@dataclass
class Conic2f:
    """约束：0.5 * x^T Q x + p^T x + r <= 0。"""
    q00: float = 0.0
    q01: float = 0.0
    q11: float = 0.0
    p0: float = 0.0
    p1: float = 0.0
    r: float = 0.0


@dataclass
class Polygon2f:
    pts: List[Vec2f] = field(default_factory=list)

    def __post_init__(self) -> None:
        self.pts = [replace(pt) for pt in self.pts[:kQcpMaxPoints]]

    @property
    def count(self) -> int:
        return len(self.pts)

    def clear(self) -> None:
        self.pts.clear()

    def add(self, pt: Vec2f) -> None:
        if self.count < kQcpMaxPoints:
            self.pts.append(replace(pt))


class QcpSolver2D:
    @staticmethod
    def evalConic(conic: Conic2f, x: Vec2f) -> float:
        return (0.5 * (conic.q00 * x.x * x.x + 2.0 * conic.q01 * x.x * x.y
                       + conic.q11 * x.y * x.y)
                + conic.p0 * x.x + conic.p1 * x.y + conic.r)

    @staticmethod
    def clipHalfPlane(poly: Polygon2f, nx: float, ny: float, c: float) -> None:
        """原地用 nx*x + ny*y <= c 裁剪多边形。"""
        if poly.count <= 0:
            return
        out = Polygon2f()
        if poly.count == 1:
            p = poly.pts[0]
            if nx * p.x + ny * p.y <= c:
                out.add(p)
        else:
            for i, cur in enumerate(poly.pts):
                prev = poly.pts[i - 1]
                d_cur = nx * cur.x + ny * cur.y - c
                d_prev = nx * prev.x + ny * prev.y - c
                if d_cur <= 0.0:
                    if d_prev > 0.0:
                        t = d_prev / (d_prev - d_cur)
                        out.add(Vec2f(prev.x + t * (cur.x - prev.x),
                                      prev.y + t * (cur.y - prev.y)))
                    out.add(cur)
                elif d_prev <= 0.0:
                    t = d_prev / (d_prev - d_cur)
                    out.add(Vec2f(prev.x + t * (cur.x - prev.x),
                                  prev.y + t * (cur.y - prev.y)))
        poly.pts = out.pts

    @staticmethod
    def pointInPolygon(poly: Polygon2f, pt: Vec2f, tol: float = 1e-5) -> bool:
        if poly.count < 3:
            return False
        sign = 0
        for i, p1 in enumerate(poly.pts):
            p2 = poly.pts[(i + 1) % poly.count]
            cross = (p2.x - p1.x) * (pt.y - p1.y) - (p2.y - p1.y) * (pt.x - p1.x)
            if abs(cross) <= tol:
                continue
            s = 1 if cross > 0.0 else -1
            if sign == 0:
                sign = s
            elif s != sign:
                return False
        return True

    @staticmethod
    def _edge_coefficients(conic: Conic2f, x1: Vec2f, dx: float, dy: float) -> Tuple[float, float]:
        qa = 0.5 * (conic.q00 * dx * dx + 2.0 * conic.q01 * dx * dy + conic.q11 * dy * dy)
        qb = (conic.q00 * x1.x * dx + conic.q01 * (x1.x * dy + x1.y * dx)
              + conic.q11 * x1.y * dy + conic.p0 * dx + conic.p1 * dy)
        return qa, qb

    @staticmethod
    def solve(poly: Polygon2f, conic: Conic2f, a: Vec2f) -> Tuple[Vec2f, bool]:
        """最大化 a^T*x，返回 (点, feasible)。

        保留 C++ 的尽力返回行为：无候选解时返回所检查顶点/边上 g 最小
        的点；这不保证是整个多边形内的 g 全局最小点。空多边形返回原点、False。
        正则化、零目标及退化约束的处理也沿用原实现，未额外修正。
        """
        if poly.count <= 0:
            return Vec2f(), False
        best = Vec2f()
        best_g = 1e30
        for pt in poly.pts:
            g = QcpSolver2D.evalConic(conic, pt)
            if g < best_g:
                best_g, best = g, pt
        if poly.count >= 2:
            for i, x1 in enumerate(poly.pts):
                x2 = poly.pts[(i + 1) % poly.count]
                dx, dy = x2.x - x1.x, x2.y - x1.y
                qa, qb = QcpSolver2D._edge_coefficients(conic, x1, dx, dy)
                if qa > 1e-9:
                    edge_t = -qb / (2.0 * qa)
                    if 0.0 < edge_t < 1.0:
                        edge_min = Vec2f(x1.x + edge_t * dx, x1.y + edge_t * dy)
                        g = QcpSolver2D.evalConic(conic, edge_min)
                        if g < best_g:
                            best_g, best = g, edge_min

        ac = replace(conic)
        det = ac.q00 * ac.q11 - ac.q01 * ac.q01
        det_scale = abs(ac.q00 * ac.q11) + ac.q01 * ac.q01
        if ac.q00 <= 0.0 or ac.q11 <= 0.0 or det <= 1e-6 * max(det_scale, 1e-12):
            q_scale = max(abs(ac.q00), abs(ac.q11), abs(ac.q01))
            epsilon = max(1e-6 * q_scale, 1e-9)
            ac.q00 += epsilon
            ac.q11 += epsilon
            det = ac.q00 * ac.q11 - ac.q01 * ac.q01
        invertible = ac.q00 > 0.0 and ac.q11 > 0.0 and det > 0.0
        area_twice = sum(p1.x * poly.pts[(i + 1) % poly.count].y
                         - poly.pts[(i + 1) % poly.count].x * p1.y
                         for i, p1 in enumerate(poly.pts))
        if invertible and abs(area_twice) > 1e-8:
            inv00, inv01, inv11 = ac.q11 / det, -ac.q01 / det, ac.q00 / det
            s = a.x * (inv00 * a.x + inv01 * a.y) + a.y * (inv01 * a.x + inv11 * a.y)
            t = (ac.p0 * (inv00 * ac.p0 + inv01 * ac.p1)
                 + ac.p1 * (inv01 * ac.p0 + inv11 * ac.p1) - 2.0 * ac.r)
            if s > 1e-9 and t > 1e-9:
                lam = sqrt(s / t)
                tangent = Vec2f(inv00 * (a.x / lam - ac.p0) + inv01 * (a.y / lam - ac.p1),
                               inv01 * (a.x / lam - ac.p0) + inv11 * (a.y / lam - ac.p1))
                if (QcpSolver2D.evalConic(conic, tangent) <= 1e-4
                        and QcpSolver2D.pointInPolygon(poly, tangent)):
                    return tangent, True

        found, best_obj = False, 0.0
        for pt in poly.pts:
            if QcpSolver2D.evalConic(conic, pt) > 1e-4:
                continue
            obj = a.x * pt.x + a.y * pt.y
            if not found or obj > best_obj:
                best_obj, best, found = obj, pt, True
        for i, x1 in enumerate(poly.pts):
            x2 = poly.pts[(i + 1) % poly.count]
            dx, dy = x2.x - x1.x, x2.y - x1.y
            qa, qb = QcpSolver2D._edge_coefficients(conic, x1, dx, dy)
            qc = QcpSolver2D.evalConic(conic, x1)
            roots = []
            if abs(qa) < 1e-9:
                if abs(qb) > 1e-9:
                    roots.append(-qc / qb)
            else:
                delta = qb * qb - 4.0 * qa * qc
                if delta >= -1e-6:
                    sq = sqrt(max(delta, 0.0))
                    roots = [(-qb + sq) / (2.0 * qa), (-qb - sq) / (2.0 * qa)]
            for root in roots:
                if root < -1e-5 or root > 1.0 + 1e-5:
                    continue
                root = min(max(root, 0.0), 1.0)
                pt = Vec2f(x1.x + root * dx, x1.y + root * dy)
                obj = a.x * pt.x + a.y * pt.y
                if not found or obj > best_obj:
                    best_obj, best, found = obj, pt, True
        return replace(best), found
