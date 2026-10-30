"""运行：python examples/python/qcp/example_qcp.py"""
from pathlib import Path
import sys

# 支持从任意工作目录直接运行，无需安装项目。
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "python"))

from qcp import Conic2f, Polygon2f, QcpSolver2D, Vec2f


def main():
    poly = Polygon2f([Vec2f(-2, -2), Vec2f(2, -2), Vec2f(2, 2), Vec2f(-2, 2)])
    # x <= 0.6；二次约束 x² + y² <= 1；目标 max(x + y)。
    QcpSolver2D.clipHalfPlane(poly, 1.0, 0.0, 0.6)
    conic = Conic2f(q00=2.0, q11=2.0, r=-1.0)
    point, feasible = QcpSolver2D.solve(poly, conic, Vec2f(1.0, 1.0))
    print(f"feasible={feasible}, x={point.x:.6f}, y={point.y:.6f}")
    print(f"objective={point.x + point.y:.6f}")


if __name__ == '__main__':
    main()
