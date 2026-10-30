"""QCP 求解性能、延迟和持续压力测试（仅使用标准库）。

快速运行：python benchmarks/python/qcp/benchmark_qcp.py
较长压测：python benchmarks/python/qcp/benchmark_qcp.py --seconds 2 --stress-seconds 60 --samples 10000 --random-cases 10000 --output reports/benchmarks/python/qcp/benchmark_results.json
计时对象为 solve()；输入构造和正确性检查均在计时之外。
"""
import argparse
from dataclasses import asdict, dataclass
import gc
import json
import math
from pathlib import Path
import platform
import random
import statistics
import sys
import time

# 支持从任意工作目录直接运行，无需安装项目。
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "python"))

from qcp import Conic2f, Polygon2f, QcpSolver2D, Vec2f, kQcpMaxPoints


@dataclass
class Case:
    name: str
    poly: Polygon2f
    conic: Conic2f
    a: Vec2f
    expected_feasible: bool = True
    expected_objective: float = None

    def solve(self):
        return QcpSolver2D.solve(self.poly, self.conic, self.a)


def polygon(n, radius=1.0):
    return Polygon2f([Vec2f(radius * math.cos(2 * math.pi * i / n),
                           radius * math.sin(2 * math.pi * i / n)) for i in range(n)])


def cases():
    disk = Conic2f(q00=2.0, q11=2.0, r=-1.0)
    box = Polygon2f([Vec2f(-2, -2), Vec2f(2, -2), Vec2f(2, 2), Vec2f(-2, 2)])
    clipped = Polygon2f(box.pts)
    QcpSolver2D.clipHalfPlane(clipped, 1, 0, .6)
    result = [
        Case('tangent_4', box, disk, Vec2f(1, 1), True, math.sqrt(2)),
        Case('edge_intersection_4', clipped, disk, Vec2f(1, 1), True, 1.4),
        Case('all_feasible_16', polygon(16), Conic2f(r=-1), Vec2f(1, 0), True, 1),
        Case('infeasible_16', polygon(16), Conic2f(q00=2, q11=2, r=1), Vec2f(1, 0), False),
        Case('near_singular_16', polygon(16), Conic2f(q00=2, q11=1e-12, r=-.25), Vec2f(0, 1), True, 1),
        Case('linear_constraint_16', polygon(16), Conic2f(p0=1, r=-.2), Vec2f(1, 0), True, .2),
        Case('segment', Polygon2f([Vec2f(-2, 0), Vec2f(2, 0)]), disk, Vec2f(1, 0), True, 1),
        Case('empty', Polygon2f(), disk, Vec2f(1, 0), False),
    ]
    # 半径介于正多边形内切圆和外接圆之间：每条边有两个交点。
    # 目标朝边的中点，圆的无约束切点位于多边形之外，迫使完整枚举。
    for n in (4, 8, kQcpMaxPoints):
        radius = (1 + math.cos(math.pi / n)) / 2
        result.append(Case(f'all_edges_intersect_{n}', polygon(n),
                           Conic2f(q00=2, q11=2, r=-radius**2),
                           Vec2f(math.cos(math.pi / n), math.sin(math.pi / n)),
                           True, math.cos(math.pi / n)))
    return result


def validate(case, check_expected=True):
    pt, feasible = case.solve()
    if not (math.isfinite(pt.x) and math.isfinite(pt.y)):
        return 'nonfinite result'
    if check_expected and feasible != case.expected_feasible:
        return f'feasible={feasible}, expected={case.expected_feasible}'
    if feasible:
        residual = QcpSolver2D.evalConic(case.conic, pt)
        if residual > 1.01e-4:
            return f'quadratic violation={residual:.6g}'
        if case.poly.count >= 3 and not QcpSolver2D.pointInPolygon(case.poly, pt, 2e-5):
            return 'outside polygon'
        objective = case.a.x * pt.x + case.a.y * pt.y
        if check_expected and case.expected_objective is not None:
            if not math.isclose(objective, case.expected_objective, rel_tol=2e-4, abs_tol=2e-4):
                return f'objective={objective}, expected={case.expected_objective}'
    return None


def percentile(values, fraction):
    return sorted(values)[max(0, math.ceil(len(values) * fraction) - 1)]


def measure(case, seconds, samples, deadline_us):
    solve = case.solve
    for _ in range(100):
        solve()
    # 批量计时用于估算稳态吞吐；逐调用计时用于观察延迟尾部。
    count = 0
    start = time.perf_counter_ns()
    end = start + seconds * 1e9
    while time.perf_counter_ns() < end:
        for _ in range(100):
            solve()
        count += 100
    elapsed = (time.perf_counter_ns() - start) / 1e9
    latency = []
    for _ in range(samples):
        t0 = time.perf_counter_ns()
        solve()
        latency.append((time.perf_counter_ns() - t0) / 1000)
    return dict(name=case.name, vertices=case.poly.count, calls=count,
                elapsed_seconds=elapsed, calls_per_second=count / elapsed,
                average_us=elapsed * 1e6 / count,
                latency_median_us=statistics.median(latency),
                latency_p95_us=percentile(latency, .95),
                latency_p99_us=percentile(latency, .99), latency_max_us=max(latency),
                deadline_us=deadline_us, deadline_misses=sum(t > deadline_us for t in latency),
                latency_samples=samples)


def random_cases(count, seed):
    rng = random.Random(seed)
    result = []
    for i in range(count):
        n = rng.randint(3, kQcpMaxPoints)
        poly = polygon(n, rng.uniform(.5, 3))
        angle = rng.uniform(-math.pi, math.pi)
        u, v = math.cos(angle), math.sin(angle)
        e1, e2 = 10**rng.uniform(-6, 3), 10**rng.uniform(-6, 3)
        q00, q01, q11 = e1*u*u+e2*v*v, (e1-e2)*u*v, e1*v*v+e2*u*u
        cx, cy = rng.uniform(-2, 2), rng.uniform(-2, 2)
        # g = .5*(x-center)^T Q*(x-center) - level；Q 保证半正定。
        conic = Conic2f(q00, q01, q11, -q00*cx-q01*cy, -q01*cx-q11*cy,
                        .5*(q00*cx*cx+2*q01*cx*cy+q11*cy*cy)-rng.uniform(-.1, 2))
        result.append(Case(f'random_{i}', poly, conic, Vec2f(rng.uniform(-1, 1), rng.uniform(-1, 1))))
    return result


def stress(pool, seconds, deadline_us):
    """持续单线程饱和运行，按批次记录吞吐以观察长时间波动。"""
    start = time.perf_counter()
    total = 0
    rates = []
    worst = (0, None)
    misses = 0
    while time.perf_counter() - start < seconds:
        t0 = time.perf_counter()
        for _ in range(1000):
            case = pool[total % len(pool)]
            before = time.perf_counter_ns()
            case.solve()
            elapsed_us = (time.perf_counter_ns() - before) / 1000
            if elapsed_us > worst[0]:
                worst = (elapsed_us, case.name)
            misses += elapsed_us > deadline_us
            total += 1
        rates.append(1000 / (time.perf_counter() - t0))
    elapsed = time.perf_counter() - start
    return dict(calls=total, elapsed_seconds=elapsed, calls_per_second=total/elapsed,
                min_batch_calls_per_second=min(rates), max_batch_calls_per_second=max(rates),
                max_us=worst[0], slowest_case=worst[1], deadline_misses=misses)


def positive(value):
    x = float(value)
    if not math.isfinite(x) or x <= 0:
        raise argparse.ArgumentTypeError('必须是有限正数')
    return x


def positive_int(value):
    x = int(value)
    if x <= 0:
        raise argparse.ArgumentTypeError('必须是正整数')
    return x


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--seconds', type=positive, default=.2, help='每个场景的吞吐计时秒数，默认 .2')
    parser.add_argument('--samples', type=positive_int, default=2000, help='每场景单次延迟样本数')
    parser.add_argument('--stress-seconds', type=positive, default=3, help='持续饱和压测秒数，默认 3')
    parser.add_argument('--random-cases', type=positive_int, default=1000, help='固定种子随机输入数量')
    parser.add_argument('--seed', type=int, default=20261031)
    parser.add_argument('--deadline-us', type=positive, default=1000, help='延迟预算，默认 1000 微秒 (1kHz)')
    parser.add_argument('--disable-gc', action='store_true', help='关闭循环垃圾回收，默认保留实际运行状态')
    parser.add_argument('--output', type=Path, help='可选 JSON 报告路径')
    args = parser.parse_args()
    fixed = cases()
    for case in fixed:
        error = validate(case)
        if error:
            raise RuntimeError(f'{case.name}: {error}')
    random_pool = random_cases(args.random_cases, args.seed)
    violations = []
    for case in random_pool:
        error = validate(case, check_expected=False)
        if error:
            point, feasible = case.solve()
            violations.append(dict(name=case.name, error=error, input=asdict(case),
                                   result=asdict(point), feasible=feasible,
                                   residual=QcpSolver2D.evalConic(case.conic, point)))
    print(f'固定场景 {len(fixed)} 项正确性检查通过；随机可行性检查异常 {len(violations)}/{len(random_pool)}。', flush=True)
    print('随机检查只验证返回可行点的约束，不认证全局最优或不可行结论。', flush=True)
    gc_was_enabled = gc.isenabled()
    if args.disable_gc:
        gc.disable()
    try:
        rows = []
        print(f'{"场景":28} {"平均us":>9} {"次/秒":>10} {"P99us":>9} {"最大us":>10} {"超预算":>8}', flush=True)
        for case in fixed:
            row = measure(case, args.seconds, args.samples, args.deadline_us)
            rows.append(row)
            print(f'{case.name:28} {row["average_us"]:9.2f} {row["calls_per_second"]:10.0f} '
                  f'{row["latency_p99_us"]:9.2f} {row["latency_max_us"]:10.2f} {row["deadline_misses"]:8}', flush=True)
        print(f'持续饱和压测 {args.stress_seconds:g} 秒……', flush=True)
        sustained = stress(fixed + random_pool, args.stress_seconds, args.deadline_us)
    finally:
        if gc_was_enabled:
            gc.enable()
    print(f'持续吞吐 {sustained["calls_per_second"]:.0f} 次/秒，最大延迟 {sustained["max_us"]:.2f} us，'
          f'超预算 {sustained["deadline_misses"]}/{sustained["calls"]}。')
    print('观测最大延迟包含系统调度及计时开销，不是最坏执行时间保证。吞吐测试不包含构造输入或半平面裁剪。')
    report = dict(environment=dict(python=sys.version, platform=platform.platform(),
                                   processor=platform.processor(), gc_enabled=not args.disable_gc,
                                   clock=time.get_clock_info('perf_counter').implementation),
                  configuration={k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()},
                  scenarios=rows, stress=sustained, random_feasibility_violations=violations)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding='utf-8')
        print(f'报告已保存：{args.output.resolve()}')
    return 1 if violations else 0


if __name__ == '__main__':
    sys.exit(main())
