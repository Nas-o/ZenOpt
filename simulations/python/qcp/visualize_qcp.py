"""交互 QCP 仿真；无桌面环境时可用 --headless 导出扫描图和 CSV。"""
import argparse
import csv
from pathlib import Path
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / 'python'))

import numpy as np
from qcp import Conic2f, Polygon2f, QcpSolver2D, Vec2f


class Simulation:
    def __init__(self, grid_size=161, repeats=30):
        self.grid_size = grid_size
        self.repeats = repeats
        self.history = []

    def evaluate(self, angle=45., radius=1., clip=.6, vertices=4, center=0., mode='Circle'):
        # 正方形使用水平/垂直边；所有输入构造均在计时区间外。
        n = int(vertices)
        theta = np.arange(n) * 2 * np.pi / n + (np.pi / 4 if n == 4 else 0)
        scale = 2 * np.sqrt(2) if n == 4 else 2.
        poly = Polygon2f([Vec2f(float(scale*np.cos(t)), float(scale*np.sin(t))) for t in theta])
        QcpSolver2D.clipHalfPlane(poly, 1., 0., clip)
        qy = 1e-8 if mode == 'Near-singular' else 2.
        conic = Conic2f(q00=2., q11=qy, p0=-2*center,
                        r=center*center + (radius*radius if mode == 'Infeasible' else -radius*radius))
        a = Vec2f(float(np.cos(np.deg2rad(angle))), float(np.sin(np.deg2rad(angle))))
        # 当前参数下预热；计时仅包围原 Python solve，绘图和网格不计入。
        for _ in range(5):
            QcpSolver2D.solve(poly, conic, a)
        latencies = []
        for _ in range(self.repeats):
            start = time.perf_counter_ns()
            point, feasible = QcpSolver2D.solve(poly, conic, a)
            latencies.append((time.perf_counter_ns()-start)/1000)
        residual = QcpSolver2D.evalConic(conic, point)
        inside = QcpSolver2D.pointInPolygon(poly, point)
        valid = feasible and inside and residual <= 1e-4
        objective = a.x*point.x + a.y*point.y
        coords = np.linspace(-2.25, 2.25, self.grid_size)
        x, y = np.meshgrid(coords, coords)
        start = time.perf_counter_ns()
        mask = np.ones(x.shape, dtype=bool)
        if poly.count < 3:
            mask[:] = False
        else:
            pts = poly.pts
            area = sum(p.x*pts[(i+1)%len(pts)].y-pts[(i+1)%len(pts)].x*p.y for i,p in enumerate(pts))
            for i, p in enumerate(pts):
                p2 = pts[(i+1)%len(pts)]
                cross = (p2.x-p.x)*(y-p.y)-(p2.y-p.y)*(x-p.x)
                mask &= cross >= -1e-10 if area > 0 else cross <= 1e-10
        g = .5*(conic.q00*x*x + 2*conic.q01*x*y + conic.q11*y*y)+conic.p0*x+conic.p1*y+conic.r
        mask &= g <= 0
        scores = np.where(mask, a.x*x+a.y*y, -np.inf)
        index = np.unravel_index(np.argmax(scores), scores.shape)
        grid_objective = float(scores[index]) if mask.any() else float('nan')
        grid_point = (float(x[index]), float(y[index])) if mask.any() else None
        grid_us = (time.perf_counter_ns()-start)/1000
        row = dict(frame=len(self.history), angle_deg=float(angle), vertices=poly.count,
                   mode=mode, radius=radius, clip=clip, center=center,
                   x=point.x, y=point.y, feasible=feasible, validated_feasible=valid,
                   residual=residual, objective=objective if valid else float('nan'),
                   grid_objective=grid_objective, median_us=float(np.median(latencies)),
                   max_us=max(latencies), grid_search_us=grid_us)
        self.history.append(row)
        return dict(row=row, poly=poly, conic=conic, a=a, point=point, x=x, y=y,
                    g=g, mask=mask, grid_point=grid_point)


def make_figure():
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(13, 8))
    layout = fig.add_gridspec(2, 2, left=.07, right=.96, bottom=.34, top=.91,
                              width_ratios=[1.15, 1], hspace=.5, wspace=.3)
    axes = [fig.add_subplot(layout[:, 0]), fig.add_subplot(layout[0, 1]), fig.add_subplot(layout[1, 1])]
    fig.suptitle('QCP | feasible geometry, solution quality and solve latency', fontsize=15)
    return fig, axes


def draw(fig, axes, sim, data):
    geo, quality, timing = axes
    for ax in axes:
        ax.clear()
    row, poly, a, pt = data['row'], data['poly'], data['a'], data['point']
    x, y, g, mask = data['x'], data['y'], data['g'], data['mask']
    if mask.any():
        geo.contourf(x, y, mask.astype(float), levels=[.5, 1.5], colors=['#b9dfcf'], alpha=.8)
    if g.min() <= 0 <= g.max():
        geo.contour(x, y, g, levels=[0], colors=['#298e6c'], linewidths=1.5)
    if poly.count:
        pts = poly.pts + poly.pts[:1]
        geo.plot([p.x for p in pts], [p.y for p in pts], color='#455d7a', label='Linear boundary')
    if data['grid_point'] is not None:
        geo.scatter(*data['grid_point'], marker='x', s=75, color='#d08a1b', label='Grid best')
    if row['validated_feasible']:
        geo.scatter(pt.x, pt.y, marker='*', s=180, color='#c13e59', label='QCP solution', zorder=5)
        # 目标等值线，法向量即最大化方向。
        direction = np.array([-a.y, a.x])
        endpoints = np.array([pt.x, pt.y]) + np.array([-4, 4])[:, None]*direction
        geo.plot(endpoints[:, 0], endpoints[:, 1], '--', color='#c13e59', alpha=.6)
    else:
        geo.scatter(pt.x, pt.y, marker='s', s=55, color='#c13e59', label='Fallback (not certified)')
    geo.arrow(-1.7, -1.7, .65*a.x, .65*a.y, width=.012, head_width=.12,
              length_includes_head=True, color='#6554a0')
    geo.set(xlim=(-2.25,2.25), ylim=(-2.25,2.25), xlabel='x', ylabel='y', title='Green = sampled feasible region; arrow = objective')
    geo.set_aspect('equal')
    geo.legend(loc='upper right', fontsize=8)
    history = sim.history[-180:]
    frames = [r['frame'] for r in history]
    quality.plot(frames, [r['objective'] for r in history], color='#c13e59', label='QCP (validated)')
    quality.plot(frames, [r['grid_objective'] for r in history], '--', color='#d08a1b', label='Grid lower bound')
    quality.set(xlabel='Parameter update / frame', ylabel='Objective value', title='Solution quality (larger is better)')
    quality.legend(fontsize=8)
    timing.plot(frames, [r['median_us'] for r in history], color='#455d7a', label='QCP median')
    timing.plot(frames, [r['max_us'] for r in history], color='#c13e59', alpha=.6, label='QCP max sample')
    timing.axhline(1000, color='#777777', linestyle=':', label='1 ms budget')
    timing.set(yscale='log', xlabel='Parameter update / frame', ylabel='Solve latency (us)', title='Measured solve time (drawing excluded)')
    timing.legend(fontsize=8, ncol=3)
    for ax in axes:
        ax.grid(alpha=.17)
    status = 'VALID' if row['validated_feasible'] else ('INFEASIBLE' if not row['feasible'] else 'CONSTRAINT VIOLATION')
    gap = row['objective']-row['grid_objective']
    fig._qcp_status.set_text(
        f'{status} | angle={row["angle_deg"]:.0f} deg | g(x)={row["residual"]:.2e} | QCP-grid={gap:.3g}\n'
        f'Solve median={row["median_us"]:.1f} us, max={row["max_us"]:.1f} us ({sim.repeats} samples) | '
        f'Grid search={row["grid_search_us"]:.0f} us ({sim.grid_size} x {sim.grid_size})')
    fig.canvas.draw_idle()


def save_results(sim, path):
    with path.open('w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=list(sim.history[0]))
        writer.writeheader()
        writer.writerows(sim.history)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--headless', action='store_true', help='不打开窗口，扫描目标方向并导出 PNG/CSV')
    parser.add_argument('--frames', type=int, default=72)
    parser.add_argument('--grid-size', type=int, default=161)
    parser.add_argument('--repeats', type=int, default=30)
    parser.add_argument('--output-dir', type=Path,
                        default=Path(__file__).resolve().parents[3]/'reports/simulations/python/qcp')
    args = parser.parse_args()
    if args.frames < 2 or args.grid_size < 11 or args.repeats < 1:
        parser.error('frames >= 2, grid-size >= 11, repeats >= 1')
    import matplotlib
    if args.headless:
        matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.widgets import Button, RadioButtons, Slider
    sim = Simulation(args.grid_size, args.repeats)
    fig, axes = make_figure()
    fig._qcp_status = fig.text(.07, .24, '', fontsize=10, va='top', family='monospace')
    args.output_dir.mkdir(parents=True, exist_ok=True)
    if args.headless:
        for angle in np.linspace(0, 360, args.frames, endpoint=False):
            data = sim.evaluate(angle=float(angle))
        draw(fig, axes, sim, data)
        fig.text(.07, .14, 'Direction sweep with fixed constraints. This is an analytic solver, not an iterative convergence trace.\n'
                 'Grid is an approximate lower bound; a small positive objective gap is expected.', fontsize=10)
        fig.savefig(args.output_dir/'direction_sweep.png', dpi=150)
        save_results(sim, args.output_dir/'direction_sweep.csv')
        plt.close(fig)
        print(f'Saved {args.frames} frames: {args.output_dir.resolve()}')
        return
    sliders = {}
    specifications = [('Angle', 0, 360, 45, None), ('Radius', .1, 2., 1., None),
                      ('Clip x', -.9, 2., .6, None), ('Vertices', 3, 16, 4, 1),
                      ('Center x', -1., 1., 0., None)]
    for i, (label, lo, hi, value, step) in enumerate(specifications):
        ax = fig.add_axes([.12, .175-i*.033, .45, .02])
        sliders[label] = Slider(ax, label, lo, hi, valinit=value, valstep=step)
    modes = RadioButtons(fig.add_axes([.65,.08,.17,.11]), ['Circle', 'Near-singular', 'Infeasible'])
    state = {'playing': False}

    def update(_=None):
        values = [s.val for s in sliders.values()]
        data = sim.evaluate(*values, mode=modes.value_selected)
        draw(fig, axes, sim, data)

    def toggle(_):
        state['playing'] = not state['playing']
        play.label.set_text('Pause' if state['playing'] else 'Play sweep')
        fig.canvas.draw_idle()

    def tick():
        if state['playing']:
            sliders['Angle'].set_val((sliders['Angle'].val+4) % 360)

    def export(_):
        fig.savefig(args.output_dir/'interactive_snapshot.png', dpi=150)
        save_results(sim, args.output_dir/'interactive_history.csv')
        print(f'Saved snapshot and history: {args.output_dir.resolve()}')

    play = Button(fig.add_axes([.65,.025,.12,.04]), 'Play sweep')
    save = Button(fig.add_axes([.80,.025,.12,.04]), 'Save PNG/CSV')
    play.on_clicked(toggle)
    save.on_clicked(export)
    for slider in sliders.values():
        slider.on_changed(update)
    modes.on_clicked(update)
    timer = fig.canvas.new_timer(interval=120)
    timer.add_callback(tick)
    timer.start()
    update()
    plt.show()


if __name__ == '__main__':
    main()
