# ZenOpt

最优化方法工程，按语言和算法组织实现，示例与性能测试独立存放。

```text
cpp/
  qcp/                         # QCP 的 C++ 实现
python/
  qcp/                       # QCP 的 Python 实现及公开接口
examples/
  python/qcp/                  # 使用示例
benchmarks/
  python/qcp/                  # 性能与压力测试脚本
reports/
  benchmarks/python/qcp/       # 实测结果，非算法源代码
docs/                          # 说明文档（见下方链接）
```

## 运行 QCP

在仓库根目录执行，Python 代码不依赖第三方库：

```bash
python examples/python/qcp/example_qcp.py
python benchmarks/python/qcp/benchmark_qcp.py
```

完整压测：

```bash
python benchmarks/python/qcp/benchmark_qcp.py --seconds 2 --stress-seconds 60 --samples 10000 --random-cases 10000 --output reports/benchmarks/python/qcp/benchmark_results.json
```

上述脚本也支持从其他工作目录通过绝对路径运行。`--output` 的相对路径相对于运行时的工作目录；每次使用同名输出路径会覆盖该报告，保留多次结果时请更换文件名。

应用代码导入时，将本仓库的 `python/` 加入 `PYTHONPATH`，使用：

```python
from qcp import Conic2f, Polygon2f, QcpSolver2D, Vec2f
```

性能指标及局限见 [QCP 性能测试说明](docs/benchmarks/qcp.md)。历史报告保留当时的运行参数，目录迁移不改变历史测量值。已保存报告中的 8 个随机可行性异常是原算法的数值问题；本次目录整理未修改算法。

## 新增算法约定

以 `new_algorithm` 为例：

- C++ 实现放在 `cpp/new_algorithm/`。
- Python 实现放在 `python/new_algorithm/`，通过 `__init__.py` 导出接口。
- 示例放在 `examples/<语言>/new_algorithm/`。
- 性能测试放在 `benchmarks/<语言>/new_algorithm/`。
- 性能报告放在 `reports/benchmarks/<语言>/new_algorithm/`，文件名可包含日期和设备名。
- 说明文档放在 `docs/`；以后新增自动化正确性测试时放在 `tests/<语言>/new_algorithm/`。

新增算法沿用目录即可，不需要为尚未实现的算法提前创建空目录。
