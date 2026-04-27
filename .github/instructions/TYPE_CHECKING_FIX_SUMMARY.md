# TYPE_CHECKING 修复总结

## 🎯 问题描述

在 Jupyter Notebook 中，输入 `What.` 后没有触发自动补全下拉列表，原因是 `ExtendedAffineWeylGroup` 类型没有在 `TYPE_CHECKING` 块中导入，导致 VSCode Pylance 无法识别变量类型。

## ✅ 已完成的修复

### 1. 修复代码文件

**`pyw/core/affine_lie_algebra.py`**
```python
if TYPE_CHECKING:
    from .affine_weight import AffineWeight
    from .weyl_group import FiniteWeylGroup, ExtendedAffineWeylGroup
```

添加了缺失的类型导入，确保所有字符串类型注解都有对应的 TYPE_CHECKING 导入。

### 2. 创建规则文档

**`.github/instructions/python-typing.md`**
- 详细说明 TYPE_CHECKING 的使用规则
- 提供错误/正确示例对比
- 解释为什么这对 IDE 很重要

### 3. 创建自动检查工具

**`scripts/check_type_checking.py`**
智能检查脚本，特性包括：
- ✅ 自动识别 `from __future__ import annotations`
- ✅ 过滤自引用类型（类返回自身类型）
- ✅ 过滤内置类型（`Optional`, `List`, `Dict` 等）
- ✅ 清晰的错误提示和修复建议

### 4. 配置 Pre-commit Hooks

**`.pre-commit-config.yaml`**
包含以下检查：
- Black (代码格式化)
- Ruff (代码检查)
- Mypy (类型检查)
- **自定义 TYPE_CHECKING 检查** ⭐
- 通用文件检查（空白字符、文件结尾等）

## 📊 检查结果

```bash
✅ All TYPE_CHECKING imports are correct!
```

所有 `pyw` 包中的 Python 文件都通过了检查。

## 🔧 使用方法

### 手动运行检查

```bash
# 检查单个文件
python scripts/check_type_checking.py pyw/core/affine_lie_algebra.py

# 检查所有文件
python scripts/check_type_checking.py pyw/**/*.py
```

### Pre-commit 自动检查

```bash
# 已安装并配置
pre-commit install  # ✅ 已完成

# 手动运行所有检查
pre-commit run --all-files

# 只运行 TYPE_CHECKING 检查
pre-commit run check-type-checking-imports --all-files
```

### Git Commit 时自动检查

现在每次 `git commit` 时，pre-commit 会自动运行所有检查，包括 TYPE_CHECKING 导入检查。

## 🎓 学到的知识

### 为什么 Kiro 能工作而 VSCode 不行？

| 特性 | VSCode Pylance | Kiro/其他编辑器 |
|------|----------------|----------------|
| 类型检查严格性 | 高（严格遵守 PEP 484） | 中/低（更宽松） |
| 需要显式 TYPE_CHECKING 导入 | 是 | 否 |
| 类型推断策略 | 保守（只看声明） | 激进（扫描整个项目） |

**结论：** VSCode 的严格性虽然要求更多，但能帮助写出更规范、更易维护的代码。

### `from __future__ import annotations` 的作用

当使用这个导入时：
- 所有类型注解自动变成字符串（延迟求值）
- 自引用类型不需要在 TYPE_CHECKING 中导入
- 例如：`AffineWeight` 类的方法返回 `AffineWeight` 不需要额外导入

## 📝 给 AI 助手的规则

在 `.github/instructions/python-typing.md` 中定义了强制规则：

**规则 1：** 所有字符串类型注解必须在 TYPE_CHECKING 块中导入

```python
# ❌ 错误
def method(self) -> "ClassName":
    ...

# ✅ 正确
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .module import ClassName

def method(self) -> "ClassName":
    ...
```

## 🚀 后续建议

1. **定期运行检查**
   ```bash
   pre-commit run --all-files
   ```

2. **更新 pre-commit hooks**
   ```bash
   pre-commit autoupdate
   ```

3. **在 CI/CD 中集成**
   ```yaml
   # .github/workflows/ci.yml
   - name: Run pre-commit
     run: pre-commit run --all-files
   ```

## 📚 参考资料

- [PEP 484 - Type Hints](https://www.python.org/dev/peps/pep-0484/)
- [PEP 563 - Postponed Evaluation of Annotations](https://www.python.org/dev/peps/pep-0563/)
- [typing.TYPE_CHECKING documentation](https://docs.python.org/3/library/typing.html#typing.TYPE_CHECKING)

---

**修复日期：** 2026-02-01  
**状态：** ✅ 完成并测试通过
