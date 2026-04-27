# Python Type Annotations Best Practices

## TYPE_CHECKING 导入规则 (强制执行)

### 规则 1: 所有用于类型注解的类必须在 TYPE_CHECKING 块中导入

当函数/方法的返回类型或参数类型使用字符串注解时，**必须**在 `TYPE_CHECKING` 块中导入对应的类型。

**❌ 错误示例：**
```python
from typing import TYPE_CHECKING

# 缺少导入！
def affine_weyl_group(self) -> "ExtendedAffineWeylGroup":
    from .weyl_group import ExtendedAffineWeylGroup
    return ExtendedAffineWeylGroup(self)
```

**✅ 正确示例：**
```python
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .weyl_group import ExtendedAffineWeylGroup  # 必须导入！

def affine_weyl_group(self) -> "ExtendedAffineWeylGroup":
    from .weyl_group import ExtendedAffineWeylGroup
    return ExtendedAffineWeylGroup(self)
```

### 规则 2: 检查所有字符串类型注解

在编写代码时，必须检查：
1. 所有返回类型注解中的字符串类型
2. 所有参数类型注解中的字符串类型
3. 所有类属性类型注解中的字符串类型

这些类型都必须在 `TYPE_CHECKING` 块中导入。

### 规则 3: 避免循环导入

使用 `TYPE_CHECKING` 的主要目的是避免循环导入：

```python
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    # 这些导入只在类型检查时执行，不会在运行时导入
    from .module_a import ClassA
    from .module_b import ClassB
```

### 规则 4: 完整的 TYPE_CHECKING 模板

```python
from typing import TYPE_CHECKING, Optional, List, Dict, Any

if TYPE_CHECKING:
    # 导入所有用于类型注解的类
    from .module1 import Class1, Class2
    from .module2 import Class3
    # 可以使用 TYPE_CHECKING 避免循环导入
```

## 为什么这很重要？

### 对 IDE 和 Language Server 的影响

1. **VSCode Pylance**: 严格要求 `TYPE_CHECKING` 导入，否则无法提供自动补全
2. **PyCharm**: 相对宽松，但最佳实践仍然要求显式导入
3. **Mypy/Pyright**: 类型检查器需要这些导入来验证类型正确性

### 实际影响

**没有 TYPE_CHECKING 导入时：**
- ❌ IDE 自动补全失败
- ❌ 类型检查器报错
- ❌ 代码导航（Go to Definition）失败
- ❌ 重构工具无法正确识别类型

**有 TYPE_CHECKING 导入时：**
- ✅ 完整的 IDE 支持
- ✅ 类型检查通过
- ✅ 代码导航正常工作
- ✅ 重构工具正常工作

## 检查清单

在提交代码前，检查：

- [ ] 所有字符串类型注解对应的类都在 `TYPE_CHECKING` 块中导入
- [ ] 没有未使用的 `TYPE_CHECKING` 导入
- [ ] 运行时导入和类型检查导入分离清晰
- [ ] 使用 `mypy` 或 `pyright` 验证类型注解正确性

## 自动检查

可以使用以下工具自动检查：

```bash
# 使用 mypy 检查类型
mypy --strict your_module.py

# 使用 pyright 检查类型
pyright your_module.py
```

## 参考资料

- [PEP 484 - Type Hints](https://www.python.org/dev/peps/pep-0484/)
- [PEP 563 - Postponed Evaluation of Annotations](https://www.python.org/dev/peps/pep-0563/)
- [typing.TYPE_CHECKING documentation](https://docs.python.org/3/library/typing.html#typing.TYPE_CHECKING)
