#!/usr/bin/env python3
"""
Check that all string type annotations have corresponding TYPE_CHECKING imports.

This script ensures that when a function/method uses string type annotations
(e.g., -> "ClassName"), the corresponding class is imported in a TYPE_CHECKING block.

Usage:
    python scripts/check_type_checking.py file1.py file2.py ...
"""

import ast
import sys
import re
from pathlib import Path
from typing import Set, List, Tuple


class TypeAnnotationVisitor(ast.NodeVisitor):
    """Extract all string type annotations from an AST."""
    
    def __init__(self, has_future_annotations: bool = False):
        self.string_annotations: Set[str] = set()
        self.has_future_annotations = has_future_annotations
        self.defined_classes: Set[str] = set()
        
    def visit_ClassDef(self, node: ast.ClassDef) -> None:
        """Track class definitions in the file."""
        self.defined_classes.add(node.name)
        self.generic_visit(node)
    def visit_FunctionDef(self, node: ast.FunctionDef) -> None:
        """Visit function definitions to extract return type annotations."""
        # Check return annotation
        if node.returns:
            self._extract_annotation(node.returns)
        
        # Check argument annotations
        for arg in node.args.args + node.args.posonlyargs + node.args.kwonlyargs:
            if arg.annotation:
                self._extract_annotation(arg.annotation)
        
        if node.args.vararg and node.args.vararg.annotation:
            self._extract_annotation(node.args.vararg.annotation)
        
        if node.args.kwarg and node.args.kwarg.annotation:
            self._extract_annotation(node.args.kwarg.annotation)
        
        self.generic_visit(node)
    
    def visit_AsyncFunctionDef(self, node: ast.AsyncFunctionDef) -> None:
        """Visit async function definitions."""
        self.visit_FunctionDef(node)  # Same logic as regular functions
    
    def visit_AnnAssign(self, node: ast.AnnAssign) -> None:
        """Visit annotated assignments (class/instance variables)."""
        self._extract_annotation(node.annotation)
        self.generic_visit(node)
    
    def _extract_annotation(self, annotation: ast.expr) -> None:
        """Extract type names from an annotation node."""
        if isinstance(annotation, ast.Constant) and isinstance(annotation.value, str):
            # String annotation like "ClassName"
            self._parse_string_annotation(annotation.value)
        elif isinstance(annotation, ast.Subscript):
            # Generic types like List["ClassName"]
            if isinstance(annotation.slice, ast.Constant) and isinstance(annotation.slice.value, str):
                self._parse_string_annotation(annotation.slice.value)
            elif isinstance(annotation.slice, ast.Tuple):
                for elt in annotation.slice.elts:
                    if isinstance(elt, ast.Constant) and isinstance(elt.value, str):
                        self._parse_string_annotation(elt.value)
    
    def _parse_string_annotation(self, annotation: str) -> None:
        """Parse a string annotation to extract class names."""
        # Remove quotes and whitespace
        annotation = annotation.strip().strip('"').strip("'")
        
        # Extract class names (simple heuristic)
        # Matches: ClassName, module.ClassName, Optional[ClassName], etc.
        # This is a simplified parser - may need enhancement for complex cases
        class_names = re.findall(r'\b([A-Z][a-zA-Z0-9_]*)\b', annotation)
        
        # If file has __future__ annotations, filter out self-referential types
        if self.has_future_annotations:
            class_names = [name for name in class_names if name not in self.defined_classes]
        
        self.string_annotations.update(class_names)


class TypeCheckingImportVisitor(ast.NodeVisitor):
    """Extract all imports from TYPE_CHECKING blocks."""
    
    def __init__(self):
        self.type_checking_imports: Set[str] = set()
        self._in_type_checking = False
    
    def visit_If(self, node: ast.If) -> None:
        """Visit if statements to find TYPE_CHECKING blocks."""
        # Check if this is "if TYPE_CHECKING:"
        if self._is_type_checking_condition(node.test):
            self._in_type_checking = True
            self.generic_visit(node)
            self._in_type_checking = False
        else:
            self.generic_visit(node)
    
    def visit_ImportFrom(self, node: ast.ImportFrom) -> None:
        """Visit from...import statements."""
        if self._in_type_checking:
            for alias in node.names:
                self.type_checking_imports.add(alias.name)
        self.generic_visit(node)
    
    def visit_Import(self, node: ast.Import) -> None:
        """Visit import statements."""
        if self._in_type_checking:
            for alias in node.names:
                # Extract the last part of the module name
                name = alias.name.split('.')[-1]
                self.type_checking_imports.add(name)
        self.generic_visit(node)
    
    def _is_type_checking_condition(self, test: ast.expr) -> bool:
        """Check if a test expression is TYPE_CHECKING."""
        if isinstance(test, ast.Name) and test.id == 'TYPE_CHECKING':
            return True
        if isinstance(test, ast.Attribute) and test.attr == 'TYPE_CHECKING':
            return True
        return False


def check_file(filepath: Path) -> Tuple[bool, List[str]]:
    """
    Check a single Python file for missing TYPE_CHECKING imports.
    
    Returns:
        (success, missing_imports): success is True if all checks pass
    """
    try:
        content = filepath.read_text(encoding='utf-8')
        tree = ast.parse(content, filename=str(filepath))
    except SyntaxError as e:
        print(f"⚠️  Syntax error in {filepath}: {e}")
        return True, []  # Don't fail on syntax errors (let other tools handle it)
    except Exception as e:
        print(f"⚠️  Error reading {filepath}: {e}")
        return True, []
    
    # Check if file uses __future__ annotations
    has_future_annotations = 'from __future__ import annotations' in content
    
    # Extract string annotations
    annotation_visitor = TypeAnnotationVisitor(has_future_annotations)
    annotation_visitor.visit(tree)
    
    # Extract TYPE_CHECKING imports
    import_visitor = TypeCheckingImportVisitor()
    import_visitor.visit(tree)
    
    # Find missing imports
    missing = annotation_visitor.string_annotations - import_visitor.type_checking_imports
    
    # Filter out built-in types and typing module types
    builtin_types = {
        'Any', 'Optional', 'Union', 'List', 'Dict', 'Set', 'Tuple', 
        'Callable', 'Type', 'TypeVar', 'Generic', 'Protocol',
        'Literal', 'Final', 'ClassVar', 'Annotated',
    }
    missing = missing - builtin_types
    
    if missing:
        return False, sorted(missing)
    
    return True, []


def main() -> int:
    """Main entry point."""
    if len(sys.argv) < 2:
        print("Usage: check_type_checking.py <file1.py> [file2.py ...]")
        return 0
    
    files = [Path(f) for f in sys.argv[1:]]
    all_passed = True
    
    for filepath in files:
        if not filepath.exists():
            print(f"⚠️  File not found: {filepath}")
            continue
        
        success, missing = check_file(filepath)
        
        if not success:
            all_passed = False
            print(f"\n❌ {filepath}")
            print(f"   Missing TYPE_CHECKING imports: {', '.join(missing)}")
            print(f"   Add to TYPE_CHECKING block:")
            print(f"   ")
            print(f"   if TYPE_CHECKING:")
            for name in missing:
                print(f"       from .module import {name}")
    
    if all_passed:
        print("✅ All TYPE_CHECKING imports are correct!")
        return 0
    else:
        print("\n" + "="*60)
        print("❌ Some files have missing TYPE_CHECKING imports.")
        print("   See: .github/instructions/python-typing.md")
        return 1


if __name__ == '__main__':
    sys.exit(main())
