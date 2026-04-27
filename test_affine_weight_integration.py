#!/usr/bin/env sage -python
"""
Quick validation script for AffineWeight integration with AffineLieAlgebra.

This script tests the new functionality without requiring pytest.
Run with: sage test_affine_weight_integration.py
"""

from sage.all import RootSystem
from pyw.core.affine_lie_algebra import AffineLieAlgebra
from pyw.core.affine_weight import AffineWeight


def test_simple_reflection_with_affine_weight():
    """Test simple_reflection() accepts AffineWeight and returns AffineWeight."""
    print("=" * 60)
    print("Test 1: simple_reflection() with AffineWeight input")
    print("=" * 60)
    
    ala = AffineLieAlgebra(["A", 2, 1])
    
    # Get FINITE weights (not AffineWeight objects)
    finite_rs = RootSystem(["A", 2])
    Lambda_finite = finite_rs.weight_space().fundamental_weights()
    
    # Create AffineWeight from finite weight
    lambda_hat = AffineWeight(ala, Lambda_finite[1], level=1, grade=0)
    print(f"Input: {lambda_hat}")
    
    # Apply simple reflection with AffineWeight input
    result = ala.simple_reflection(1, lambda_hat)
    
    # Check result type
    assert isinstance(result, AffineWeight), f"Expected AffineWeight, got {type(result)}"
    print(f"✓ Result is AffineWeight: {result}")
    
    # Check finite part is reflected (compare as strings)
    expected_str = str(-Lambda_finite[1] + Lambda_finite[2])
    result_str = str(result.finite_part)
    assert result_str == expected_str, f"Expected {expected_str}, got {result_str}"
    print(f"✓ Finite part correctly reflected: {result.finite_part}")
    
    # Level and grade should be unchanged
    assert result.level == 1, f"Expected level=1, got {result.level}"
    assert result.grade == 0, f"Expected grade=0, got {result.grade}"
    print(f"✓ Level and grade preserved: level={result.level}, grade={result.grade}")
    
    print("✅ Test 1 PASSED\n")


def test_weyl_reflection_with_affine_weight():
    """Test weyl_reflection() accepts AffineWeight and returns AffineWeight."""
    print("=" * 60)
    print("Test 2: weyl_reflection() with AffineWeight input")
    print("=" * 60)
    
    ala = AffineLieAlgebra(["A", 2, 1])
    
    # Get finite weights
    finite_rs = RootSystem(["A", 2])
    Lambda_finite = finite_rs.weight_space().fundamental_weights()
    
    # Get AFFINE roots (not finite roots!)
    alpha_affine = ala.simple_roots()
    
    # Create AffineWeight
    lambda_hat = AffineWeight(ala, Lambda_finite[1], level=2, grade=1)
    print(f"Input: {lambda_hat}")
    
    # Apply Weyl reflection with AffineWeight input
    result = ala.weyl_reflection(alpha_affine[1], lambda_hat)
    
    # Check result type
    assert isinstance(result, AffineWeight), f"Expected AffineWeight, got {type(result)}"
    print(f"✓ Result is AffineWeight: {result}")
    
    # Check finite part is reflected (compare as strings)
    expected_str = str(-Lambda_finite[1] + Lambda_finite[2])
    result_str = str(result.finite_part)
    assert result_str == expected_str, f"Expected {expected_str}, got {result_str}"
    print(f"✓ Finite part correctly reflected: {result.finite_part}")
    
    # Level and grade should be preserved
    assert result.level == 2, f"Expected level=2, got {result.level}"
    assert result.grade == 1, f"Expected grade=1, got {result.grade}"
    print(f"✓ Level and grade preserved: level={result.level}, grade={result.grade}")
    
    print("✅ Test 2 PASSED\n")


def test_backward_compatibility():
    """Test that old tuple-style calls still work."""
    print("=" * 60)
    print("Test 3: Backward compatibility (tuple style)")
    print("=" * 60)
    
    ala = AffineLieAlgebra(["A", 2, 1])
    
    # Get finite weights
    finite_rs = RootSystem(["A", 2])
    Lambda_finite = finite_rs.weight_space().fundamental_weights()
    
    # Old style: pass (weight, k, n) separately
    result = ala.simple_reflection(1, Lambda_finite[1], k=1, n=0)
    
    # Should return tuple
    assert isinstance(result, tuple), f"Expected tuple, got {type(result)}"
    assert len(result) == 3, f"Expected 3-tuple, got {len(result)}-tuple"
    print(f"✓ Result is tuple: {result}")
    
    finite_part, k, n = result
    expected_str = str(-Lambda_finite[1] + Lambda_finite[2])
    result_str = str(finite_part)
    assert result_str == expected_str, f"Expected {expected_str}, got {result_str}"
    assert k == 1, f"Expected k=1, got {k}"
    assert n == 0, f"Expected n=0, got {n}"
    print(f"✓ Tuple contents correct: ({finite_part}, {k}, {n})")
    
    print("✅ Test 3 PASSED\n")


def test_finite_reflection_unchanged():
    """Test that finite reflections (no k, n) still work."""
    print("=" * 60)
    print("Test 4: Finite reflection (no k, n)")
    print("=" * 60)
    
    ala = AffineLieAlgebra(["A", 2, 1])
    
    # Get finite weights
    finite_rs = RootSystem(["A", 2])
    Lambda_finite = finite_rs.weight_space().fundamental_weights()
    
    # Finite reflection (no k, n provided)
    result = ala.simple_reflection(1, Lambda_finite[1])
    
    # Should return finite weight (not AffineWeight, not tuple)
    assert not isinstance(result, AffineWeight), f"Expected finite weight, got AffineWeight"
    assert not isinstance(result, tuple), f"Expected finite weight, got tuple"
    print(f"✓ Result is finite weight: {result}")
    
    expected_str = str(-Lambda_finite[1] + Lambda_finite[2])
    result_str = str(result)
    assert result_str == expected_str, f"Expected {expected_str}, got {result_str}"
    print(f"✓ Reflection correct: {result}")
    
    print("✅ Test 4 PASSED\n")


def main():
    """Run all tests."""
    print("\n" + "=" * 60)
    print("AffineWeight Integration Tests")
    print("=" * 60 + "\n")
    
    try:
        test_simple_reflection_with_affine_weight()
        test_weyl_reflection_with_affine_weight()
        test_backward_compatibility()
        test_finite_reflection_unchanged()
        
        print("=" * 60)
        print("🎉 ALL TESTS PASSED!")
        print("=" * 60)
        return 0
        
    except AssertionError as e:
        print(f"\n❌ TEST FAILED: {e}")
        return 1
    except Exception as e:
        print(f"\n❌ UNEXPECTED ERROR: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    exit(main())
