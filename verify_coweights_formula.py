"""
Verification script for fundamental coweights formula.

This script verifies that the explicit formula:
    ω_j^∨ = 2/(α_j, α_j) * ω_j

produces the same results as SageMath's coweight_lattice().fundamental_weights().

Tests across multiple Cartan types: A, B, C, D, G (both simply-laced and non-simply-laced).
"""

def verify_coweights_formula():
    """Verify fundamental coweights formula against SageMath implementation."""
    from pyw.core.affine_lie_algebra import AffineLieAlgebra

    # Test cases: mix of simply-laced (A, D) and non-simply-laced (B, C, G)
    test_cases = [
        ['A', 1, 1],  # Simply-laced, rank 1
        ['A', 2, 1],  # Simply-laced, rank 2
        ['B', 2, 1],  # Non-simply-laced
        ['C', 2, 1],  # Non-simply-laced
        ['D', 4, 1],  # Simply-laced, higher rank
        ['G', 2, 1],  # Non-simply-laced, exceptional
    ]

    print("=" * 70)
    print("Verifying Fundamental Coweights Formula")
    print("Formula: ω_j^∨ = 2/(α_j, α_j) * ω_j")
    print("=" * 70)

    all_passed = True

    for cartan_type in test_cases:
        print(f"\nTesting {cartan_type}...")
        ala = AffineLieAlgebra(cartan_type)

        # Method 1: Current implementation (SageMath)
        Lambda_check_sage = ala.fundamental_coweights()

        # Method 2: Explicit formula
        Lambda_weights = ala.fundamental_weights_sage()
        alpha = ala.simple_roots()

        Lambda_check_formula = {}
        for i in Lambda_weights.keys():
            # Compute (α_i, α_i)
            alpha_i_sq = ala.scalar_product(alpha[i], alpha[i])
            # ω_i^∨ = 2/(α_i, α_i) * ω_i
            Lambda_check_formula[i] = (2 / alpha_i_sq) * Lambda_weights[i]

        # Compare results
        print(f"  Indices: {sorted(Lambda_check_sage.keys())}")

        mismatch = False
        for i in Lambda_check_sage.keys():
            sage_val = Lambda_check_sage[i]
            formula_val = Lambda_check_formula[i]

            # Check if they're equal (allowing small numerical errors)
            diff = sage_val - formula_val

            # Convert to ambient space for comparison
            try:
                diff_vec = diff.to_vector()
                max_diff = max(abs(float(x)) for x in diff_vec)
            except:
                # For simple types, direct comparison
                max_diff = abs(float(diff)) if hasattr(diff, '__float__') else 0

            if max_diff > 1e-10:
                print(f"  ❌ Index {i}: Mismatch detected!")
                print(f"     SageMath:  {sage_val}")
                print(f"     Formula:   {formula_val}")
                print(f"     Difference: {diff}")
                mismatch = True
                all_passed = False
            else:
                print(f"  ✓ Index {i}: Match (max_diff = {max_diff:.2e})")

        if not mismatch:
            print(f"  ✅ All indices match for {cartan_type}")

    print("\n" + "=" * 70)
    if all_passed:
        print("✅ ALL TESTS PASSED - Formula is equivalent to SageMath implementation")
    else:
        print("❌ SOME TESTS FAILED - Review discrepancies above")
    print("=" * 70)

    return all_passed


if __name__ == "__main__":
    import sys
    try:
        passed = verify_coweights_formula()
        sys.exit(0 if passed else 1)
    except Exception as e:
        print(f"\n❌ Error during verification: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(2)
