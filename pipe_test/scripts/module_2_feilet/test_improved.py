#!/usr/bin/env python3
# scripts/module_2/test_improved.py
"""
Quick verification script for improved Module 2 files.

This script checks that all imports work correctly and functions are accessible.
"""

import sys
from pathlib import Path


def test_imports():
    """Test that all modules can be imported."""
    print("Testing imports...")
    
    try:
        from scripts.module_2 import config_improved
        print("  ✓ config_improved")
    except ImportError as e:
        print(f"  ✗ config_improved: {e}")
        return False
    
    try:
        from scripts.module_2 import uniprot_fetch_improved
        print("  ✓ uniprot_fetch_improved")
    except ImportError as e:
        print(f"  ✗ uniprot_fetch_improved: {e}")
        return False
    
    try:
        from scripts.module_2 import cazy_fetch_improved
        print("  ✓ cazy_fetch_improved")
    except ImportError as e:
        print(f"  ✗ cazy_fetch_improved: {e}")
        return False
    
    try:
        from scripts.module_2 import ncbi_fetch_improved
        print("  ✓ ncbi_fetch_improved")
    except ImportError as e:
        print(f"  ✗ ncbi_fetch_improved: {e}")
        return False
    
    try:
        from scripts.module_2 import merge_sequences_improved
        print("  ✓ merge_sequences_improved")
    except ImportError as e:
        print(f"  ✗ merge_sequences_improved: {e}")
        return False
    
    try:
        from scripts.module_2 import run_module2_improved
        print("  ✓ run_module2_improved")
    except ImportError as e:
        print(f"  ✗ run_module2_improved: {e}")
        return False
    
    return True


def test_config():
    """Test config module functions."""
    print("\nTesting config module...")
    
    try:
        from scripts.module_2.config_improved import (
            check_dependencies,
            validate_config,
            PROJECT_ROOT,
            DATA_RAW_DIR,
            CAZY_FAMILIES
        )
        
        print(f"  Project root: {PROJECT_ROOT}")
        print(f"  Data directory: {DATA_RAW_DIR}")
        print(f"  CAZy families: {CAZY_FAMILIES}")
        
        # Check dependencies
        deps_ok = check_dependencies()
        if deps_ok:
            print("  ✓ All dependencies installed")
        else:
            print("  ✗ Missing dependencies (see above)")
            return False
        
        # Validate config
        config_ok = validate_config()
        if config_ok:
            print("  ✓ Configuration is valid")
        else:
            print("  ✗ Configuration is invalid")
            return False
        
        return True
        
    except Exception as e:
        print(f"  ✗ Error testing config: {e}")
        return False


def test_function_signatures():
    """Test that all main functions exist and have correct signatures."""
    print("\nTesting function signatures...")
    
    try:
        from scripts.module_2.uniprot_fetch_improved import fetch_uniprot_lpmos
        from scripts.module_2.cazy_fetch_improved import fetch_all_cazy_genbank_ids
        from scripts.module_2.ncbi_fetch_improved import fetch_cazy_lpmos_to_fasta
        from scripts.module_2.merge_sequences_improved import merge_uniprot_and_cazy
        from scripts.module_2.run_module2_improved import run_module2_pipeline
        
        # Check that functions are callable
        assert callable(fetch_uniprot_lpmos), "fetch_uniprot_lpmos is not callable"
        assert callable(fetch_all_cazy_genbank_ids), "fetch_all_cazy_genbank_ids is not callable"
        assert callable(fetch_cazy_lpmos_to_fasta), "fetch_cazy_lpmos_to_fasta is not callable"
        assert callable(merge_uniprot_and_cazy), "merge_uniprot_and_cazy is not callable"
        assert callable(run_module2_pipeline), "run_module2_pipeline is not callable"
        
        print("  ✓ All main functions are callable")
        return True
        
    except Exception as e:
        print(f"  ✗ Error testing functions: {e}")
        return False


def test_paths():
    """Test that paths are correctly configured."""
    print("\nTesting path configuration...")
    
    try:
        from scripts.module_2.config_improved import (
            PROJECT_ROOT,
            DATA_RAW_DIR,
            UNIPROT_FASTA,
            CAZY_FASTA,
            MERGED_FASTA
        )
        
        print(f"  Project root exists: {PROJECT_ROOT.exists()}")
        print(f"  Data dir path: {DATA_RAW_DIR}")
        print(f"  UniProt FASTA: {UNIPROT_FASTA}")
        print(f"  CAZy FASTA: {CAZY_FASTA}")
        print(f"  Merged FASTA: {MERGED_FASTA}")
        
        # Check that all paths are Path objects
        assert isinstance(PROJECT_ROOT, Path), "PROJECT_ROOT is not a Path object"
        assert isinstance(DATA_RAW_DIR, Path), "DATA_RAW_DIR is not a Path object"
        assert isinstance(UNIPROT_FASTA, Path), "UNIPROT_FASTA is not a Path object"
        assert isinstance(CAZY_FASTA, Path), "CAZY_FASTA is not a Path object"
        assert isinstance(MERGED_FASTA, Path), "MERGED_FASTA is not a Path object"
        
        print("  ✓ All paths are Path objects")
        
        # Check that project root exists
        if not PROJECT_ROOT.exists():
            print(f"  ✗ Project root does not exist: {PROJECT_ROOT}")
            return False
        
        print("  ✓ Project root exists")
        return True
        
    except Exception as e:
        print(f"  ✗ Error testing paths: {e}")
        return False


def main():
    """Run all tests."""
    print("=" * 70)
    print("MODULE 2 IMPROVED - VERIFICATION TESTS")
    print("=" * 70)
    print()
    
    results = []
    
    # Run tests
    results.append(("Imports", test_imports()))
    results.append(("Config", test_config()))
    results.append(("Function signatures", test_function_signatures()))
    results.append(("Paths", test_paths()))
    
    # Summary
    print("\n" + "=" * 70)
    print("TEST SUMMARY")
    print("=" * 70)
    
    all_passed = True
    for test_name, passed in results:
        status = "✓ PASS" if passed else "✗ FAIL"
        print(f"{status}: {test_name}")
        if not passed:
            all_passed = False
    
    print()
    
    if all_passed:
        print("✓ ALL TESTS PASSED")
        print("\nThe improved module is ready to use!")
        print("\nTo run the pipeline:")
        print("  export NCBI_EMAIL='your.email@institution.no'")
        print("  python -m scripts.module_2.run_module2_improved")
        return 0
    else:
        print("✗ SOME TESTS FAILED")
        print("\nPlease check the errors above and fix any issues.")
        return 1


if __name__ == "__main__":
    sys.exit(main())
