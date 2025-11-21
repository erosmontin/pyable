#!/usr/bin/env python
"""
PHASE 4: Regression Test Report
Tests backward compatibility and validates all fixes
"""

import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import subprocess
import json
from datetime import datetime

def run_test_suite():
    """Run comprehensive test suite and generate report"""
    
    print("=" * 80)
    print("PHASE 4: REGRESSION TESTS & COMPATIBILITY VERIFICATION")
    print("=" * 80)
    print(f"Timestamp: {datetime.now()}")
    print()
    
    results = {
        "timestamp": str(datetime.now()),
        "phase1_fixes": [],
        "phase2_unit_tests": [],
        "phase3_integration_tests": [],
        "backward_compatibility": [],
        "summary": {}
    }
    
    # Test PHASE 1 fixes
    print("\n" + "=" * 80)
    print("PHASE 1: VALIDATING STATIC ANALYSIS FIXES")
    print("=" * 80)
    
    phase1_fixes = [
        ("Fix duplicate vtk2sitk() in meshable.py", validate_no_duplicates),
        ("Fix f-string formatting errors", validate_f_strings),
        ("Fix operator precedence (& vs and)", validate_operators),
        ("Fix spelling/typos", validate_spelling),
        ("Fix variable scope in mergeLabels()", validate_mergelabels_scope),
        ("Fix NaN comparison", validate_nan_comparison),
        ("Fix concatenated method name", validate_method_names),
    ]
    
    for fix_name, validator in phase1_fixes:
        try:
            result = validator()
            status = "✓ PASS" if result else "✗ FAIL"
            print(f"{status}: {fix_name}")
            results["phase1_fixes"].append({"name": fix_name, "status": "pass" if result else "fail"})
        except Exception as e:
            print(f"✗ ERROR: {fix_name} - {str(e)}")
            results["phase1_fixes"].append({"name": fix_name, "status": "error", "error": str(e)})
    
    # Test imports and basic functionality
    print("\n" + "=" * 80)
    print("PHASE 2-3: TESTING PACKAGE IMPORTS & BASIC FUNCTIONALITY")
    print("=" * 80)
    
    import_tests = [
        ("Import Imaginable", test_import_imaginable),
        ("Import ROIable", test_import_roiable),
        ("Import LabelMapable", test_import_labelmapable),
        ("Import VTK converters", test_import_vtk_converters),
        ("Create basic image", test_basic_image_creation),
        ("Test image transformations", test_image_operations),
        ("Test ROI operations", test_roi_operations),
    ]
    
    for test_name, test_func in import_tests:
        try:
            result = test_func()
            status = "✓ PASS" if result else "✗ FAIL"
            print(f"{status}: {test_name}")
            results["phase2_unit_tests"].append({"name": test_name, "status": "pass" if result else "fail"})
        except Exception as e:
            print(f"✗ ERROR: {test_name} - {str(e)}")
            results["phase2_unit_tests"].append({"name": test_name, "status": "error", "error": str(e)})
    
    # Backward compatibility tests
    print("\n" + "=" * 80)
    print("BACKWARD COMPATIBILITY VERIFICATION")
    print("=" * 80)
    
    compat_tests = [
        ("Core API unchanged", test_core_api_compatibility),
        ("Image class hierarchy intact", test_class_hierarchy),
        ("Existing methods callable", test_existing_methods),
    ]
    
    for test_name, test_func in compat_tests:
        try:
            result = test_func()
            status = "✓ PASS" if result else "✗ FAIL"
            print(f"{status}: {test_name}")
            results["backward_compatibility"].append({"name": test_name, "status": "pass" if result else "fail"})
        except Exception as e:
            print(f"✗ ERROR: {test_name} - {str(e)}")
            results["backward_compatibility"].append({"name": test_name, "status": "error", "error": str(e)})
    
    # Summary
    print("\n" + "=" * 80)
    print("TEST SUMMARY")
    print("=" * 80)
    
    total_tests = len(results["phase1_fixes"]) + len(results["phase2_unit_tests"]) + len(results["backward_compatibility"])
    passed_tests = sum(1 for r in results["phase1_fixes"] + results["phase2_unit_tests"] + results["backward_compatibility"] 
                      if r["status"] == "pass")
    failed_tests = sum(1 for r in results["phase1_fixes"] + results["phase2_unit_tests"] + results["backward_compatibility"] 
                      if r["status"] == "fail")
    error_tests = sum(1 for r in results["phase1_fixes"] + results["phase2_unit_tests"] + results["backward_compatibility"] 
                     if r["status"] == "error")
    
    results["summary"] = {
        "total_tests": total_tests,
        "passed": passed_tests,
        "failed": failed_tests,
        "errors": error_tests,
        "success_rate": f"{100 * passed_tests / max(total_tests, 1):.1f}%"
    }
    
    print(f"Total Tests: {total_tests}")
    print(f"Passed: {passed_tests} ✓")
    print(f"Failed: {failed_tests} ✗")
    print(f"Errors: {error_tests} ⚠")
    print(f"Success Rate: {results['summary']['success_rate']}")
    
    return results


# ============================================================================
# VALIDATION FUNCTIONS FOR PHASE 1 FIXES
# ============================================================================

def validate_no_duplicates():
    """Check that duplicate vtk2sitk function is removed"""
    with open('/home/erosm/pyable/pyable_eros_montin/meshable.py', 'r') as f:
        content = f.read()
    # Count occurrences of function definition
    count = content.count('def vtk2sitk(')
    return count == 1

def validate_f_strings():
    """Check that f-strings are used for formatting"""
    with open('/home/erosm/pyable/pyable_eros_montin/imaginable.py', 'r') as f:
        content = f.read()
    # Check that bad pattern is gone and f-strings are used
    return 'f"Can\'t {message}"' in content or "f'Can't {message}'" in content

def validate_operators():
    """Check that bitwise & is replaced with 'and' in boolean context"""
    with open('/home/erosm/pyable/pyable_eros_montin/utilizers.py', 'r') as f:
        content = f.read()
    # Check that fix is applied
    return "(r is not None) and (t is not None)" in content

def validate_spelling():
    """Check that spelling errors are fixed"""
    with open('/home/erosm/pyable/pyable_eros_montin/imaginable.py', 'r') as f:
        content = f.read()
    # Check typos are fixed
    has_typos = "tyoe" in content or "shuld" in content or "pleae" in content
    has_fixes = "type" in content and "should" in content and "please" in content
    return (not has_typos) and has_fixes

def validate_mergelabels_scope():
    """Check that LABELMAP is initialized before use"""
    with open('/home/erosm/pyable/pyable_eros_montin/imaginable.py', 'r') as f:
        content = f.read()
    # Check for proper initialization
    return "LABELMAP = None" in content and "if LABELMAP is None:" in content

def validate_nan_comparison():
    """Check that NaN comparison uses np.isnan()"""
    with open('/home/erosm/pyable/pyable_eros_montin/imaginable.py', 'r') as f:
        content = f.read()
    # Check that bad pattern is gone and good one is present
    return "==np.NaN" not in content and "np.isnan(" in content

def validate_method_names():
    """Check that method name is fixed"""
    with open('/home/erosm/pyable/pyable_eros_montin/imaginable.py', 'r') as f:
        content = f.read()
    # Check that bad name is gone and good name is present
    return "getCenterOfGravitygetCenterOfGravityIndex" not in content and "getCenterOfGravityIndex" in content


# ============================================================================
# TEST FUNCTIONS FOR FUNCTIONALITY
# ============================================================================

def test_import_imaginable():
    """Test importing Imaginable class"""
    try:
        from pyable_eros_montin.imaginable import Imaginable, SITKImaginable
        return True
    except Exception as e:
        print(f"  Error: {e}")
        return False

def test_import_roiable():
    """Test importing Roiable class (lowercase 'oi')"""
    try:
        from pyable_eros_montin.imaginable import Roiable
        return True
    except Exception as e:
        print(f"  Error: {e}")
        return False

def test_import_labelmapable():
    """Test importing LabelMapable class"""
    try:
        from pyable_eros_montin.imaginable import LabelMapable
        return True
    except Exception as e:
        print(f"  Error: {e}")
        return False

def test_import_vtk_converters():
    """Test importing VTK converter functions"""
    try:
        from pyable_eros_montin.meshable import vtk2sitk, sitk2vtk
        return True
    except Exception as e:
        print(f"  Error: {e}")
        return False

def test_basic_image_creation():
    """Test creating basic image"""
    try:
        import SimpleITK as sitk
        from pyable_eros_montin.imaginable import SITKImaginable
        
        img = sitk.Image([50, 50, 50], sitk.sitkUInt8)
        imaginable = SITKImaginable(image=img)
        return imaginable.getImage() is not None
    except Exception as e:
        print(f"  Error: {e}")
        return False

def test_image_operations():
    """Test basic image operations"""
    try:
        import SimpleITK as sitk
        from pyable_eros_montin.imaginable import SITKImaginable
        
        img = sitk.Image([50, 50, 50], sitk.sitkFloat32)
        imaginable = SITKImaginable(image=img)
        
        # Test getting image
        ret_img = imaginable.getImage()
        
        # Test duplicate
        dup = imaginable.getDuplicate()
        
        return ret_img is not None and dup is not None
    except Exception as e:
        print(f"  Error: {e}")
        return False

def test_roi_operations():
    """Test ROI-specific operations"""
    try:
        import SimpleITK as sitk
        from pyable_eros_montin.imaginable import Roiable
        
        # Create non-empty image with label
        img = sitk.Image([50, 50, 50], sitk.sitkUInt8)
        for i in range(20, 30):
            for j in range(20, 30):
                for k in range(20, 30):
                    img.SetPixel([i, j, k], 1)
        
        roi = Roiable(image=img)
        
        # Test duplicate
        roi_dup = roi.getDuplicate()
        
        return roi_dup is not None and roi_dup.getImage() is not None
    except Exception as e:
        print(f"  Error: {e}")
        return False

def test_core_api_compatibility():
    """Test that core API remains unchanged"""
    try:
        from pyable_eros_montin.imaginable import (
            Imaginable, SITKImaginable, Roiable, 
            LabelMapable, Fieldable
        )
        
        # Check essential methods exist
        img_obj = SITKImaginable(image=__import__('SimpleITK').Image([10,10,10], __import__('SimpleITK').sitkUInt8))
        
        methods = [
            'getImage', 'getDuplicate', 'setImage', 'changePixelType',
            'rotateImage', 'scaleImage', 'translateImage',
            'getImageDimension', 'getImageSpacing'
        ]
        
        for method in methods:
            if not hasattr(img_obj, method):
                print(f"  Missing method: {method}")
                return False
        
        return True
    except Exception as e:
        print(f"  Error: {e}")
        return False

def test_class_hierarchy():
    """Test class hierarchy is intact"""
    try:
        from pyable_eros_montin.imaginable import (
            Imaginable, SITKImaginable, Roiable,
            LabelMapable, Fieldable
        )
        
        # Test inheritance
        import SimpleITK as sitk
        img = sitk.Image([10,10,10], sitk.sitkUInt8)
        
        # Create with non-empty data
        for i in range(2, 8):
            for j in range(2, 8):
                for k in range(2, 8):
                    img.SetPixel([i, j, k], 1)
        
        # All should instantiate without error
        sititk_obj = SITKImaginable(image=img)
        roi_obj = Roiable(image=img)
        
        # Roiable inherits from Imaginable (not directly SITKImaginable)
        return isinstance(roi_obj, Imaginable) and isinstance(sititk_obj, Imaginable)
    except Exception as e:
        print(f"  Error: {e}")
        return False

def test_existing_methods():
    """Test that existing methods still work"""
    try:
        import SimpleITK as sitk
        from pyable_eros_montin.imaginable import SITKImaginable
        
        img = sitk.Image([30, 30, 30], sitk.sitkUInt8)
        obj = SITKImaginable(image=img)
        
        # Test chaining operations
        dup = obj.getDuplicate()
        dup.changePixelType(sitk.sitkFloat32)
        
        # Test information retrieval
        dim = dup.getImageDimension()
        spacing = dup.getImage().GetSpacing()
        
        return dim == 3 and len(spacing) == 3
    except Exception as e:
        print(f"  Error: {e}")
        return False


if __name__ == '__main__':
    results = run_test_suite()
    
    print("\n" + "=" * 80)
    print("REGRESSION TEST COMPLETE")
    print("=" * 80)
    
    # Exit with appropriate code
    if results["summary"]["failed"] > 0 or results["summary"]["errors"] > 0:
        sys.exit(1)
    else:
        print("\n✓ All tests passed! Package is ready for deployment.")
        sys.exit(0)
