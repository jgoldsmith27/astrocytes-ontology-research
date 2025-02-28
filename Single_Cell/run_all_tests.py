import os
import time

def run_tests():
    """Run all tests with YOLO enthusiasm"""
    print("\n🚀🚀🚀 INITIATING YOLO TEST SEQUENCE 🚀🚀🚀")
    print("=" * 60)
    
    # Install dependencies if needed
    print("🔧 Ensuring dependencies are installed...")
    os.system("pip install scanpy anndata pandas rdflib numpy scipy matplotlib")
    
    # Run individual tests
    tests = [
        "test_read_data.py",
        "test_ontology_builder.py",
        "test_integration.py"
    ]
    
    results = []
    for test in tests:
        print("\n" + "=" * 60)
        print(f"🧪 Running {test}...")
        start_time = time.time()
        exit_code = os.system(f"python {test}")
        duration = time.time() - start_time
        success = exit_code == 0
        results.append((test, success, duration))
        print(f"{'✅' if success else '❌'} {test} {'PASSED' if success else 'FAILED'} in {duration:.2f}s")
    
    # Print summary
    print("\n" + "=" * 60)
    print("📊 TEST SUMMARY:")
    all_passed = all(r[1] for r in results)
    for test, success, duration in results:
        status = "✅ PASSED" if success else "❌ FAILED"
        print(f"{status} - {test} ({duration:.2f}s)")
    
    if all_passed:
        print("\n🎉🎉🎉 ALL TESTS PASSED! YOLO MODE ACTIVATED! 🎉🎉🎉")
        print("Your code is ready to process single-cell data without coexpression!")
    else:
        print("\n💣 SOME TESTS FAILED! CHECK THE LOGS!")
    
    print("=" * 60)

if __name__ == "__main__":
    run_tests() 