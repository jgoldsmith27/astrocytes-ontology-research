import os
import sys
from test_read_data import create_test_h5ad, test_read_data
from test_ontology_builder import test_ontology_builder

def test_full_pipeline():
    """Test the entire pipeline from data reading to ontology creation"""
    print("\n🔥🔥🔥 FULL YOLO INTEGRATION TEST 🔥🔥🔥")
    
    # Step 1: Test read_data
    read_data_success = test_read_data()
    if not read_data_success:
        print("❌ Integration test failed at read_data stage")
        return False
    
    # Step 2: Test ontology builder
    ontology_success = test_ontology_builder()
    if not ontology_success:
        print("❌ Integration test failed at ontology builder stage")
        return False
    
    # Step 3: Test full pipeline with actual main.py
    print("\n🧪 Testing full pipeline execution...")
    try:
        # Create test data
        test_file = create_test_h5ad()
        
        # Modify INPUT_PATH in main.py temporarily
        with open("main.py", "r") as f:
            main_content = f.read()
        
        # Replace the input path with our test file
        modified_content = main_content.replace(
            'INPUT_PATH = "/Users/jacob/Desktop/Karolinska Research/astrocytes-ontology-research/single_cell_data_cleaned.h5ad"',
            f'INPUT_PATH = "{test_file}"'
        )
        
        with open("main_test.py", "w") as f:
            f.write(modified_content)
        
        # Run the modified main
        print("🧪 Running full pipeline with test data...")
        os.system("python main_test.py")
        
        # Verify outputs exist
        assert os.path.exists("YOLO_EXPRESSION.csv"), "CSV output not created"
        assert os.path.exists("YOLO_ONTOLOGY.ttl"), "TTL output not created"
        
        print("✅ Full pipeline executed successfully!")
        print("✅ CSV output created: YOLO_EXPRESSION.csv")
        print("✅ Ontology output created: YOLO_ONTOLOGY.ttl")
        
        print("\n🎉🎉🎉 INTEGRATION TEST PASSED! YOLO MODE SUCCESSFUL! 🎉🎉🎉")
        return True
        
    except Exception as e:
        print(f"❌ INTEGRATION TEST FAILED: {str(e)}")
        return False
    finally:
        # Clean up
        for file in [test_file, "main_test.py", "YOLO_EXPRESSION.csv", "YOLO_ONTOLOGY.ttl"]:
            if os.path.exists(file):
                os.remove(file)
        print("🧹 Cleaned up test files")

if __name__ == "__main__":
    test_full_pipeline() 