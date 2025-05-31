import pygedm

print("Attempting a single call to pygedm.dist_to_dm with ymw16...")
try:
    # Using values that previously showed the error
    dm_info = pygedm.dist_to_dm(gl=0, gb=0, dist=1000, method='ymw16')
    print(f"Call successful. DM Info: {dm_info}")
except Exception as e:
    print(f"Call failed. Error type: {type(e)}")
    print(f"Error message: {e}")
    import traceback
    traceback.print_exc()

print("-" * 30)
print("Attempting a single call to pygedm.dist_to_dm with ne2001...")
try:
    # Using values that previously showed the error
    dm_info = pygedm.dist_to_dm(gl=0, gb=0, dist=1000, method='ne2001')
    print(f"Call successful. DM Info: {dm_info}")
except Exception as e:
    print(f"Call failed. Error type: {type(e)}")
    print(f"Error message: {e}")
    import traceback
    traceback.print_exc()

print("-" * 30)
print("Debug script finished.")
