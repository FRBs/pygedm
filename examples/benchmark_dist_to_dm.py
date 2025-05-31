import pygedm
import time
import sys # Added for sys.exit()

gl_values = [0, 45, 90, 135, 180, 225, 270, 315]  # degrees
gb_values = [-90, -45, 0, 45, 90]  # degrees
dist_values = [1000, 5000, 10000, 20000, 50000]  # pc

methods = ['ymw16', 'ne2001']
results = {method: [] for method in methods}

print("Starting benchmark for pygedm.dist_to_dm...")
print("Input values:")
print(f"  gl: {gl_values}")
print(f"  gb: {gb_values}")
print(f"  dist: {dist_values}")
print("-" * 30)

try:
    for gl_val in gl_values:
        for gb_val in gb_values:
            for dist_val in dist_values:
                for method_val in methods:
                    start_time = time.perf_counter()
                    try:
                        dm = pygedm.dist_to_dm(gl_val, gb_val, dist_val, method=method_val)
                        end_time = time.perf_counter()
                        execution_time = end_time - start_time
                        results[method_val].append(execution_time)
                        # Ensure dm is a tuple of two quantities for the print statement
                        if isinstance(dm, tuple) and len(dm) == 2 and hasattr(dm[0], 'value') and hasattr(dm[1], 'value'):
                            print(f"gl={gl_val:<3}, gb={gb_val:<3}, dist={dist_val:<5}, method='{method_val}', time={execution_time:.6f}s, DM={dm[0].value:.2f}, Tau_sc={dm[1].value:.6e}")
                        else:
                            print(f"gl={gl_val:<3}, gb={gb_val:<3}, dist={dist_val:<5}, method='{method_val}', time={execution_time:.6f}s, DM_INFO={dm!r}")

                    except TypeError as te:
                        end_time = time.perf_counter()
                        execution_time = end_time - start_time
                        error_message = str(te)
                        print(f"gl={gl_val:<3}, gb={gb_val:<3}, dist={dist_val:<5}, method='{method_val}', time={execution_time:.6f}s, ERROR: {error_message}")
                        if "unsupported format string passed to tuple.__format__" in error_message:
                            print(f"TARGET ERROR CAUGHT: gl={gl_val}, gb={gb_val}, dist={dist_val}, method='{method_val}'")
                            print("Exiting benchmark early.")
                            raise # Re-raise to stop the script
                        # For other TypeErrors, just report, they won't be added to results.
                    except Exception as e:
                        end_time = time.perf_counter()
                        execution_time = end_time - start_time
                        print(f"gl={gl_val:<3}, gb={gb_val:<3}, dist={dist_val:<5}, method='{method_val}', time={execution_time:.6f}s, ERROR: {e}")
                        import traceback
                        traceback.print_exc()
                        # Store a None or a specific error marker if needed, or just skip for averaging
                        # For now, errors will not be included in average time calculation
except Exception as final_e:
    print(f"Benchmark loop terminated due to an error: {final_e}")
    # This ensures that if we re-raise from the TypeError catch, it's handled here.

print("-" * 30)
print("Benchmark finished.")
print("-" * 30)
print("Average execution times:")

for method, times in results.items():
    if times:
        average_time = sum(times) / len(times)
        print(f"  Method '{method}': {average_time:.6f}s (based on {len(times)} successful calls)")
    else:
        print(f"  Method '{method}': No successful calls")

print("-" * 30)
print("Note: Errors during dist_to_dm calculation are reported but not included in average time.")
