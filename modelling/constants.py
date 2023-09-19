import os
import inspect
# Get main directory
frame = inspect.currentframe()
MAIN_DIR = os.path.abspath(os.path.join(os.path.dirname(
    inspect.getfile(frame)), '..'))

# Create results directory
RESULT_DIR = os.path.join(MAIN_DIR, 'simulation_data')
if not os.path.isdir(RESULT_DIR):
    os.makedirs(RESULT_DIR)

# Define parameters directory
PARAM_DIR = os.path.join(MAIN_DIR, 'parameters')
del(os, inspect, frame)

ABS_TOL = 1e-7
REL_TOL = 1e-8
