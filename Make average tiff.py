import numpy as np
import tifffile
import os

# Folder path for the original and averaged stacks
folder_path = r"E:\MyFiles\Experimental data\smIVTT\20231101 - smIVT 10 pM 37 degrees\10 pM 37 degrees 2 aliquots"

# File names
file_name1 = '500nM_MB-ATTO647N_10pM_BTx32DNA_T7RNAPolMix_HiScribe_50min_37_degrees_OD1,5-0,2_gain_50_2_MMStack_Pos0.ome'
file_name2 = '500nM_MB-ATTO647N_10pM_BTx32DNA_T7RNAPolMix_HiScribe_80min_37_degrees_OD1,5-0,2_gain_50_2_MMStack_Pos0.ome'
output_file_name = '500nM_MB-ATTO647N_10pM_BTx32DNA_T7RNAPolMix_HiScribe_60min_37_degrees_OD1,5-0,2_gain_50_2_MMStack_Pos0.ome'

# Full file paths
file_path1 = os.path.join(folder_path, file_name1 + ".tif")
file_path2 = os.path.join(folder_path, file_name2 + ".tif")
output_file_path = os.path.join(folder_path, output_file_name + ".tif")

# Load the original .tif stacks
stack1 = tifffile.imread(file_path1)
stack2 = tifffile.imread(file_path2)

# Check shapes
if stack1.shape != stack2.shape:
    raise ValueError("Stacks must have the same shape")

# Average the stacks
avg_stack = np.mean([stack1, stack2], axis=0).astype(stack1.dtype)

# Save the averaged stack
tifffile.imsave(output_file_path, avg_stack)