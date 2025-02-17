# The script cuts .tif files saved on the F-Trap setup out of their folders and puts them into the parental folder.
# This way the file name is not too long for Fiji to open

import os
import shutil

folder_name = r'E:\MyFiles\Experimental data\smIVTT\20230808 - smIVTT over time GUV\500nM_MB-ATTO647N_T7RNAPolMix_HiScribe_1to11SybrGold_IE_GUV_0min_37d_OD1,0-0,2_gain_50'

# Loop through all folders inside folder_name
for subdir in os.listdir(folder_name):
    subdir_path = os.path.join(folder_name, subdir)

    # Only continue if it's a directory (i.e. not a file)
    if os.path.isdir(subdir_path):
        # print(subdir_path)
        # Find the .tif file inside the folder
        for filename in os.listdir(subdir_path):
            if filename.endswith('.tif'):
                # Copy the file to folder_name and delete the folder
                # print(filename)
                src_path = os.path.join(subdir_path, filename)
                dst_path = os.path.join(folder_name, filename)
                shutil.move("\\\\?\\" + src_path, dst_path)
                os.rmdir(subdir_path)
                break
