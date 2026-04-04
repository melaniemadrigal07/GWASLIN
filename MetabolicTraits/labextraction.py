import tifffile
import numpy as np
import matplotlib.pyplot as plt
from skimage import color
import os
import pandas as pd

# INPUT DIRECTORY (should be changed to where your merged images are stored) ===
input_dir = "/Volumes/MMMHD/GWASLin/Plate_1/merged_rgb/flattened_rgb_16bit_python"
output_csv = "mean_lab_summary.csv"

results = []

# LOOP THROUGH FILES
for filename in os.listdir(input_dir):
    if filename.endswith(".tif"):
        file_path = os.path.join(input_dir, filename)
        print(f"Processing {filename}...")

        #  Load and normalize image
        rgb_img = tifffile.imread(file_path)
        rgb_norm = rgb_img / 65535.0
        lab_img = color.rgb2lab(rgb_norm)

        # L* contrast check, to exclude hyphae from being used to extract values
        L_channel = lab_img[:, :, 0]
        L_contrast = L_channel.max() - L_channel.min()

        if L_contrast < 3:
            mean_lab_clean = lab_img.mean(axis=(0, 1))
            hyphae_mask = np.ones_like(L_channel, dtype=bool)
            hyphae_coverage = 0.0
        else:
            dark_threshold = np.percentile(L_channel, 10)
            hyphae_mask = L_channel >= dark_threshold
            lab_clean = lab_img[hyphae_mask]
            mean_lab_clean = lab_clean.mean(axis=0)
            hyphae_coverage = 100 * (~hyphae_mask).sum() / hyphae_mask.size

        #  Save result
        results.append({
            "Filename": filename,
            "L": mean_lab_clean[0],
            "A": mean_lab_clean[1],
            "B": mean_lab_clean[2],
        })

# EXPORT TO CSV 
df = pd.DataFrame(results)
df.to_csv(output_csv, index=False)
print("Done. Results saved to:", output_csv)
