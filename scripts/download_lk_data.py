import lightkurve as lk
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

print("Downloading data...")
search_result = lk.search_lightcurve("Kepler-10", author="Kepler", exptime=1800)
lc_collection = search_result.download_all()
lc_stitched = lc_collection.stitch().flatten(window_length=101).remove_nans()


CHUNK_SIZE_DAYS = 90 

time_array = lc_stitched.time.value
start_time = time_array[0]
relative_time = (time_array - start_time) % CHUNK_SIZE_DAYS
quarter_index = (time_array - start_time) // CHUNK_SIZE_DAYS

df = pd.DataFrame({
    'time_absolute': time_array,
    'time_relative': relative_time,
    'quarter_id': quarter_index.astype(int),
    'flux': lc_stitched.flux.value,
    'flux_err': lc_stitched.flux_err.value
})


df.to_csv("C:/Users/jlama/Documents/ms-spectra-analysis/data/raw/lightkurve/kepler10_quarters.csv", index=False)
print("Data saved to kepler10_quarters.csv")