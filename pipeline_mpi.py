import os
import useful_funcs 
import bagpipes_burst
import bagpipes_continuity
import matplotlib as mpl
import bagpipes_results
import bagpipes_plotting
import plot_sort
from datetime import datetime
from multiprocessing import Pool, cpu_count

# Set matplotlib settings
mpl.rcParams['font.family'] = 'serif'
mpl.rc('text', usetex=True)
useful_funcs.plotting_style()

# Input parameters
catalogue_path = "/Users/katherineormerod/Documents/LJMU/WIDE/catalogues/EGS_GOODSN_UDS_combined.fits"
burst = True
continuity = True
burst_version = 1.3
continuity_version = 1.3
run = True
results = True
plot = True
sort = True

# Timing
start_time = datetime.now()

# Define burst function
def run_burst():
    model = "burst"
    run_type = f"burst_v{burst_version}"
    save_path = f"/Users/katherineormerod/Documents/LJMU/WIDE/catalogues/bagpipes_burst_v{burst_version}.fits"
    res_catalogue_path = save_path

    if run:
        bagpipes_burst.bagpipes_burst_fit(burst_version, catalogue_path)

    if results:
        bagpipes_results.bagpipes_results(model, burst_version, catalogue_path, save_path)

    if plot:
        bagpipes_plotting.plotting(res_catalogue_path, run_type=run_type, model=model)

    if sort:
        plot_sort.bagpipes_sort(run_type)

# Define continuity function
def run_continuity():
    model = "continuity"
    run_type = f"continuity_v{continuity_version}"
    save_path = f"/Users/katherineormerod/Documents/LJMU/WIDE/catalogues/bagpipes_continuity_v{continuity_version}.fits"
    res_catalogue_path = save_path

    if run:
        bagpipes_continuity.bagpipes_continuity(continuity_version, catalogue_path)

    if results:
        bagpipes_results.bagpipes_results(model, continuity_version, catalogue_path, save_path)

    if plot:
        bagpipes_plotting.plotting(res_catalogue_path, run_type=run_type, model=model)

    if sort:
        plot_sort.bagpipes_sort(run_type)

if __name__ == "__main__":
    # Get the number of CPU cores available
    num_cores = cpu_count()

    # Create a pool of worker processes
    with Pool(processes=num_cores) as pool:
        # Map the functions to the pool
        if burst:
            pool.apply_async(run_burst)

        if continuity:
            pool.apply_async(run_continuity)

        # Close the pool and wait for all tasks to complete
        pool.close()
        pool.join()

    # Timing
    end_time = datetime.now()
    time_taken = end_time - start_time
    print(f"Time taken: {time_taken}")
