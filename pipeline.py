import useful_funcs 
import bagpipes_burst
import bagpipes_continuity
import matplotlib as mpl
import bagpipes_results
import bagpipes_plotting
import plot_sort
from datetime import datetime


mpl.rcParams['font.family'] = 'serif'
mpl.rc('text', usetex=True)
useful_funcs.plotting_style()

catalogue_path = "/Users/katherineormerod/Documents/LJMU/WIDE/catalogues/EGS_GOODSN_UDS_combined.fits"

burst = False
continuity = True
burst_version = 1.1
continuity_version = 1.1

run = False
results = False
plot = False
sort = True

start_time = datetime.now()

if burst:
    model = "burst"
    run_type = f"burst_v{burst_version}"
    save_path = f"/Users/katherineormerod/Documents/LJMU/WIDE/catalogues/EGS_GOODSN_UDS_burst_v{burst_version}.fits"
    res_catalogue_path = save_path

    if run:
        bagpipes_burst.bagpipes_burst_fit(burst_version, catalogue_path)

    if results:
        bagpipes_results.bagpipes_results(model, burst_version, catalogue_path, save_path)

    if plot:
        bagpipes_plotting.plotting(res_catalogue_path, run_type = run_type, model = model)

    if sort:
        plot_sort.bagpipes_sort(run_type)


if continuity:
    model = "continuity"
    run_type = f"continuity_v{continuity_version}"
    save_path = f"/Users/katherineormerod/Documents/LJMU/WIDE/catalogues/EGS_GOODSN_UDS_continuity_v{continuity_version}.fits"
    res_catalogue_path = save_path

    if run:
        bagpipes_continuity.bagpipes_continuity(continuity_version, catalogue_path)

    if results:
        bagpipes_results.bagpipes_results(model, continuity_version, catalogue_path, save_path)

    if plot:
        bagpipes_plotting.plotting(res_catalogue_path, run_type = run_type, model = model)
    
    if sort:
        plot_sort.bagpipes_sort(run_type)


end_time = datetime.now()

time_taken = datetime.strptime(end_time, "%H:%M:%S") - datetime.strptime(start_time, "%H:%M:%S")

print(f"Time taken: {time_taken}")