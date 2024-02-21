from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib as mpl
from astropy import units as u
import os
from useful_funcs import plotting_style
from astropy.table import Table
from datetime import datetime
from astropy import constants as const
import numpy as np

def bagpipes_results(model, version, catalogue_path, save_path):

    version = str(version)

    run_type = f"{model}_v{version}"

    posterior_path = f"/Users/katherineormerod/Documents/LJMU/WIDE/Bagpipes/{run_type}/posterior/"

    catalogue = Table.read(catalogue_path)

    """
    Create the lists of each property
    Create the lists for the upper and lower limit of each set by the 16th and 84th percentile
    Iterate over each object and append the values to the lists
    """

    if model == "burst":

        burst_age = []
        burst_age_upper = []
        burst_age_lower = []

        burst_massformed = []
        burst_massformed_upper = []
        burst_massformed_lower = []

        burst_metalllicity = []
        burst_metalllicity_upper = []
        burst_metalllicity_lower = []

        dust_av = []
        dust_av_upper = []
        dust_av_lower = []

        nebular_logu = []
        nebular_logu_upper = []
        nebular_logu_lower = []

        noise_scaling = []
        noise_scaling_upper = []
        noise_scaling_lower = []

        redshift = []
        redshift_upper = []
        redshift_lower = []

        veldisp = []
        veldisp_upper = []
        veldisp_lower = []

        stellar_mass = []
        stellar_mass_upper = []
        stellar_mass_lower = []

        formed_mass = []
        formed_mass_upper = []
        formed_mass_lower = []

        sfr = []
        sfr_upper = []
        sfr_lower = []

        ssfr = []
        ssfr_upper = []
        ssfr_lower = []

        nsfr = []
        nsfr_upper = []
        nsfr_lower = []

        mass_weighted_age = []
        mass_weighted_age_upper = []
        mass_weighted_age_lower = []

        tform = []
        tform_upper = []
        tform_lower = []

        tquench = []
        tquench_upper = []
        tquench_lower = []

        mass_weighted_zmet = []
        mass_weighted_zmet_upper = []
        mass_weighted_zmet_lower = []

    elif model == "continuity":

        continuity_dsfr1 = []
        continuity_dsfr1_upper = []
        continuity_dsfr1_lower = []

        continuity_dsfr2 = []
        continuity_dsfr2_upper = []
        continuity_dsfr2_lower = []

        continuity_massformed = []
        continuity_massformed_upper = []
        continuity_massformed_lower = []

        continuity_metallicity = []
        continuity_metallicity_upper = []
        continuity_metallicity_lower = []

        dust_av = []
        dust_av_upper = []
        dust_av_lower = []

        nebular_logu = []
        nebular_logu_upper = []
        nebular_logu_lower = []

        stellar_mass = []
        stellar_mass_upper = []
        stellar_mass_lower = []

        formed_mass = []
        formed_mass_upper = []
        formed_mass_lower = []

        sfr = []
        sfr_upper = []
        sfr_lower = []

        ssfr = []
        ssfr_upper = []
        ssfr_lower = []

        nsfr = []
        nsfr_upper = []
        nsfr_lower = []

        mass_weighted_age = []
        mass_weighted_age_upper = []
        mass_weighted_age_lower = []

        tform = []
        tform_upper = []
        tform_lower = []

        tquench = []
        tquench_upper = []
        tquench_lower = []

        mass_weighted_zmet = []
        mass_weighted_zmet_upper = []
        mass_weighted_zmet_lower = []

    else:
        print("Please select a model.")



    for row in catalogue:
        ID = ID = str(int(row["NIRSpec ID"])).zfill(6)

        object_results_path = posterior_path + f"{ID}_posterior.txt"

        if not os.path.exists(object_results_path):
            print(f"Results for {ID} not found")
        else:
            print(f"Adding results for {ID}")

        object_data = np.genfromtxt(object_results_path, skip_header=1)

        print(object_data)


        """
        Median = row[1]
        Upper = row[3]
        Lower = row[2]
        """

        if model == "burst":

            burst_age.append(object_data[0][1])
            burst_age_upper.append(object_data[0][3])
            burst_age_lower.append(object_data[0][2])

            burst_massformed.append(object_data[1][1])
            burst_massformed_upper.append(object_data[1][3])
            burst_massformed_lower.append(object_data[1][2])

            burst_metalllicity.append(object_data[2][1])
            burst_metalllicity_upper.append(object_data[2][3])
            burst_metalllicity_lower.append(object_data[2][2])

            dust_av.append(object_data[3][1])
            dust_av_upper.append(object_data[3][3])
            dust_av_lower.append(object_data[3][2])

            nebular_logu.append(object_data[4][1])
            nebular_logu_upper.append(object_data[4][3])
            nebular_logu_lower.append(object_data[4][2])

            noise_scaling.append(object_data[5][1])
            noise_scaling_upper.append(object_data[5][3])
            noise_scaling_lower.append(object_data[5][2])

            redshift.append(object_data[6][1])
            redshift_upper.append(object_data[6][3])
            redshift_lower.append(object_data[6][2])

            veldisp.append(object_data[7][1])
            veldisp_upper.append(object_data[7][3])
            veldisp_lower.append(object_data[7][2])

            stellar_mass.append(object_data[8][1])
            stellar_mass_upper.append(object_data[8][3])
            stellar_mass_lower.append(object_data[8][2])

            formed_mass.append(object_data[9][1])
            formed_mass_upper.append(object_data[9][3])
            formed_mass_lower.append(object_data[9][2])

            sfr.append(object_data[10][1])
            sfr_upper.append(object_data[10][3])
            sfr_lower.append(object_data[10][2])

            ssfr.append(object_data[11][1])
            ssfr_upper.append(object_data[11][3])
            ssfr_lower.append(object_data[11][2])

            nsfr.append(object_data[12][1])
            nsfr_upper.append(object_data[12][3])
            nsfr_lower.append(object_data[12][2])

            mass_weighted_age.append(object_data[13][1])
            mass_weighted_age_upper.append(object_data[13][3])
            mass_weighted_age_lower.append(object_data[13][2])

            tform.append(object_data[14][1])
            tform_upper.append(object_data[14][3])
            tform_lower.append(object_data[14][2])

            tquench.append(object_data[15][1])
            tquench_upper.append(object_data[15][3])
            tquench_lower.append(object_data[15][2])

            mass_weighted_zmet.append(object_data[16][1])
            mass_weighted_zmet_upper.append(object_data[16][3])
            mass_weighted_zmet_lower.append(object_data[16][2])

        elif model == "continuity":

            continuity_dsfr1.append(object_data[0][1])
            continuity_dsfr1_upper.append(object_data[0][3])
            continuity_dsfr1_lower.append(object_data[0][2])

            continuity_dsfr2.append(object_data[1][1])
            continuity_dsfr2_upper.append(object_data[1][3])
            continuity_dsfr2_lower.append(object_data[1][2])

            continuity_massformed.append(object_data[2][1])
            continuity_massformed_upper.append(object_data[2][3])
            continuity_massformed_lower.append(object_data[2][2])

            continuity_metallicity.append(object_data[3][1])
            continuity_metallicity_upper.append(object_data[3][3])
            continuity_metallicity_lower.append(object_data[3][2])

            dust_av.append(object_data[4][1])
            dust_av_upper.append(object_data[4][3])
            dust_av_lower.append(object_data[4][2])

            nebular_logu.append(object_data[5][1])
            nebular_logu_upper.append(object_data[5][3])
            nebular_logu_lower.append(object_data[5][2])

            stellar_mass.append(object_data[6][1])
            stellar_mass_upper.append(object_data[6][3])
            stellar_mass_lower.append(object_data[6][2])

            formed_mass.append(object_data[7][1])
            formed_mass_upper.append(object_data[7][3])
            formed_mass_lower.append(object_data[7][2])

            sfr.append(object_data[8][1])
            sfr_upper.append(object_data[8][3])
            sfr_lower.append(object_data[8][2])

            ssfr.append(object_data[9][1])
            ssfr_upper.append(object_data[9][3])
            ssfr_lower.append(object_data[9][2])

            nsfr.append(object_data[10][1])
            nsfr_upper.append(object_data[10][3])
            nsfr_lower.append(object_data[10][2])

            mass_weighted_age.append(object_data[11][1])
            mass_weighted_age_upper.append(object_data[11][3])
            mass_weighted_age_lower.append(object_data[11][2])

            tform.append(object_data[12][1])
            tform_upper.append(object_data[12][3])
            tform_lower.append(object_data[12][2])

            tquench.append(object_data[13][1])
            tquench_upper.append(object_data[13][3])
            tquench_lower.append(object_data[13][2])

            mass_weighted_zmet.append(object_data[14][1])
            mass_weighted_zmet_upper.append(object_data[14][3])
            mass_weighted_zmet_lower.append(object_data[14][2])

        else:
            print("Please select a model.")



        print(f"Results for {ID} added to list. Moving on to next object")


    """
    Now update the catalogue
    If the column already exists, it will be overwritten
    """

    def process_column(hdul, column_name, new_column_data, output_filename, update_existing=True):
        """
        This function will either replace an existing column, add a new one,
        or create a new table. Specify if the update mode is desired 
        or the addition of a brand new table with a new HDUList.

        Args:
            hdul (astropy.io.fits.HDUList): The input FITS HDUList.
            column_name (str): The name of the column to process.
            new_column_data (array-like): The new data for the column.
            output_filename (str): The name of the output FITS file.
            update_existing (bool): If True, updates the table in the existing file.
                                    If False, creates a new table in the output file.
        """
        orig_table = hdul[1].data 
        new_column_data = np.array(new_column_data)

        column_names = orig_table.columns.names
        column_index = column_names.index(column_name) if column_name in column_names else -1

        if update_existing:  # Update existing table
            if column_index != -1:
                orig_table[column_name] = new_column_data
            else:
                new_col = fits.Column(name=column_name, format='D', array=new_column_data)
                hdul[1].data = fits.BinTableHDU.from_columns(orig_table.columns + fits.ColDefs([new_col])).data
            hdul.writeto(output_filename, overwrite=True) 

        else:  # Create new table
            # Create a new FITS BinTableHDU (assume only one table extension)
            if column_index != -1: 
                new_table = fits.BinTableHDU(orig_table) 
                new_table.data[column_name] = new_column_data
            else:
                new_col = fits.Column(name=column_name, format='D', array=new_column_data)
                new_table = fits.BinTableHDU.from_columns(fits.ColDefs([new_col])) 
    
            # Optional: Copy header from original file, adding appropriate updates
            new_table.header = hdul[1].header.copy() 
            # ...(update header elements as needed) 

            # Create a new HDUList and save
            new_hdul = fits.HDUList([hdul[0], new_table])  # Assuming primary HDU at index 0
            new_hdul.writeto(output_filename, overwrite=True)

    with fits.open(catalogue_path) as hdul:

        if model == "burst":
            process_column(hdul, 'burst_age', burst_age, output_filename=save_path, update_existing=False)
            process_column(hdul, 'burst_age_upper', burst_age_upper, output_filename=save_path, update_existing=False)
            process_column(hdul, 'burst_age_lower', burst_age_lower, output_filename=save_path, update_existing=False)

            process_column(hdul, 'burst_massformed', burst_massformed, output_filename=save_path, update_existing=False)
            process_column(hdul, 'burst_massformed_upper', burst_massformed_upper, output_filename=save_path, update_existing=False)
            process_column(hdul, 'burst_massformed_lower', burst_massformed_lower, output_filename=save_path, update_existing=False)

            process_column(hdul, 'burst_metalllicity', burst_metalllicity, output_filename=save_path, update_existing=False)
            process_column(hdul, 'burst_metalllicity_upper', burst_metalllicity_upper, output_filename=save_path, update_existing=False)
            process_column(hdul, 'burst_metalllicity_lower', burst_metalllicity_lower, output_filename=save_path, update_existing=False)

            process_column(hdul, 'dust_av', dust_av, output_filename=save_path, update_existing=False)
            process_column(hdul, 'dust_av_upper', dust_av_upper, output_filename=save_path, update_existing=False)
            process_column(hdul, 'dust_av_lower', dust_av_lower, output_filename=save_path, update_existing=False)

            process_column(hdul, 'nebular_logu', nebular_logu, output_filename=save_path, update_existing=False)
            process_column(hdul, 'nebular_logu_upper', nebular_logu_upper, output_filename=save_path, update_existing=False)
            process_column(hdul, 'nebular_logu_lower', nebular_logu_lower, output_filename=save_path, update_existing=False)

            process_column(hdul, 'noise_scaling', noise_scaling, output_filename=save_path, update_existing=False)
            process_column(hdul, 'noise_scaling_upper', noise_scaling_upper, output_filename=save_path, update_existing=False)
            process_column(hdul, 'noise_scaling_lower', noise_scaling_lower, output_filename=save_path, update_existing=False)

            process_column(hdul, 'redshift', redshift, output_filename=save_path, update_existing=False)
            process_column(hdul, 'redshift_upper', redshift_upper, output_filename=save_path, update_existing=False)
            process_column(hdul, 'redshift_lower', redshift_lower, output_filename=save_path, update_existing=False)

            process_column(hdul, 'veldisp', veldisp, output_filename=save_path, update_existing=False)
            process_column(hdul, 'veldisp_upper', veldisp_upper, output_filename=save_path, update_existing=False)
            process_column(hdul, 'veldisp_lower', veldisp_lower, output_filename=save_path, update_existing=False)

            process_column(hdul, 'stellar_mass', stellar_mass, output_filename=save_path, update_existing=False)
            process_column(hdul, 'stellar_mass_upper', stellar_mass_upper, output_filename=save_path, update_existing=False)
            process_column(hdul, 'stellar_mass_lower', stellar_mass_lower, output_filename=save_path, update_existing=False)

            process_column(hdul, 'formed_mass', formed_mass, output_filename=save_path, update_existing=False)
            process_column(hdul, 'formed_mass_upper', formed_mass_upper, output_filename=save_path, update_existing=False)
            process_column(hdul, 'formed_mass_lower', formed_mass_lower, output_filename=save_path, update_existing=False)

            process_column(hdul, 'sfr', sfr, output_filename=save_path, update_existing=False)
            process_column(hdul, 'sfr_upper', sfr_upper, output_filename=save_path, update_existing=False)
            process_column(hdul, 'sfr_lower', sfr_lower, output_filename=save_path, update_existing=False)

            process_column(hdul, 'ssfr', ssfr, output_filename=save_path, update_existing=False)
            process_column(hdul, 'ssfr_upper', ssfr_upper, output_filename=save_path, update_existing=False)
            process_column(hdul, 'ssfr_lower', ssfr_lower, output_filename=save_path, update_existing=False)

            process_column(hdul, 'nsfr', nsfr, output_filename=save_path, update_existing=False)
            process_column(hdul, 'nsfr_upper', nsfr_upper, output_filename=save_path, update_existing=False)
            process_column(hdul, 'nsfr_lower', nsfr_lower, output_filename=save_path, update_existing=False)

            process_column(hdul, 'mass_weighted_age', mass_weighted_age, output_filename=save_path, update_existing=False)
            process_column(hdul, 'mass_weighted_age_upper', mass_weighted_age_upper, output_filename=save_path, update_existing=False)
            process_column(hdul, 'mass_weighted_age_lower', mass_weighted_age_lower, output_filename=save_path, update_existing=False)

            process_column(hdul, 'tform', tform, output_filename=save_path, update_existing=False)
            process_column(hdul, 'tform_upper', tform_upper, output_filename=save_path, update_existing=False)
            process_column(hdul, 'tform_lower', tform_lower, output_filename=save_path, update_existing=False)

            process_column(hdul, 'tquench', tquench, output_filename=save_path, update_existing=False)
            process_column(hdul, 'tquench_upper', tquench_upper, output_filename=save_path, update_existing=False)
            process_column(hdul, 'tquench_lower', tquench_lower, output_filename=save_path, update_existing=False)

            process_column(hdul, 'mass_weighted_zmet', mass_weighted_zmet, output_filename=save_path, update_existing=False)
            process_column(hdul, 'mass_weighted_zmet_upper', mass_weighted_zmet_upper, output_filename=save_path, update_existing=False)
            process_column(hdul, 'mass_weighted_zmet_lower', mass_weighted_zmet_lower, output_filename=save_path, update_existing=False)

        elif model == "continuity":

            process_column(hdul, 'continuity_dsfr1', continuity_dsfr1, output_filename=save_path, update_existing=False)
            process_column(hdul, 'continuity_dsfr1_upper', continuity_dsfr1_upper, output_filename=save_path, update_existing=False)
            process_column(hdul, 'continuity_dsfr1_lower', continuity_dsfr1_lower, output_filename=save_path, update_existing=False)

            process_column(hdul, 'continuity_dsfr2', continuity_dsfr2, output_filename=save_path, update_existing=False)
            process_column(hdul, 'continuity_dsfr2_upper', continuity_dsfr2_upper, output_filename=save_path, update_existing=False)
            process_column(hdul, 'continuity_dsfr2_lower', continuity_dsfr2_lower, output_filename=save_path, update_existing=False)

            process_column(hdul, 'continuity_massformed', continuity_massformed, output_filename=save_path, update_existing=False)
            process_column(hdul, 'continuity_massformed_upper', continuity_massformed_upper, output_filename=save_path, update_existing=False)
            process_column(hdul, 'continuity_massformed_lower', continuity_massformed_lower, output_filename=save_path, update_existing=False)

            process_column(hdul, 'continuity_metallicity', continuity_metallicity, output_filename=save_path, update_existing=False)
            process_column(hdul, 'continuity_metallicity_upper', continuity_metallicity_upper, output_filename=save_path, update_existing=False)
            process_column(hdul, 'continuity_metallicity_lower', continuity_metallicity_lower, output_filename=save_path, update_existing=False)

            process_column(hdul, 'dust_av', dust_av, output_filename=save_path, update_existing=False)
            process_column(hdul, 'dust_av_upper', dust_av_upper, output_filename=save_path, update_existing=False)
            process_column(hdul, 'dust_av_lower', dust_av_lower, output_filename=save_path, update_existing=False)

            process_column(hdul, 'nebular_logu', nebular_logu, output_filename=save_path, update_existing=False)
            process_column(hdul, 'nebular_logu_upper', nebular_logu_upper, output_filename=save_path, update_existing=False)
            process_column(hdul, 'nebular_logu_lower', nebular_logu_lower, output_filename=save_path, update_existing=False)

            process_column(hdul, 'stellar_mass', stellar_mass, output_filename=save_path, update_existing=False)
            process_column(hdul, 'stellar_mass_upper', stellar_mass_upper, output_filename=save_path, update_existing=False)
            process_column(hdul, 'stellar_mass_lower', stellar_mass_lower, output_filename=save_path, update_existing=False)

            process_column(hdul, 'formed_mass', formed_mass, output_filename=save_path, update_existing=False)
            process_column(hdul, 'formed_mass_upper', formed_mass_upper, output_filename=save_path, update_existing=False)
            process_column(hdul, 'formed_mass_lower', formed_mass_lower, output_filename=save_path, update_existing=False)

            process_column(hdul, 'sfr', sfr, output_filename=save_path, update_existing=False)
            process_column(hdul, 'sfr_upper', sfr_upper, output_filename=save_path, update_existing=False)
            process_column(hdul, 'sfr_lower', sfr_lower, output_filename=save_path, update_existing=False)

            process_column(hdul, 'ssfr', ssfr, output_filename=save_path, update_existing=False)
            process_column(hdul, 'ssfr_upper', ssfr_upper, output_filename=save_path, update_existing=False)
            process_column(hdul, 'ssfr_lower', ssfr_lower, output_filename=save_path, update_existing=False)

            process_column(hdul, 'nsfr', nsfr, output_filename=save_path, update_existing=False)
            process_column(hdul, 'nsfr_upper', nsfr_upper, output_filename=save_path, update_existing=False)
            process_column(hdul, 'nsfr_lower', nsfr_lower, output_filename=save_path, update_existing=False)

            process_column(hdul, 'mass_weighted_age', mass_weighted_age, output_filename=save_path, update_existing=False)
            process_column(hdul, 'mass_weighted_age_upper', mass_weighted_age_upper, output_filename=save_path, update_existing=False)
            process_column(hdul, 'mass_weighted_age_lower', mass_weighted_age_lower, output_filename=save_path, update_existing=False)

            process_column(hdul, 'tform', tform, output_filename=save_path, update_existing=False)
            process_column(hdul, 'tform_upper', tform_upper, output_filename=save_path, update_existing=False)
            process_column(hdul, 'tform_lower', tform_lower, output_filename=save_path, update_existing=False)

            process_column(hdul, 'tquench', tquench, output_filename=save_path, update_existing=False)
            process_column(hdul, 'tquench_upper', tquench_upper, output_filename=save_path, update_existing=False)
            process_column(hdul, 'tquench_lower', tquench_lower, output_filename=save_path, update_existing=False)

            process_column(hdul, 'mass_weighted_zmet', mass_weighted_zmet, output_filename=save_path, update_existing=False)
            process_column(hdul, 'mass_weighted_zmet_upper', mass_weighted_zmet_upper, output_filename=save_path, update_existing=False)
            process_column(hdul, 'mass_weighted_zmet_lower', mass_weighted_zmet_lower, output_filename=save_path, update_existing=False)

        else:
            print("Please select a model.")



    print("Catalogue updated")