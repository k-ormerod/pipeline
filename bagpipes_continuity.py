import useful_funcs 
import bagpipes as pipes
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib as mpl
from astropy import units as u
import os
import matplotlib
from astropy.table import Table
from datetime import datetime
from astropy import constants as const
from astropy.cosmology import FlatLambdaCDM

mpl.rcParams['font.family'] = 'serif'
mpl.rc('text', usetex=True)
useful_funcs.plotting_style()

def bagpipes_continuity(version, catalogue_path, mask_lya = True):

    version = str(version)

    now = datetime.now()

    start_time = now.strftime("%H:%M:%S")

    model = "leja"
    disp_factor = 2

    catalogue = Table.read(catalogue_path, format='fits') #read in catalogue    
    run_type = f"continuity_v{version}" #name of run
    no_residuals = [] # list of galaxies that did not have residuals plotted

    def bin(spectrum, binn):
        """ Bins up two or three column spectral data by a specified factor. """

        binn = int(binn)
        nbins = len(spectrum) // binn
        binspec = np.zeros((nbins, spectrum.shape[1]))

        for i in range(binspec.shape[0]):
            spec_slice = spectrum[i*binn:(i+1)*binn, :]
            binspec[i, 0] = np.mean(spec_slice[:, 0])
            binspec[i, 1] = np.mean(spec_slice[:, 1])

            if spectrum.shape[1] == 3:
                binspec[i,2] = (1./float(binn)
                                *np.sqrt(np.sum(spec_slice[:, 2]**2)))

        return binspec


    def load_spec(ID):
        """ Loads spectroscopic data from file. """
        try:
            hdulist = fits.open(f"/Users/katherineormerod/Documents/LJMU/WIDE/prism_clear/{ID}_prism_clear_v3.1_1D.fits")
        except:
            hdulist = fits.open(f"/Users/katherineormerod/Documents/LJMU/WIDE/prism_clear/{ID}_prism_clear_v3.0_1D.fits")

        wavelength = np.array(hdulist[9].data) * 1e10 * u.Angstrom
        flux = np.array(hdulist[1].data) * u.W / (u.m * u.m * u.m)
        flux_cgs = flux.to(u.erg / (u.s * u.cm * u.cm * u.Angstrom))
        flux_err = np.array(hdulist[2].data) * u.W / (u.m * u.m * u.m)
        flux_err_cgs = flux_err.to(u.erg / (u.s * u.cm * u.cm * u.Angstrom))


        spectrum = np.c_[wavelength.value,
                        flux_cgs.value,
                        flux_err_cgs.value]

        spectrum = spectrum[~np.isnan(spectrum).any(axis=1)]

        if mask_lya:
            indices = np.where((spectrum[:,0] > ((1200)*(1+z))) & (spectrum[:,0] < ((1400)*(1+z))))
            spectrum[indices, 2] = 999
        else:
            pass

        print(spectrum)
        return bin(spectrum, 1)


    for row in catalogue:
        ID = str(int(row["NIRSpec ID"])).zfill(6) #make sure the format is 6 digits long
        print(ID)
        z = row["z_spec"]
        field = str(row["FIELD"])

        galaxy = pipes.galaxy(ID, load_spec, photometry_exists=False)

        path_to_delete = f"/Users/katherineormerod/pipes/posterior/{run_type}/{ID}.h5"
        if os.path.exists(path_to_delete):
            os.remove(path_to_delete)
            print(f".h5 file deleted for {ID}")
        else:
            print(f".h5 file does not exist for {ID}, continuing...")
        """
        Load prism dispersion and priors
        """


        dust = {}
        dust["type"] = "Calzetti"
        # dust["eta"] = 1.
        dust["Av"] = (0., 5.)

        nebular = {}
        # nebular["logU"] = -3.
        nebular['logU'] = (-3, -0.5)

        fit_instructions = {}
        fit_instructions["dust"] = dust
        fit_instructions["nebular"] = nebular
        # fit_instructions["t_bc"] = 0.01
        fit_instructions["redshift"] = z

        # bin_edges = calculate_time_bins(z, 20, 3, 10, 4)

        metallicity_prior = 'uniform'
        sfh = 'continuity'

        def calculate_bins(redshift, redshift_sfr_start=20, log_time=True, output_unit = 'yr', return_flat = False, num_bins=6, first_bin=10*u.Myr, second_bin=None, cosmo = FlatLambdaCDM(H0=70, Om0=0.3, Tcmb0=2.725)):
            time_observed = cosmo.lookback_time(redshift)
            time_sfr_start = cosmo.lookback_time(redshift_sfr_start)
            time_dif = abs(time_observed - time_sfr_start)
            if second_bin is not None:
                assert second_bin > first_bin, "Second bin must be greater than first bin"

            if second_bin is None:
                diff = np.linspace(np.log10(first_bin.to(output_unit).value), np.log10(time_dif.to(output_unit).value), num_bins)
            else:
                diff = np.linspace(np.log10(second_bin.to(output_unit).value), np.log10(time_dif.to(output_unit).value), num_bins-1)
            
            if not log_time:
                diff = 10**diff

            if return_flat:
                if second_bin is None:
                    return np.concatenate(([0],diff))
                else:
                    if log_time:
                        return np.concatenate([[0, np.log10(first_bin.to(output_unit).value)], diff])
                    else:
                        return np.concatenate([[0, first_bin.to(output_unit).value], diff])
            bins = []
            bins.append([0, np.log10(first_bin.to('year').value) if log_time else first_bin.to('year').value])
            if second_bin is not None:
                bins.append([np.log10(first_bin.to('year').value) if log_time else first_bin.to('year').value, np.log10(second_bin.to('year').value) if log_time else second_bin.to('year').value])
            
            for i in range(1, len(diff)):
                bins.append([diff[i-1], diff[i]])
            
            return  bins
        
        bin_list = list(calculate_bins(z, 20, False, 'Myr', True, 4, 3*u.Myr, 10*u.Myr))
            

        continuity = {}
        #continuity["age"] = (0.01, 15) # Gyr
        #continuity['age_prior'] = age_prior
        continuity["massformed"] = (5., 12.)  # Log_10 total stellar mass formed: M_Solar
        continuity['metallicity_prior'] = metallicity_prior
        if metallicity_prior == 'log_10':
            continuity["metallicity"] = (1e-06, 10.0)
        elif metallicity_prior == 'uniform':
            continuity["metallicity"] = (0.0, 0.5)

        continuity['bin_edges'] = bin_list[:-1]
        scale = 0
        if sfh == 'continuity':
            scale = 0.3
        if sfh == 'continuity_bursty':
            scale = 1.0

        for i in range(1, len(continuity["bin_edges"])-1):
            continuity["dsfr" + str(i)] = (-10., 10.)
            continuity["dsfr" + str(i) + "_prior"] = "student_t"
            continuity["dsfr" + str(i) + "_prior_scale"] = scale  # Defaults to this value as in Leja19, but can be set
            continuity["dsfr" + str(i) + "_prior_df"] = 2       # Defaults to this value as in Leja19, but can be set

        # for i in range(1, len(continuity["bin_edges"])-1):
        #     continuity["dsfr" + str(i)] = (-10., 10.)
        #     continuity["dsfr" + str(i) + "_prior"] = "student_t"
            #continuity["dsfr" + str(i) + "_prior_scale"] = 0.3  # Defaults to this value as in Leja19, but can be set
            #continuity["dsfr" + str(i) + "_prior_df"] = 2       # Defaults to this value as in Leja19, but can be set

        fit_instructions["continuity"] = continuity

        new_lsf = np.genfromtxt("/Users/katherineormerod/Documents/LJMU/point_source_lsf_clear_prism_QD3_i185_j85.csv", delimiter=",", skip_header=1)

        wave = new_lsf[:,0] * u.Angstrom
        delta_wave = new_lsf[:,1] * u.Angstrom


        resolution = (wave / delta_wave)/2.355

        fit_instructions["R_curve"] = np.c_[wave.value, resolution.value]

        fit = pipes.fit(galaxy, fit_instructions, run=f"{run_type}")
        now = datetime.now()

        current_time = now.strftime("%H:%M:%S")
        print(f"""Fitting spectrum for galaxy {ID}.
        The time is {current_time}.
        """)
        fit.fit(verbose = True) # don't need updates

        """
        Create and save plots
        """

        # bagpipes plotting

        fig = fit.plot_spectrum_posterior(save = True, show = False)
        fig = fit.plot_corner(save = True, show = False)
        fig = fit.plot_sfh_posterior(save = True, show = False)

        fig = plt.figure(figsize=(12, 7))
        gs = mpl.gridspec.GridSpec(7, 4, hspace=3., wspace=0.1)

        ax1 = plt.subplot(gs[:4, :])

        pipes.plotting.add_spectrum_posterior(fit, ax1)

        labels = ["sfr", "dust:Av", "stellar_mass", "ssfr"]

        post_quantities = dict(zip(labels, [fit.posterior.samples[l] for l in labels]))
        # plt.title(f"ID: {ID} z_spec = {z}, z_bagpipes = {np.median(fit.posterior.samples['redshift']):.3f}")
        plt.title(f"ID: {ID} z_spec = {z}")
        axes = []
        for i in range(4):
            axes.append(plt.subplot(gs[4:, i]))
            pipes.plotting.hist1d(post_quantities[labels[i]], axes[-1], smooth=True, label=labels[i])

        directory_path = f"/Users/katherineormerod/Documents/LJMU/WIDE/Bagpipes/{run_type}/plots/"
        file_name = f"{ID}_spec_and_properties.pdf"
        file_path = os.path.join(directory_path, file_name)


        if not os.path.exists(directory_path):
            os.makedirs(directory_path)

        
        plt.savefig(file_path)
        plt.close()

        print(f"Saving spec and properties plot for {ID}")


        #own plotting
        #load in posterior spectrum

        spec_post = np.percentile(fit.posterior.samples["spectrum"], 50, axis=0)
        spec_upper = np.percentile(fit.posterior.samples["spectrum"], 16, axis=0)
        spec_lower = np.percentile(fit.posterior.samples["spectrum"], 84, axis=0)

        wave_post = fit.posterior.galaxy.spec_wavs


        # open the observed spectrum
        try:
            hdulist = fits.open(f"/Users/katherineormerod/Documents/LJMU/WIDE/prism_clear/{ID}_prism_clear_v3.1_1D.fits")
        except:
            hdulist = fits.open(f"/Users/katherineormerod/Documents/LJMU/WIDE/prism_clear/{ID}_prism_clear_v3.0_1D.fits")

        # conversions of the observed data

        wavelength = np.array(hdulist[9].data) * 1e10 * u.Angstrom
        flux = np.array(hdulist[1].data) * u.W / (u.m * u.m * u.m)
        flux_cgs = flux.to(u.erg / (u.s * u.cm * u.cm * u.Angstrom))
        err = hdulist[2].data
        err = err * u.W / (u.m * u.m * u.m)
            
        err_cgs = err.to(u.erg / (u.s * u.cm * u.cm * u.Angstrom))

        wave = wavelength

        err_fnu = (wave ** 2) * err_cgs / const.c

        err_fnu = err_fnu.to(u.erg / (u.s * u.cm * u.cm * u.Hz))

        plt.plot(wavelength, flux_cgs, alpha = 0.5, label = 'Observed spectrum', color = 'mediumblue')
        plt.plot(wave_post, spec_post, label = r'{\sc{Bagpipes}} fit', color = 'orange')

        plt.title(f"ID: {ID}, z_spec = {z}")
        plt.xlabel(r"Wavelength / \AA")
        plt.ylabel(r"Flux / erg s$^{-1}$ cm$^{-2}$ \AA$^{-1}$")
        plt.legend()

        plt.xlim(8000, 53000)
        plt.savefig(f"/Users/katherineormerod/Documents/LJMU/WIDE/Bagpipes/{run_type}/plots/{ID}_spec_flam.pdf")
        plt.close()
        """
        Save posterior spectra and data
        """
        output_params = []
        output_medians = []
        output_lower = []
        output_upper = []

        for key in fit.posterior.samples.keys():
            if key == "sfh" or key == "spectrum" or key == "spectrum_full" or key == "uvj" or key == "dust_curve":
                print(f"{key} not added to list")
            else:
                parameter_median = np.median(fit.posterior.samples[key])
                parameter_lower = np.percentile(fit.posterior.samples[key], 16, axis=0)
                parameter_upper = np.percentile(fit.posterior.samples[key], 84, axis=0)
                output_params.append(key)
                output_medians.append(parameter_median)
                output_lower.append(parameter_lower)
                output_upper.append(parameter_upper)
                print(f"{key} added to list")
                print(f"median: {parameter_median}, lower: {parameter_lower}, upper: {parameter_upper}")

        headers = "parameter, median, lower, upper"

        output_data = list(zip(output_params, output_medians, output_lower, output_upper))



        data_filepath = f"/Users/katherineormerod/Documents/LJMU/WIDE/Bagpipes/{run_type}/posterior/"


        file_name = f"{ID}_posterior.txt"
        file_path = os.path.join(data_filepath, file_name)


        if not os.path.exists(data_filepath):
            os.makedirs(data_filepath)

        post_array_filepath = file_path

        # Save the data as a text file
        with open(post_array_filepath, 'w') as file:
            file.write("{:<30} {:<30} {:<30} {:<30}\n".format("parameter", "median", "lower", "upper"))
            for row in output_data:
                file.write("{:<30} {:<30} {:<30} {:<30}\n".format(str(row[0]), str(row[1]), str(row[2]), str(row[3])))


        spec_post_list = list(spec_post)
        spec_post_lower_list = list(spec_lower)
        spec_post_upper_list = list(spec_upper)
        wave_post_list = list(np.array(wave_post))

        data_array = np.c_[wave_post_list, spec_post_list, spec_post_lower_list, spec_post_upper_list]

        data_filepath = f"/Users/katherineormerod/Documents/LJMU/WIDE/Bagpipes/{run_type}/spectra/"
        headers = "wavelength, flux, flux_lower, flux_upper"


        file_name = f"{ID}_spectra.txt"
        file_path = os.path.join(data_filepath, file_name)


        if not os.path.exists(data_filepath):
            os.makedirs(data_filepath)

        np.savetxt(file_path, data_array, header=headers)


        print(f"{ID} complete")



    finish_time = now.strftime("%H:%M:%S")

    time_taken = datetime.strptime(finish_time, "%H:%M:%S") - datetime.strptime(start_time, "%H:%M:%S")


    print(f"Finished bagpipes fitting for field {field}. {len(catalogue)} galaxies fitted.")
    print(f"Time taken: {time_taken}")
    print(f"The following galaxies did not have plots saved with residuals:")
    for i in no_residuals:
        print(i)