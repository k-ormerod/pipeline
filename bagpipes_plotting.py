import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib as mpl
from astropy import units as u
import os
from useful_funcs import plotting_style
from astropy.table import Table
from datetime import datetime
from astropy import constants as const

mpl.rc('text', usetex=True)

def plotting(catalogue_path, run_type, model):
    catalogue = Table.read(catalogue_path)

    for row in catalogue:
        ID = ID = str(int(row["NIRSpec ID"])).zfill(6)
        z = row["redshift"]

        z_spec = row["z_spec"]

        if not model == "continuity":
            z = z_spec
    

        try: 
            hdulist = fits.open(f"/Users/katherineormerod/Documents/LJMU/WIDE/prism_clear/{ID}_prism_clear_v3.1_1D.fits")
        except: 
            hdulist = fits.open(f"/Users/katherineormerod/Documents/LJMU/WIDE/prism_clear/{ID}_prism_clear_v3.0_1D.fits")

        wavelength = np.array(hdulist[9].data) * 1e10 * u.Angstrom
        flux = np.array(hdulist[1].data) * u.W / (u.m * u.m * u.m)
        flux_cgs = flux.to(u.erg / (u.s * u.cm * u.cm * u.Angstrom))

        posterior_path = f"/Users/katherineormerod/Documents/LJMU/WIDE/Bagpipes/{run_type}/spectra/{ID}_spectra.txt"

        if not os.path.exists(posterior_path):
            print(f"{ID} has no posterior spectrum file")
            continue

        posterior = np.genfromtxt(posterior_path, skip_header=1)

        # print(posterior)

        wave_post = posterior[:,0]
        spec_post = posterior[:,1]

        lower_lim = 8000
        upper_lim = 53000

        wave_post = wave_post * u.Angstrom
        spec_post_cgs = spec_post * (u.erg / (u.s * u.cm * u.cm * u.Angstrom))
        spec_post_fnu = (wave_post ** 2) * spec_post_cgs / const.c
        spec_post_fnu = spec_post_fnu.to(u.erg / (u.s * u.cm * u.cm * u.Hz))

        wave = wavelength
        flux_fnu = (wave ** 2) * flux_cgs / const.c   
        flux_fnu = flux_fnu.to(u.erg / (u.s * u.cm * u.cm * u.Hz))
        err = hdulist[2].data
        err = err * u.W / (u.m * u.m * u.m)
        err_cgs = err.to(u.erg / (u.s * u.cm * u.cm * u.Angstrom))
        err_fnu = (wave ** 2) * err_cgs / const.c
        err_fnu = err_fnu.to(u.erg / (u.s * u.cm * u.cm * u.Hz))




        indices = np.where((wave_post.value > lower_lim) & (wave_post.value < upper_lim))
        indices_obs = np.where((wavelength.value > lower_lim) & (wavelength.value < upper_lim))

        flux_regional = flux_fnu[indices_obs]
        spec_post_regional = spec_post_fnu[indices]
        err_region = err_fnu[indices_obs]
        wave_regional = wavelength[indices_obs]

        # print(len(flux_regional), len(spec_post_regional))


        min_length = min(len(flux_regional), len(spec_post_regional), len(err_region), len(wave_regional))

        # Truncate or interpolate lines to match the minimum length
        flux_regional = flux_regional[:min_length]
        spec_post_regional = spec_post_regional[:min_length]
        err_region = err_region[:min_length]
        wave_regional = wave_regional[:min_length]

        rest_wavelength = np.linspace(6000,53000, 5000) / (1 + z)

        fig = plt.figure(figsize=(8, 8))
        gs = mpl.gridspec.GridSpec(2, 1, height_ratios=[3, 1], hspace=0.)

        ax1 = plt.subplot(gs[0])
        ax2 = plt.subplot(gs[1], sharex=ax1)

        ax1.plot(wavelength, flux_fnu, alpha = 0.5, label = r'\rm{Observed spectrum}', color = 'mediumblue')
        ax1.plot(wave_post, spec_post_fnu, label = r'{\sc{Bagpipes}} \rm{fit}', color = 'orange')

        ax1.set_xlim(lower_lim, upper_lim)
        ax1.legend()
        if not model == "continuity":
            ax1.set_title(rf"${ID}$, $z_{{\rm{{spec}}}} = {z_spec}$, $z_{{\rm{{bagpipes}}}} = {z}$")
        else:
            ax1.set_title(rf"${ID}$, $z_{{\rm{{spec}}}} = {z_spec}$")

        ax1.set_ylabel(r"$\rm{Flux / erg s}^{-1} \rm{cm}^{-2} \rm{Hz}^{-1}$")

        ax1_top = ax1.twiny()
        ax1_top.set_xlim((8000/(1+z)), (53000/(1+z)))  # Set the same limits as the observed wavelength axis
        # ax1_top.set_xticks([])  # Hide the ticks for the secondary x-axis

        ax1_top.set_xlabel(r"$\rm{Rest\ Wavelength}\ /\ \AA$")

        # Plotting the rest wavelengths on the top x-axis
        rest_tick_locations = [2000, 3000, 4000, 5000, 6000]
        # rest_tick_labels = [str(int(label)) for label in rest_tick_locations]
        rest_tick_labels = [r"$\mathrm{" + str(int(label)) + "}$" for label in rest_tick_locations]
        ax1_top.set_xticks(rest_tick_locations)
        ax1_top.set_xticklabels(rest_tick_labels, fontsize=10)

        residuals = (flux_regional - spec_post_regional)/err_region
        ax2.plot(wave_regional, residuals, color='black', alpha = 0.7, lw=1)
        ax2.axhline(0, color='red', linestyle='--', alpha = 1,lw=1)
        ax2.fill_between(np.arange(0,100000,5000), -1, 1, color='blue', alpha=0.3, label=r"$\pm 1\sigma$", lw=0)
        ax2.fill_between(np.arange(0,100000,5000), -2, -1, color='blue', alpha=0.15, label=r"$\pm 2\sigma$",lw=0)
        ax2.fill_between(np.arange(0,100000,5000), 1, 2, color='blue', alpha=0.15,lw=0)

        ax2.set_xlabel(r"$\rm{Observed Wavelength}$ / \AA")
        ax2.set_ylabel(r"$\rm{Residual}$")
        ax2.legend()

        ax1.set_xlim(8000,53000)

        plt.subplots_adjust(hspace=0)
        plt.savefig(f"/Users/katherineormerod/Documents/LJMU/WIDE/Bagpipes/{run_type}/plots/{ID}_spec_fnu_residuals.pdf")
        plt.close()



    for row in catalogue:
        ID  = str(int(row["NIRSpec ID"])).zfill(6)
        z = row["redshift"]
        z_spec = row["z_spec"]


        try: 
            hdulist = fits.open(f"/Users/katherineormerod/Documents/LJMU/WIDE/prism_clear/{ID}_prism_clear_v3.1_1D.fits")
        except: 
            hdulist = fits.open(f"/Users/katherineormerod/Documents/LJMU/WIDE/prism_clear/{ID}_prism_clear_v3.0_1D.fits")

        wavelength = np.array(hdulist[9].data) * 1e10 * u.Angstrom
        flux = np.array(hdulist[1].data) * u.W / (u.m * u.m * u.m)
        flux_cgs = flux.to(u.erg / (u.s * u.cm * u.cm * u.Angstrom))

        posterior_path = f"/Users/katherineormerod/Documents/LJMU/WIDE/Bagpipes/{run_type}/spectra/{ID}_spectra.txt"
        posterior = np.genfromtxt(posterior_path, skip_header=1)

        wave_post = posterior[:,0]
        spec_post = posterior[:,1]

        lower_lim = 8000
        upper_lim = 53000

        wave_post = wave_post * u.Angstrom
        spec_post_cgs = spec_post * (u.erg / (u.s * u.cm * u.cm * u.Angstrom))
        spec_post_fnu = (wave_post ** 2) * spec_post_cgs / const.c
        spec_post_fnu = spec_post_fnu.to(u.erg / (u.s * u.cm * u.cm * u.Hz))

        wave = wavelength
        flux_fnu = (wave ** 2) * flux_cgs / const.c   
        flux_fnu = flux_fnu.to(u.erg / (u.s * u.cm * u.cm * u.Hz))
        err = hdulist[2].data
        err = err * u.W / (u.m * u.m * u.m)
        err_cgs = err.to(u.erg / (u.s * u.cm * u.cm * u.Angstrom))
        err_fnu = (wave ** 2) * err_cgs / const.c
        err_fnu = err_fnu.to(u.erg / (u.s * u.cm * u.cm * u.Hz))




        indices = np.where((wave_post.value > lower_lim) & (wave_post.value < upper_lim))
        indices_obs = np.where((wavelength.value > lower_lim) & (wavelength.value < upper_lim))

        flux_regional = flux_fnu[indices_obs]
        spec_post_regional = spec_post_fnu[indices]
        err_region = err_fnu[indices_obs]
        wave_regional = wavelength[indices_obs]

        # print(len(flux_regional), len(spec_post_regional))


        min_length = min(len(flux_regional), len(spec_post_regional), len(err_region), len(wave_regional))

        # Truncate or interpolate lines to match the minimum length
        flux_regional = flux_regional[:min_length]
        spec_post_regional = spec_post_regional[:min_length]
        err_region = err_region[:min_length]
        wave_regional = wave_regional[:min_length]


        fig = plt.figure(figsize=(8, 8))
        gs = mpl.gridspec.GridSpec(2, 1, height_ratios=[3, 1], hspace=0.)

        ax1 = plt.subplot(gs[0])
        ax2 = plt.subplot(gs[1], sharex=ax1)

        ax1.plot(wavelength, flux_fnu, alpha = 0.5, label = r'\rm{Observed spectrum}', color = 'mediumblue')
        ax1.plot(wave_post, spec_post_fnu, label = r'{\sc{Bagpipes}} \rm{fit}', color = 'orange')

        ax1.set_xlim(lower_lim, upper_lim)
        ax1.legend()
        if not model == "continuity":
            ax1.set_title(rf"${ID}$, $z_{{\rm{{spec}}}} = {z_spec}$, $z_{{\rm{{bagpipes}}}} = {z}$")
        else:
            ax1.set_title(rf"${ID}$, $z_{{\rm{{spec}}}} = {z_spec}$")

        ax1.set_ylabel(r"$\rm{Flux / erg s}^{-1} \rm{cm}^{-2} \rm{Hz}^{-1}$")

        residuals = (flux_regional - spec_post_regional)/err_region
        ax2.plot(wave_regional, residuals, color='black', alpha = 0.7, lw=1)
        ax2.axhline(0, color='red', linestyle='--', alpha = 1,lw=1)
        ax2.fill_between(np.arange(0,100000,5000), -1, 1, color='blue', alpha=0.3, label=r"$\pm 1\sigma$", lw=0)
        ax2.fill_between(np.arange(0,100000,5000), -2, -1, color='blue', alpha=0.15, label=r"$\pm 2\sigma$",lw=0)
        ax2.fill_between(np.arange(0,100000,5000), 1, 2, color='blue', alpha=0.15,lw=0)

        ax2.set_xlabel(r"$\rm{Observed Wavelength}$ / \AA")
        ax2.set_ylabel(r"$\rm{Residual}$")
        ax2.legend()

        spec_region = flux_regional[np.where((wave_post.value > 4507 * (1+z)) & (wave_post.value < 5507 * (1+z)))]

        ax1.set_xlim(4507 * (1+z), 5507 * (1+z))

        plt.subplots_adjust(hspace=0)
        plt.savefig(f"/Users/katherineormerod/Documents/LJMU/WIDE/Bagpipes/{run_type}/plots/{ID}_OIII.pdf")
        plt.close()

    for row in catalogue:
        ID = ID = str(int(row["NIRSpec ID"])).zfill(6)
        z = row["redshift"]
        z_spec = row["z_spec"]

        try: 
            hdulist = fits.open(f"/Users/katherineormerod/Documents/LJMU/WIDE/prism_clear/{ID}_prism_clear_v3.1_1D.fits")
        except: 
            hdulist = fits.open(f"/Users/katherineormerod/Documents/LJMU/WIDE/prism_clear/{ID}_prism_clear_v3.0_1D.fits")

        wavelength = np.array(hdulist[9].data) * 1e10 * u.Angstrom
        flux = np.array(hdulist[1].data) * u.W / (u.m * u.m * u.m)
        flux_cgs = flux.to(u.erg / (u.s * u.cm * u.cm * u.Angstrom))

        posterior_path = f"/Users/katherineormerod/Documents/LJMU/WIDE/Bagpipes/{run_type}/spectra/{ID}_spectra.txt"
        posterior = np.genfromtxt(posterior_path, skip_header=1)

        wave_post = posterior[:,0]
        spec_post = posterior[:,1]

        lower_lim = 8000
        upper_lim = 53000

        wave_post = wave_post * u.Angstrom
        spec_post_cgs = spec_post * (u.erg / (u.s * u.cm * u.cm * u.Angstrom))
        spec_post_fnu = (wave_post ** 2) * spec_post_cgs / const.c
        spec_post_fnu = spec_post_fnu.to(u.erg / (u.s * u.cm * u.cm * u.Hz))

        wave = wavelength
        flux_fnu = (wave ** 2) * flux_cgs / const.c   
        flux_fnu = flux_fnu.to(u.erg / (u.s * u.cm * u.cm * u.Hz))
        err = hdulist[2].data
        err = err * u.W / (u.m * u.m * u.m)
        err_cgs = err.to(u.erg / (u.s * u.cm * u.cm * u.Angstrom))
        err_fnu = (wave ** 2) * err_cgs / const.c
        err_fnu = err_fnu.to(u.erg / (u.s * u.cm * u.cm * u.Hz))




        indices = np.where((wave_post.value > lower_lim) & (wave_post.value < upper_lim))
        indices_obs = np.where((wavelength.value > lower_lim) & (wavelength.value < upper_lim))

        flux_regional = flux_fnu[indices_obs]
        spec_post_regional = spec_post_fnu[indices]
        err_region = err_fnu[indices_obs]
        wave_regional = wavelength[indices_obs]

        # print(len(flux_regional), len(spec_post_regional))


        min_length = min(len(flux_regional), len(spec_post_regional), len(err_region), len(wave_regional))

        # Truncate or interpolate lines to match the minimum length
        flux_regional = flux_regional[:min_length]
        spec_post_regional = spec_post_regional[:min_length]
        err_region = err_region[:min_length]
        wave_regional = wave_regional[:min_length]


        fig = plt.figure(figsize=(8, 8))
        gs = mpl.gridspec.GridSpec(2, 1, height_ratios=[3, 1], hspace=0.)

        ax1 = plt.subplot(gs[0])
        ax2 = plt.subplot(gs[1], sharex=ax1)

        ax1.plot(wavelength, flux_fnu, alpha = 0.5, label = r'\rm{Observed spectrum}', color = 'mediumblue')
        ax1.plot(wave_post, spec_post_fnu, label = r'{\sc{Bagpipes}} \rm{fit}', color = 'orange')

        ax1.set_xlim(lower_lim, upper_lim)
        ax1.legend()
        if not model == "continuity":
            ax1.set_title(rf"${ID}$, $z_{{\rm{{spec}}}} = {z_spec}$, $z_{{\rm{{bagpipes}}}} = {z}$")
        else:
            ax1.set_title(rf"${ID}$, $z_{{\rm{{spec}}}} = {z_spec}$")

        ax1.set_ylabel(r"$\rm{Flux / erg s}^{-1} \rm{cm}^{-2} \rm{Hz}^{-1}$")

        residuals = (flux_regional - spec_post_regional)/err_region
        ax2.plot(wave_regional, residuals, color='black', alpha = 0.7, lw=1)
        ax2.axhline(0, color='red', linestyle='--', alpha = 1,lw=1)
        ax2.fill_between(np.arange(0,100000,5000), -1, 1, color='blue', alpha=0.3, label=r"$\pm 1\sigma$", lw=0)
        ax2.fill_between(np.arange(0,100000,5000), -2, -1, color='blue', alpha=0.15, label=r"$\pm 2\sigma$",lw=0)
        ax2.fill_between(np.arange(0,100000,5000), 1, 2, color='blue', alpha=0.15,lw=0)

        ax2.set_xlabel(r"$\rm{Observed Wavelength}$ / \AA")
        ax2.set_ylabel(r"$\rm{Residual}$")
        ax2.legend()
        

        spec_region = flux_regional[np.where((wave_post.value > 1100 * (1+z)) & (wave_post.value < 1500 * (1+z)))]

        ax1.set_xlim(1100 * (1+z), 1500 * (1+z))
        try:
            ax1.set_ylim(-1.1 * abs(np.min(spec_region.value)), 1.1 * abs(np.max(spec_region.value)))
        except:
            continue

        plt.subplots_adjust(hspace=0)
        plt.savefig(f"/Users/katherineormerod/Documents/LJMU/WIDE/Bagpipes/{run_type}/plots/{ID}_break.pdf")
        plt.close()

    for row in catalogue:
        ID = ID = str(int(row["NIRSpec ID"])).zfill(6)
        z = row["redshift"]
        z_spec = row["z_spec"]

        try: 
            hdulist = fits.open(f"/Users/katherineormerod/Documents/LJMU/WIDE/prism_clear/{ID}_prism_clear_v3.1_1D.fits")
        except: 
            hdulist = fits.open(f"/Users/katherineormerod/Documents/LJMU/WIDE/prism_clear/{ID}_prism_clear_v3.0_1D.fits")

        wavelength = np.array(hdulist[9].data) * 1e10 * u.Angstrom
        flux = np.array(hdulist[1].data) * u.W / (u.m * u.m * u.m)
        flux_cgs = flux.to(u.erg / (u.s * u.cm * u.cm * u.Angstrom))

        posterior_path = f"/Users/katherineormerod/Documents/LJMU/WIDE/Bagpipes/{run_type}/spectra/{ID}_spectra.txt"
        posterior = np.genfromtxt(posterior_path, skip_header=1)

        wave_post = posterior[:,0]
        spec_post = posterior[:,1]

        lower_lim = 8000
        upper_lim = 53000

        wave_post = wave_post * u.Angstrom
        spec_post_cgs = spec_post * (u.erg / (u.s * u.cm * u.cm * u.Angstrom))
        spec_post_fnu = (wave_post ** 2) * spec_post_cgs / const.c
        spec_post_fnu = spec_post_fnu.to(u.erg / (u.s * u.cm * u.cm * u.Hz))

        wave = wavelength
        flux_fnu = (wave ** 2) * flux_cgs / const.c   
        flux_fnu = flux_fnu.to(u.erg / (u.s * u.cm * u.cm * u.Hz))
        err = hdulist[2].data
        err = err * u.W / (u.m * u.m * u.m)
        err_cgs = err.to(u.erg / (u.s * u.cm * u.cm * u.Angstrom))
        err_fnu = (wave ** 2) * err_cgs / const.c
        err_fnu = err_fnu.to(u.erg / (u.s * u.cm * u.cm * u.Hz))




        indices = np.where((wave_post.value > lower_lim) & (wave_post.value < upper_lim))
        indices_obs = np.where((wavelength.value > lower_lim) & (wavelength.value < upper_lim))

        flux_regional = flux_fnu[indices_obs]
        spec_post_regional = spec_post_fnu[indices]
        err_region = err_fnu[indices_obs]
        wave_regional = wavelength[indices_obs]

        # print(len(flux_regional), len(spec_post_regional))


        min_length = min(len(flux_regional), len(spec_post_regional), len(err_region), len(wave_regional))

        # Truncate or interpolate lines to match the minimum length
        flux_regional = flux_regional[:min_length]
        spec_post_regional = spec_post_regional[:min_length]
        err_region = err_region[:min_length]
        wave_regional = wave_regional[:min_length]


        fig = plt.figure(figsize=(8, 8))
        gs = mpl.gridspec.GridSpec(2, 1, height_ratios=[3, 1], hspace=0.)

        ax1 = plt.subplot(gs[0])
        ax2 = plt.subplot(gs[1], sharex=ax1)

        ax1.plot(wavelength, flux_fnu, alpha = 0.5, label = r'\rm{Observed spectrum}', color = 'mediumblue')
        ax1.plot(wave_post, spec_post_fnu, label = r'{\sc{Bagpipes}} \rm{fit}', color = 'orange')

        ax1.set_xlim(lower_lim, upper_lim)
        ax1.legend()
        if not model == "continuity":
            ax1.set_title(rf"${ID}$, $z_{{\rm{{spec}}}} = {z_spec}$, $z_{{\rm{{bagpipes}}}} = {z}$")
        else:
            ax1.set_title(rf"${ID}$, $z_{{\rm{{spec}}}} = {z_spec}$")

        ax1.set_ylabel(r"$\rm{Flux / erg s}^{-1} \rm{cm}^{-2} \rm{Hz}^{-1}$")

        residuals = (flux_regional - spec_post_regional)/err_region
        ax2.plot(wave_regional, residuals, color='black', alpha = 0.7, lw=1)
        ax2.axhline(0, color='red', linestyle='--', alpha = 1,lw=1)
        ax2.fill_between(np.arange(0,100000,5000), -1, 1, color='blue', alpha=0.3, label=r"$\pm 1\sigma$", lw=0)
        ax2.fill_between(np.arange(0,100000,5000), -2, -1, color='blue', alpha=0.15, label=r"$\pm 2\sigma$",lw=0)
        ax2.fill_between(np.arange(0,100000,5000), 1, 2, color='blue', alpha=0.15,lw=0)

        ax2.set_xlabel(r"$\rm{Observed Wavelength}$ / \AA")
        ax2.set_ylabel(r"$\rm{Residual}$")
        ax2.legend()
        

        

        ax1.set_xlim(4000 * (1+z), 7000 * (1+z))

        plt.subplots_adjust(hspace=0)
        plt.savefig(f"/Users/katherineormerod/Documents/LJMU/WIDE/Bagpipes/{run_type}/plots/{ID}_optical.pdf")
        plt.close()