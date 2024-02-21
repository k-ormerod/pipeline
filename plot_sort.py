import os
import glob 
import shutil
from PyPDF2 import PdfMerger

def merge_pdf(directory, output):
 
    # Create a PdfFileMerger object
    merger = PdfMerger()

    # Get a list of PDF files in the directory 
    pdf_files = [f for f in os.listdir(directory) if f.endswith('.pdf')]

    # Sort the PDF files (optional, but often helpful)
    pdf_files.sort() 

    # Add PDFs to the merger 
    for filename in pdf_files:
        filepath = os.path.join(directory, filename)
        merger.append(open(filepath, 'rb'))

    # Write the merged PDF file
    with open(f'{output}.pdf', 'wb') as output_file:
        merger.write(output_file)

def bagpipes_sort(run_type):
    bagpipes_path = f"/Users/katherineormerod/Documents/LJMU/Codes/pipeline/pipes/plots/{run_type}/"
    own_plot_path = f"/Users/katherineormerod/Documents/LJMU/WIDE/Bagpipes/{run_type}/plots/"

    # list all files in the bagpipes path
    files = glob.glob(bagpipes_path + "*")
    for file in files:
        # move the files to the own plot path
        shutil.copy(file, own_plot_path)


    plot_types = ["corner", "fit", "sfh", "spec_and_properties", "spec_flam", "OIII", "optical", "spec_fnu_residual", "break"]

    # create subdirs for each plot type

    for plot_type in plot_types:
        if not os.path.exists(own_plot_path + plot_type):
            os.makedirs(own_plot_path + plot_type)

    # move the files to the correct subdirs
    
    files = glob.glob(own_plot_path + "*")
    for file in files:
    
        if "corner" in file:
            try:
                shutil.move(file, own_plot_path + "corner/")
            except:
                os.remove(file)
        elif "fit" in file:
            try:
                shutil.move(file, own_plot_path + "fit/")
            except:
                os.remove(file)  
        elif "sfh" in file:
            try:
                shutil.move(file, own_plot_path + "sfh/")
            except:
                os.remove(file)
        elif "spec_and_properties" in file:
            try:
                shutil.move(file, own_plot_path + "spec_and_properties/")
            except:
                os.remove(file)
        elif "spec_flam" in file:
            try:
                shutil.move(file, own_plot_path + "spec_flam/")
            except:
                os.remove(file)
        elif "OIII" in file:
            try:
                shutil.move(file, own_plot_path + "OIII/")
            except:
                os.remove(file)
        elif "optical" in file:
            try:
                shutil.move(file, own_plot_path + "optical/")
            except:
                os.remove(file)
        elif "spec_fnu_residual" in file:
            try:
                shutil.move(file, own_plot_path + "spec_fnu_residual/")
            except:
                os.remove(file)
        elif "break" in file:
            try:
                shutil.move(file, own_plot_path + "break/")
            except:
                os.remove(file)
        else:
            pass
    
    # merge the files in each subdir into one pdf
    for plot_type in plot_types:
        merge_pdf(own_plot_path + plot_type, own_plot_path + plot_type + "_compilation")

    # move the compilation files to comp subdir
        
    if not os.path.exists(own_plot_path + "comp"):
        os.makedirs(own_plot_path + "comp")

    files = glob.glob(own_plot_path + "*compilation.pdf")
    
    for file in files:
        try:
            shutil.move(file, own_plot_path + "comp/")
        except: 
            os.remove(file)



