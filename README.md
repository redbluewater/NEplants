# NEplants
Plant project with Northeastern University\
Getting started, Krista Longnecker, 31 May 2023\

# 
4 June 2023\
Full set of commands in /KujawinskiLaboratory/UntargCode ...this is just highlights related to this project (which is only negative ion mode data) so I can easily get the commands I need. The full description remaines in UntargCode.

## Some steps before getting into the R/XCMS work
1. Convert the .RAW files from the mass spectrometer into mzML files using msConvert (use this script in R: ``ms_convert_tool_Lumos_v2.r``)
2. Use SCP to transfer those files to Poseidon (we are putting the files into this folder: /vortexfs1/omics/kujawinski/data)
3. Make a CSV file that contains the file names, ion mode, and good data markers (mtab_JBowen_PlantMetabolomics_052423.KL.csv).
4. Put this CSV file into the folder with the mzML files on Poseidon (again with SCP). 

## How to access Poseidon, WHOI's HPC computing environment
I used a Git Bash terminal window to log into poseidon:
```ssh username@poseidon.whoi.edu```
The password is my WHOI Active Directory password. You have to be logged into the WHOI VPN for this to work. 

Once you are logged into Poseidon, activate the conda module with ```module load anaconda/5.1```

## Moving around code - Windows 10 - GitHub - Poseidon (Krista's setup)
Edit files on my local computer and then use git bash to move them back into GitHub:

```git add -A```\
```git commit -am "Brief description goes here"``` (can use the bit in quotes to describe the update)\
```git push```\
(enter the passcode I use to get files to GitHub)

For later updates, just change to the folder for this repository (UntargCode) and then use this command to move the files from GitHub to the HPC:\
```git pull https://github.com/KujawinskiLaboratory/UntargCode.git``` or just ```git pull```

## Step 1: Create metadata
Set this up to send in ionMode as a variable so I don't have to edit all the slurm scripts each time I change ion mode\
```sbatch --export=ionMode="neg" scripts_dir/step1-metadata.slurm```

Check how many files you have 
```wc -l metadata_neg.txt```

Use this number in Step 2 to set the total number of array jobs that will be run.

## Step 2: peak picking and peak shape evaluation
```sbatch --export=ionMode="neg" scripts_dir/step2-xcms1.slurm```

## Step 3: combine picked peaks
To speed up peak picking, we performed peak picking as an array. Now combine into a single MS OnDisk object

```sbatch --export=ionMode="neg" scripts_dir/step3-xcms_combine.slurm```

## Step 4: perform retention time correction, grouping and fill peaks
```sbatch --export=ionMode="neg" scripts_dir/step4-xcms2.slurm```

## Step 5: Create an xset object 
```sbatch --export=ionMode="neg" scripts_dir/step5-create_xset.slurm```

## Step 6: Use CAMERA to create pseudospectra
In this repositroy, thie following has been modified to only process neg mode data.\
```sbatch scripts_dir/step6-camera.slurm```

#
## Description of the files
The final set of files generated by this analysis can be analyzed in R/MATLAB/Python, and/or submitted to GNPS for molecular networking analysis. 
The final set of files will in the /output_dir/xcms2/ folder and are described as follows:

* 'SargPatch_untarg_neg_aligned.csv'
  * This is the aligned set of mzRT features from the analysis, with one row per mzRT feature. In addition to information about mz value and retention time, there is one column for each sample with the peak areafor the sample
* 'SargPatch_untarg_neg_picked.csv'
  * This is an unaligned list of mzRT features, these are not grouped into peaks by samples but can be useful for troubleshooting
* 'camera_neg.csv'
  * This is the aligned set of mzRT features, with the adducts and isotopes determined from CAMERA. This is the best file to use in R/MATLAB/Python for downstream analysis
* 'edgelist_neg.csv'
  * This is a list of edges that is needed for the GNPS analysis.
* 'fName_featureQuant_afterCAMERA_neg.txt'
  * This is the aligned list of mzRT features, with CAMERA information, in a format that can be submitted to GNPS.
* 'fName_featureQuant_all_neg.txt'
  * This is a feature table, mzRT information and peak areas for samples, in a format that can be submitted to GNPS. This version includes all mzRT features.
* 'fName_featureQuant_combined_neg.txt'
  * This is feature table, mzRT information and peak areas for samples, in a format that can be submitted to GNPS. In this datafile, only features with MS2 information are included
* 'ms2spectra_all_neg.mgf'
  * This is an mgf file, which can be submitted to GNPS; this version contains every MS2 spectra 
* 'ms2spectra_combined_neg.mgf'
  * This is an mgf file, which can be submitted to GNPS; this version combines MS2 spectra based on the metric set at step 4.
 
#
Updated 12 October 2023 to start pushing commits from @redbluewater

#
Updated 13 November 2023 to calculate the elemental formulas.
