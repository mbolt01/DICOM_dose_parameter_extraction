# DICOM dose parameter extraction
Script for automated extraction of dosimetric parameters from multiple RT DICOM files.
Note: The code is not optimised for speed and the directory structure must be as follows:

- Parent directory
     - Patient Folder 1
         - RS file
         - RD file
         - RP file
     - Patient Folder 2
         - RS file
         - RD file
         - RP file....

## Instructions for use
Contained within this repository is the Python module "dvh_analysis.py" which can be run on a directory to extract dosimetric parameters from DICOM RT files.
Dependencies are:
- pandas
- numpy
- tkinter
- dicompyler-core



## Example usage
```python
# import the module
import dvh_analysis

# generate and store results
my_res = dvh_analysis.dicom_dvh_stats_multi(stats='import', # allow import of file with list of dosimetric statistics to extract
                                            include_body=False, # excludes the body contour from DVH statistic generation (for speed)
                                            struct_labels=True,# allow import of file containing user-defined stucture labels to add to output
                                            save_df=True, # automatically save the generated pandas dataframe
                                            user_structures=None) # list of structures to analyse. None ==> analyses all structures.
```
