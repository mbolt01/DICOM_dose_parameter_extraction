
# coding: utf-8

# ## Process developed to analyse multiple DICOM files and extract DVH information
# All files must be stored as follows:
# - Parent directory
#     - Patient Folder 1
#         - RS file
#         - RD file
#         - RP file
#     - Patient Folder 2
#         - RS file
#         - RD file
#         - RP file....
# 
# The DVH analysis is basaed on the dicompyler-core module (which is based on pydicom).
# Notebook with examples: https://github.com/bastula/dicom-notebooks/blob/master/dicompyler-core_usage.ipynb
# 
# - The user can specify the parent directory and then all sub directories are analysed.
# - The user must specify which DVH statistics are required (Can be from an excel/text/csv file).
# - Additional structure labels can be added to allow simpler comparison between differnet plans which may have different original structure names. *Note that all structure names are converted to lowercase*.
# - Results are output in CSV format and contain patient id, plan name, structrure name and the dvh statistics requested.
# 

# In[1]:

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
#get_ipython().magic('matplotlib inline')

import itertools as it
from dicompylercore import dicomparser, dvhcalc#, dvh
import os
import pydicom

## for choosing a directory
import tkinter as tk
from tkinter import filedialog

## for exiting on error (no oncentra dose grid)
import sys


# In[2]:

def dir_open(title = 'Select Directory...'):
    """Get a directory string using the OS file chooser dialog."""
    root = tk.Tk()
    my_dir =  filedialog.askdirectory(title = title)
    #print (root.directory)
    root.withdraw() ## for some reason it freezes if the empty window is not open
    return my_dir

def file_open(title = 'Select File...'):
    """Get a file string using the OS file chooser dialog."""
    root = tk.Tk()
    my_file =  filedialog.askopenfilename(title = title)
    #print (root.filename)
    root.withdraw()
    return my_file

def file_save(title = 'Save File As...',filetypes = (("csv","*.csv"),("all files","*.*")), initialfile=''):
    root = tk.Tk()
    ftypes = [('CSV', '.csv'), ('All files', '*')]
    my_file =  filedialog.asksaveasfilename(title = title, filetypes = ftypes, defaultextension='.csv',initialfile=initialfile)
    root.withdraw()
    return my_file


# In[3]:

def save_df_to_csv(df,prompt=False,save_dir='',save_name = None):
    """Function to save the dataframe as csv.
    Option to prompt the user where to save.
    Defualt save location is within the parent directory of analysis"""
    
    if save_name == None:
        the_date = pd.Timestamp("now").strftime("%Y%m%d-%H%M") ## no need to inport date module as can use pandas
        default_save_name = 'results-' + str(the_date) + '.csv'
    
    if prompt==True:
        save_name = file_save(title = 'Save Results As...', initialfile = default_save_name)
        #print(save_name)
        if save_name == '': ## user closed before saving
            return
    else:
        save_name = os.path.join(save_dir,default_save_name)
        
    df.to_csv(save_name)
    print('Results saved:', save_name)


# In[4]:

def swap_dict(the_dict):
    """Swaps key,value pairs in dict. Must be unique values to work"""
    return {y:x for x,y in the_dict.items()}


# In[5]:

def dicom_type(file):
    """ will return 'rtss', 'rtdose', 'rtplan','ct' dependant on type of dicom file"""
    the_type = dicomparser.DicomParser(file).GetSOPClassUID()
    return the_type


# In[6]:

def get_rt_files(directory):
    """provide a directory and the rtdose, rtss, rtplan files will be returned as a dict"""
    ## get list of all files in directory
    all_files = os.listdir(directory)
    
    ## get the rt files and put tehm in a dict
    ## note: key will not exist if the correct file does not exist
    pt_rt_file_dict = {}
    for file in all_files:
        full_file_dir = os.path.join(directory,file)
        #print(full_file_dir)
        d_type = dicom_type(full_file_dir)
        if d_type == 'rtss':
            pt_rt_file_dict['rtss']=full_file_dir
        if d_type == 'rtdose':
            pt_rt_file_dict['rtdose']=full_file_dir
        if d_type == 'rtplan':
            pt_rt_file_dict['rtplan']=full_file_dir
    
    return pt_rt_file_dict


# In[7]:

def structure_dict(file):
    """Loop through (rtss) file to get {struct_name:struct_id}.
    The loop will go through the IDs to get the structure names.
    file is the complete file directory"""
    ## check to ensure rtss filetype before trying to get structure names
    file_type = dicom_type(file)
    if file_type != 'rtss':## alternative-check length of the returned dict. if 0 then not structs - best to check before
        print('File is not an RT Structure set:',file)
        return
    else:
        ## get the raw structure info from the dicom file - this returns a dict {number:{'color'...id....name etc}}
        struct_info = dicomparser.DicomParser(file).GetStructures()
        
        ## create dict to store the {strucure_name,:structure_id}
        struct_name_dict = {}
        
        ## go through structures and store name and id in dict
        for keys,values in struct_info.items():
            struct_name = struct_info[keys]['name'].lower() ## make all lowercase?
            struct_id = struct_info[keys]['id']
            #print(keys,struct_id,struct_name)
            #print(abc[keys]['id'])
            struct_name_dict[struct_name] = struct_id
        
    return struct_name_dict


# In[8]:

def get_dicom_plan_info(file):
    """Get info from the dicom plan file.
    Return results as a dictionary so can be easily expanded as required"""
    
    if dicom_type(file)!= 'rtplan':
        print('Should provide an RTPLAN filetype')
        dicom_info = None
    else:
        ds = pydicom.read_file(file)
    
        ## Try and keep this fuinction the same for all manufacturers. USe seperate functinos to get specific info.
        ## this uses the functionality of pydicom, rather than specifying the dicom tags.
        ## pydicom essentially looks up the names from a dictionary and then looks up the dicom tag
        ## https://github.com/pydicom/pydicom/blob/master/pydicom/_dicom_dict.py
        ## using names is more understanable to most people, but can also use dicom tags
        
        dicom_info = {} ## store results in a dictionary
        
        ## values specified by name
        dicom_info['study_date'] = ds.StudyDate ## study is the CT info
        dicom_info['study_time'] = ds.StudyTime
        
        if 'study_description' in ds: ## this doesnt seem to appear in Variseed, so have to check it exists before setting it
            dicom_info['study_description'] = ds.StudyDescription
        else:
            dicom_info['study_description'] = None

        dicom_info['patient_id'] = ds.PatientID
        dicom_info['patient_name'] = ds.PatientName
        dicom_info['plan_label'] = ds.RTPlanLabel ## plan name
        dicom_info['plan_date'] = ds.RTPlanDate
        dicom_info['plan_time'] = ds.RTPlanTime
        
        ## values can also be specified by dicom tag - need to use the .value to retrieve the actual value
        dicom_info['manufacturer'] = ds[0x0008, 0x0070].value
    
    return dicom_info


# In[9]:

def get_pt_id(file):
    """Get patient ID from dicom file"""
    pt_info = dicomparser.DicomParser(file).GetDemographics()
    return pt_info['id']


# In[10]:

def all_equal(iterable):
    """Returns True if all the elements are equal to each other.
    Based on itertools recipe: https://docs.python.org/3/library/itertools.html#recipes
    Used to ensure same patient IDs for generating DVH from files"""
    g = it.groupby(iterable)
    return next(g, True) and not next(g, False)


# In[11]:

def all_rt_files_exist(directory):
    """Checks if all required DICOM file types exist within suipplied directory.
    Should have rtss, rtplan, rtdose.
    Returns (True/False, List of files)"""
    rt_files = get_rt_files(directory)
    rt_files_keys = rt_files.keys()
   
    ## check if each filetype exists
    rt_non_exist = [] ## store ones which do not exist
    rt_exist = [] ## store ones which do exist
    if 'rtss' not in rt_files_keys:
        rt_non_exist.append('rtss')
    else:
        rt_exist.append('rtss')
    if 'rtdose' not in rt_files_keys:
        rt_non_exist.append('rtdose')
    else:
        rt_exist.append('rtdose')
    if 'rtplan' not in rt_files_keys:
        rt_non_exist.append('rtplan')
    else:
        rt_exist.append('rtplan')
    
    if len(rt_non_exist) != 0:
        print('The following filetypes do not exist:',rt_non_exist)
        all_exist = False
    else:
        all_exist = True
        
    ## extract the ids for the files which exist
    rt_files_exist_pt_ids = []
    for key in rt_exist:
        rt_pt_id = get_pt_id(rt_files[key])
        rt_files_exist_pt_ids.append(rt_pt_id)
    
    ## put in value to force failure for testing
    #rt_files_exist_pt_ids[2] = '333'
        
    ## check all the ids match
    if all_equal(rt_files_exist_pt_ids):
        all_ids_match = True
    else:
        all_ids_match = False
        print('Mis-matched patient IDs identified:',rt_files_exist_pt_ids)
        
    ## if all exist and match then return True
    if all_exist == True and all_ids_match == True:
        all_exist_match = True
    else:
        all_exist_match = False
    
    return all_exist_match,rt_files


# In[12]:

def get_dicom_dose(file):
    """Will return the prescribed dose (in Gy) from an rtplan file. This is needed for DVH analysis.
    An error will be shown to the user if a non RTPlan file is supplied"""
    if dicom_type(file) != 'rtplan':
        print('Prescription cannot be obtained from non RTPlan DICOM file. Returned None.')
        dicom_dose = None
    else:
        if is_oncentra(file) == True: ## oncentra seems to not store the prescribed dose in the same way, so deal with seperately
            print('Oncentra File - Attempting to determine prescription...')
            dicom_dose = get_oncentra_dose(file)
            ## get oncentra dose...
        else:
            dicom_dose = dicomparser.DicomParser(file).GetPlan()['rxdose']/100 ## convert to Gy for use
    return dicom_dose


# In[13]:

#my_plan = r'C:\Users\mb22\OneDrive\PhD\Quasar Shared\Data\DICOM\Oncentra\RP1.3.6.1.4.1.2452.6.3168945037.1244793830.2214002341.3012008818.dcm'
#my_dose = r'C:\Users\mb22\OneDrive\PhD\Quasar Shared\Data\DICOM\Oncentra\RD1.3.6.1.4.1.2452.6.2820584682.1258535538.2684740226.4060889106.dcm'
#my_ss = r'C:\Users\mb22\OneDrive\PhD\Quasar Shared\Data\DICOM\Oncentra\RS1.3.6.1.4.1.2452.6.2392950960.1269671022.97853884.238715541.dcm'


# In[14]:

def check_oncentra_dose_exists(file):
    dose_scaling_factor_tag = [0x3004,0x000e] ## use this to see if there is dose within an Oncentra file.
    try:
        dose_scale_exists = get_dicom_tag_info(file,dose_scaling_factor_tag)
        #print('Dose Grid Exists in Oncentra File')
        dose_exists = True
    except KeyError:
        #print('No Dose Grid within Oncentra File - do 3D Dose Grid Calculation')
        dose_exists = False
    return dose_exists


# In[15]:

def get_dicom_tag_info(file,tag):
    """Get info from dicom file using a dicom tag specified as a tuple (0x0000,0x0000)
    Reading the file each time a piece of info is needed is not particularly efficient, so try not to use this to much"""
    tag_value = pydicom.read_file(file)[tag].value
    return tag_value


# In[16]:

def is_oncentra(rt_plan):
    """Determine if the RTPLAN file indicates that it is from Oncentra.
    Used to allow prescribed dose to be determined from private tag"""
    manufacturer = pydicom.read_file(rt_plan)[0x00081090].value
    if manufacturer == 'Oncentra':
        oncentra = True ## probably want to look up this tag: (0008, 1090) Manufacturer's Model Name           LO: 'Oncentra'
    else:
        oncentra = False
    return oncentra

def get_oncentra_dose(my_file):
    """Get dose from private Oncentra Tag (3007, 1000)"""
    omp_dose_tag = (0x3007,0x1000)
    oncentra_dose = get_dicom_tag_info(my_file,omp_dose_tag)
    dose_exists = check_oncentra_dose_exists(my_file)
    if dose_exists == False:
        #print('3D dose grid was not calculated in Oncentra')
        sys.exit('No Dose Grid within Oncentra File - do 3D Dose Grid Calculation. Aborting Calculation')
    return oncentra_dose



# In[18]:

def calc_dvh(rtss,rtdose,structure_id,px):
    """Calculate a dvh object from input files and structure id, and assign the prescription to the dvh object.
    All 3 files must be provided."""
    dvh = dvhcalc.get_dvh(rtss,rtdose,structure_id)
    dvh.rx_dose = px
    return dvh


# In[19]:

def non_int_stat(dvh,stat_string):
    """Function to return the dvh statistic for integer and non-integer values"""
        
    ## get the stat type based on first letter (V/D)
    stat_string = stat_string.lower() ## make all lowercase for simplicity
    stat_type = stat_string[0]
    if stat_type not in ['v','d']:
        sys.exit('Unknown dose statistic type: ' + stat_type)
    
    ## get the units based on last letters (if they exist) (cc/Gy)
    ## get the end string (must reverse through the lsit to do this, then reverse it back to correct way)
    end_string = "".join(it.takewhile(str.isalpha, reversed(stat_string)))[::-1]
    if len(end_string) == 0:
        units = None
    else:
        units = end_string
        if units not in ['cc','gy']:
            sys.exit('Unknown dose statistic units: ' + units)
            
    ## get the value required - remove the stat type and the units to give the numerical value
    stat_numeric = float(stat_string[len(stat_type):len(stat_string)-len(end_string)])
    
    ## now if D.. then use dose_constraint, if V... use volume_constraint with the obtained units
    if stat_type =='v':
        stat = dvh.volume_constraint(stat_numeric, units)
        
    if stat_type =='d':
        stat = dvh.dose_constraint(stat_numeric, units)

    return stat


# In[20]:

def get_statistics(dvh,stat_list):
    """get multiple statistics from a dvh and store in a dict for use"""
    dvh_stats = {}
    
    ## record the volume
    dvh_vol = (dvh.volume, dvh.volume_units)
    dvh_stats['volume'] = (dvh_vol[0], dvh_vol[1])
    
    ## check if valid volume from the calcs
    if dvh_vol[0] == 0 or np.isnan(dvh_vol[0]):
        valid_vol = False
    else:
        valid_vol = True
    #print(dvh_vol,valid_vol)

    ## record the dosimetric statistics
    for stat in stat_list:
        if valid_vol == True:
            
            special_stats = {'mean':(dvh.mean, dvh.dose_units),
                'max':(dvh.max, dvh.dose_units),
                'min':(dvh.min, dvh.dose_units),
                'median':(dvh.statistic('D50').value,dvh.statistic('D50').units)}
            
            if stat.lower() in special_stats:
                dvh_stats[stat] = special_stats[stat.lower()]
                
            else:
                #dvh_stat = dvh.statistic(stat)
                dvh_stat = non_int_stat(dvh,stat)
                dvh_stats[stat] = (dvh_stat.value,dvh_stat.units)
        else:
            dvh_stats[stat] = None
    return dvh_stats


# In[65]:

def dicom_dvh_stats_single(directory=None, stats=None, structures=None, output='df', include_body=True,verbose=True,save_df=False,user_structures=None):
    """Get dvh statistics for the specified structures from the DICOM files in the given directory.
    directory = string :'C:\\.....'
    structures = list of structure IDs: [1,2,5....]. Default is None which will get stats from all structures.
    stats = list of stats: ['D90', 'V100'....]
    This will probably be wrapped by another funciton which will allow selection of the directory and the structure mapping.
    Raw results are as a dict, but default is a dataframe. Specifying anything other than 'df' will output a dict.
    output='df' must be used for processing multiple dicom sets"""
    
    if directory==None:
        directory = dir_open()
        
    ## list of default stats. HAs to be set as None in the function as called from multi patient version too
    default_stats = ['D98','D90','V100','V10', 'V20', 'D2cc', 'max', 'min', 'mean', 'median']
    if stats==None:
        stats = default_stats
    
    ## 1, 2, 3
    ## get the DICOM file paths into dictionary for use
    ## function also checks if the IDs match        
    all_exist_and_match, dir_files = all_rt_files_exist(directory)
    
    ## 4
    ## get structures which exist - dont need this if specifying structure ID
    ## get the structure names {name:id}
    structure_names = structure_dict(dir_files['rtss'])
    
    ## remove body structure if required (default)
    if include_body != True:
        if 'body' in structure_names:
            del structure_names['body']
    
    ## remove couch structures from analysis
    for the_struct in ['couchsurface','couchinterior']:
        if the_struct in structure_names:
            del structure_names[the_struct]
    
    ## remove any obvious pseudos:
    structure_names = {k:v for k,v in structure_names.items() if 'pseudo' not in k}
                      
    ## limit to only specified structure names by removing non-matching structures
    ## **********need to reconstruct dict by checking names are in list supplied by user... This below line is untested!******
    print('***')
    print(structure_names)
    if user_structures == 'mapping':
        the_structs = import_structure_list(file_open(title='Select File Containing Structure List...'))
        the_structs = [i.lower() for i in the_structs]
        structure_names =  {k:v for k,v in structure_names.items() if k.lower() in the_structs}
    elif user_structures is not None:
        ## if passed a list
        the_structs = [i.lower() for i in user_structures]
        structure_names =  {k:v for k,v in structure_names.items() if k.lower() in the_structs}
    print(structure_names)
    print('*********')
    
    ## swap the structure dict to give {id:name}
    all_structures = swap_dict(structure_names)
    #print (structure_names)
    
    ## if structure IDs not specified, then get data for all structures
    ## ideally user should be able specify the structure names rather than IDs.
    ## would need a funciton to get teh required IDs from the names - see above proposal
    if structures==None:
        structures = list(all_structures.keys())

    ## ******* I think this below chunk is now obsolete as list of structures is constructed differently.
#    if not all(structure in all_structures for structure in structures):
#        ## tell user which structures IDs do not match if there are any which do not match
#        ## then remove this from the list of structures to analyse
#        struct_id_mismatch = list_exist_in_dict(structures,all_structures)
#        print('The structure IDs', struct_id_mismatch, 'do not exist within the structure set.')
#        print(structures)
#        for struct_id in struct_id_mismatch:
#            structures.remove(struct_id)
    ## *********
    
    ## get the dvh for each structure and extract the required statistics
    pt_id = get_pt_id(dir_files['rtplan'])
    plan_name = get_dicom_plan_info(dir_files['rtplan'])['plan_label']
    #print(plan_name)
    prescription = get_dicom_dose(dir_files['rtplan'])
    
    print('Patient id:',pt_id, '\nPlan name:',plan_name, '\nGetting DVH statistics', stats,)
    
    all_dvh_stats = {}
    
    struct_num = 1
    tot_structs = len(structures)
    for structure_id in structures:
        structure_name = all_structures[structure_id]
        ## display ful progress if required
        if verbose == True:
            print('Patient id: ' + str(pt_id) + ', Plan name: ' + str(plan_name) +
                  ', Processing structure:' + str(struct_num) + ' of ' + str(tot_structs) + ' (' + structure_name + ')')
        struct_dvh = calc_dvh(dir_files['rtss'],dir_files['rtdose'],structure_id,prescription)
        ## get required stats from DVH
        structure_dvh_stats = get_statistics(struct_dvh,stats)
        all_dvh_stats[structure_name]=structure_dvh_stats
        #print(structure_id,structure_name,structure_dvh_stats)
        struct_num += 1
    
    print('Patient Completed')
    
    dir_string_output = os.path.split(directory)[1]
    
    results_dict = {(dir_string_output,pt_id,plan_name): all_dvh_stats}
    
    if output == 'df':
        results_output = pt_stats_df(results_dict)
        if save_df == True:
            pickle_save_dir = os.path.split(directory)[0]+ '\\' +os.path.split(directory)[1] + '.pkl'
            print(pickle_save_dir)
            results_output.to_pickle(pickle_save_dir)
    else:
        results_output = results_dict
    
    return results_output


# In[66]:

def pt_stats_df(results_dict):
    """Formats the results (which are in a dict) into a more user friendly dataframe"""
    
    ## get the patient results into a df
    results_df = pd.DataFrame.from_dict(results_dict[list(results_dict.keys())[0]],orient='index')
    #print(results_df)

    ## rename the index as 'structure' and then re-index for consistency of data access (i.e. index is seperate)
    results_df.index.rename('structure',inplace=True)
    results_df.reset_index(inplace=True,drop=False)
    
    ## insert the patient id as a seperate column at the start of the df (i.d. = first (only) key in results dict)
    results_df.insert(0,'plan_name',list(results_dict.keys())[0][2])
    results_df.insert(0,'patient_id',list(results_dict.keys())[0][1])
    results_df.insert(0,'sub_dir',list(results_dict.keys())[0][0])
    
    ## split the statistic results into values and units for simpler analysis/formatting
    for heading in results_df.columns.values:
        if heading not in ['sub_dir','structure', 'patient_id', 'plan_name']:
            head_val = heading + '_val'
            head_unit = heading + '_unit'

            results_df[[head_val,head_unit]] = results_df[heading].apply(pd.Series)
            
    ## convert pt_id, plan_name, structure columns to lowercase and remove trailing whitespace.
    results_df['patient_id'] = results_df['patient_id'].str.strip().str.lower()
    results_df['plan_name'] = results_df['plan_name'].str.strip().str.lower()
    results_df['structure'] = results_df['structure'].str.strip().str.lower()
    
    ## set all headers to lowercase
    results_df.columns = results_df.columns.str.lower()
    
    return results_df


# In[67]:

def get_sub_dirs(parent_dir=None):
    """Get a list of sub-directories (only first level) from the provided (or usually chosen form prompt) directory"""
    ## can allow user to specify path in code if tehy want to, otherwise it will prompt
    if parent_dir == None:
        parent_dir = dir_open("Select Directory Containing Subdirectories to Analyse...")

    ## extract the info from the directory
    path,dirs,files = next(os.walk(parent_dir))
    
    ## combine the path and the dirs into the required strings
    sub_dirs = [path + '/' + dirs[i] for i in range(len(dirs))]
    
    return(path,sub_dirs)


# In[68]:

def dicom_dvh_stats_multi(parent_directory=None, stats='prompt', structures=None,
                           save_results=True, save_prompt=False, include_body=True, limit=None,
                           struct_labels=False, verbose=True,save_df=False,user_structures=None):
    """Function to allow entire sub directory of folders to be analysed and the results combined.
    The number of directories can be limited using e.g. limit=2.
    The results can be saved as a csv file using save_results=True
    Save_prompt will aim to allow user specified save location and file name in the future.
    If verbose=True, then will show structures being processed."""
    if parent_directory == None:
        path, sub_dirs = get_sub_dirs()
    else:
        path, sub_dirs = get_sub_dirs(parent_dir = parent_directory)
    
    ## can get list of pseudo structure labels from user to add to results.
    if struct_labels == True:
        struct_label_file = file_open(title='Select File Containing Structure Labels...')
    
    ## can limit number of directories checked. Default is no limit.
    if limit is not None:
        sub_dirs = sub_dirs[:limit]
    
    ## get the dch stats from the user if required. (Only not required if they set stats=None)
    if stats == 'import':
        stats = import_dvh_stats(file_open(title='Select File Containing DVH Stats to Return...'))
    if stats == 'prompt':
        stats = prompt_dvh_stat_list()
        
    if user_structures == 'mapping':
        user_struct_list = import_structure_list(file_open(title='Select File Containing Structure List...'))
    else:
        user_struct_list = None
    
    ## get the results from each directory and append the results into a list which then becomes a df
    all_data_list = []
    num_sub_dirs = len(sub_dirs)
    sub_dir_count = 1
    for sub_dir in sub_dirs:
        print('Processing Directory ' + str(sub_dir_count) + ' of ' + str(num_sub_dirs) + ', ' + sub_dir)
        sub_dir_results = dicom_dvh_stats_single(directory=sub_dir, stats=stats, structures=structures, output='df',
                                           include_body=include_body, verbose=verbose,save_df=save_df,user_structures=user_struct_list)

        all_data_list.append(sub_dir_results)
        sub_dir_count += 1
    all_results = pd.concat(all_data_list,axis=0)
    print('Dosimetric data extracted from all directories')
    if struct_labels == True:
        all_results = add_struct_label_plan_id(all_results, struct_label_file)
        print('Adding struct_label structure names')
    else:
        all_results['struct_label'] = np.nan
        
    ## move the struct_label col to a better location in the df
    cols = all_results.columns.tolist() ## get a list of the columns
    cols.insert(3, cols.pop(cols.index('struct_label'))) ## get the required column (pop) and move it to new location
    all_results = all_results.reindex(columns= cols) ## reindex the df with new column order
    
    ## save the results
    if save_results == True:
        save_df_to_csv(all_results,prompt=save_prompt,save_dir=path)
    
    ## df of results is returned
    return all_results

    
#%%

def prompt_dvh_stat_list(prompt="Enter DVH Stats Seperate by a comma. e.g. D90,V100,D2cc: "):
    """ Get the stats input by the user and output as a list for use """
    print('\a')
    return input(prompt).replace(' ', '').split(',')

# In[69]:

def import_dvh_stats(file):
    """Import the DVH stats to calculate and return a list of these.
    The stats should be on seperate rows in the file, with no header."""
    
    ## get file extension to determine import method
    filetype = os.path.splitext(file)[1]
    #print(filetype)
    ## read in the struct_label structure file and ensure patient IDs are strings to match the DICOM data
    if filetype in ['.xlsx','.xls','.xlm','.xlsm']:
        df = pd.read_excel(file, index_col=None, header=None)
    elif filetype in ['.csv']:
        df = pd.read_csv(file, index_col=False, header=None)
    elif filetype in ['.txt']:
        df = pd.read_table(file, index_col=False, header=None)
    else:
        print('Can currently only import .csv, .txt, or excel filetypes for DVH stats')
    
    ## convert df to list
    dvh_stat_list = list(df[0].values)
    
    return dvh_stat_list
# %%
def import_structure_list(file):
    """Import the structures as a list.
    The stats should be on seperate rows in the file, with no header."""
    
    ## get file extension to determine import method
    filetype = os.path.splitext(file)[1]
    #print(filetype)
    ## read in the struct_label structure file and ensure patient IDs are strings to match the DICOM data
    if filetype in ['.xlsx','.xls','.xlm','.xlsm']:
        df = pd.read_excel(file, index_col=None, header=None)
    elif filetype in ['.csv']:
        df = pd.read_csv(file, index_col=False, header=None)
    elif filetype in ['.txt']:
        df = pd.read_table(file, index_col=False, header=None)
    else:
        print('Can currently only import .csv, .txt, or excel filetypes for structure names')
    
    ## convert df to list
    structure_list = list(df[0].values)
    
    return structure_list


# In[70]:

def import_struct_label_file(file):
    """Import file and get into correctly formatted dataframe for use.
    Supported filetypes are: excel, csv, txt"""
    
    ## get file extension to determine import method
    filetype = os.path.splitext(file)[1]
    #print(filetype)
    ## read in the struct_label structure file and ensure patient IDs are strings to match the DICOM data
    if filetype in ['.xlsx','.xls','.xlm','.xlsm']:
        df = pd.read_excel(file, index_col=None)
    elif filetype in ['.csv']:
        df = pd.read_csv(file, index_col=False)
    elif filetype in ['.txt']:
        df = pd.read_table(file, index_col=False)
    else:
        print('Can currently only import .csv, .txt, or excel filetypes for structure labels')
    
    ## do some tidying - lower case and remove whitespace
    df['patient_id'] = df['patient_id'].astype(str).str.strip().str.lower()
    df['structure'] = df['structure'].str.strip().str.lower()
    
    ## plan_name might not exist if all unique patients
    if 'plan_name' in df.columns:
        df['plan_name'] = df['plan_name'].str.strip().str.lower()
    
    return df


# In[71]:

def add_struct_label_plan_id(results_df,struct_label_file):
    """Function to check the struct_label structure file has the plan ID included.
    If it does not, then this is inferred from the patient ID from the results.
    Only if unique plan ids for each pateint can the struct_label strucutre names be reliably used.
    The struct_label_df is returned including the plan names"""
    
    ## read in the struct_label structure file and ensure patient IDs are strings to match the DICOM data
    df_structs = import_struct_label_file(struct_label_file)
    unique_pts = results_df['patient_id'].unique() ## get pateint ids from results
    df_structs = df_structs[df_structs['patient_id'].isin(unique_pts)] ## remove structure mapping info from uneeded pts
    

    ## if plan names do not exist, then add them from the results if they are unique for each patient
    if 'plan_name' not in df_structs.columns:
        print('Determing plan names from results df')
        #unique_pts = results_df['patient_id'].unique()
        #df_structs = df_structs[df_structs['patient_id'] in unique_pts] ## filter out structure mapping from pts not in results
        ## need to actually remove the rows completely so the lengths match properly...??

        ## store the plan names in a dict {pt_id, plan_name} for each patient.
        plan_names_dict = {}

        for patient in unique_pts:
            pt_plan_names = results_df[results_df['patient_id']==patient]['plan_name'].unique()
            plan_names_dict[patient] = pt_plan_names

        ## check length of all items in the dict. Should be 1 if unique plan names. Warn user if not.
        unique_plan_names = True
        for item in plan_names_dict:
            if len(plan_names_dict[item]) != 1 :
                print('Patient', item,
                      'does not have unique plan names. These must be specified within the struct_label structures file.')
                unique_plan_names = False

        ## if all unique then add to the dataframe
        if unique_plan_names == True:
            all_struct_pts = df_structs['patient_id'].values
            ## create list of pateint plan names and then add to df
            plan_names_to_add = []
            for patient in all_struct_pts:
                plan_names_to_add.append(plan_names_dict[patient][0])
            df_structs['plan_name'] = plan_names_to_add
            if len(df_structs) != 0: ## if nothing matches then cant set values whcih dont exist
                df_structs['plan_name'] = df_structs['plan_name'].str.strip().str.lower() ## lower case and strip whitespace
                print('Unique plan names added to struct_label structure data for use')
            else:
                df_structs['plan_name'] = 'No matching plan names'
                print('No matching plan names to add from structure mapping file')
            
    ## check psudo_struct df only contains unique rows to ensure no abiguity in results
    if df_structs.equals(df_structs.drop_duplicates()) == False:
        print('struct_label Structure Dataframe contains duplicates. Results may be ambiguous.')
    
    ## add the struct_label names to the results using the merge functionality.
    merged_df = results_df.merge(df_structs, how='outer', on=['patient_id', 'plan_name', 'structure'])
    merged_df = merged_df[~merged_df['sub_dir'].isnull()] ## remove any blank values which have appeared due to merge?
    
    return merged_df


# In[72]:

#my_results = dicom_dvh_stats_multi(save_results=True, save_prompt=True, stats='import',
#                                    limit=2 ,struct_labels=True, verbose=True, include_body=False)


# In[74]:

#my_results.head()


# In[59]:

#simple_results = dicom_dvh_stats_multi()


# ## Analysis of Results.
# - Make use of pandas to demonstrate simplicity of extracting data.

# In[75]:

## produce a boxplot of each measured dvh statistic

## produce a list of stats to plot (only want the values, and not volume)
#cols = list(my_results.columns)
## keep only items with '_val' at the end
#cols = [x for x in cols if '_val' in x]
## remove any other unnesessary thigns from the list
#to_remove = ['volume_val']
#cols = [x for x in cols if x not in to_remove]
#print(cols)

#my_results[cols].plot.box(figsize=(10,4))
#plt.show()


# ## Get non-integer stats
# e.g D0.1cc
# - probably best to check the string passed and see if there is a decimal.
#     - Then shoudl eb abel to determine method to use: https://groups.google.com/d/msg/dicompyler/EMnyhcEg4_Y/4P1wIcJ3AQAJ

# In[265]:

#plan = r'C:\Users\mb22\OneDrive\PhD\Quasar Shared\Data\Trials\PARSPORT\Tidied_Dicom_Files\1001\rtplan.dcm'
#ss = r'C:\Users\mb22\OneDrive\PhD\Quasar Shared\Data\Trials\PARSPORT\Tidied_Dicom_Files\1001\rtss.dcm'
#dose = r'C:\Users\mb22\OneDrive\PhD\Quasar Shared\Data\Trials\PARSPORT\Tidied_Dicom_Files\1001\rtdose.dcm'


# In[268]:

#abc = calc_dvh(ss,dose,3,get_dicom_dose(plan))


# In[351]:

## D90 (D90Gy is odd...? Does it even mean anything? Have I ever used it...?)
#print(abc.D9)
#print(abc.statistic('D9'))
#print(abc.dose_constraint(9.6))
#print('-')

#print(abc.D2cc)
#print(abc.statistic('D2cc'))
#print(abc.dose_constraint(2.2, volume_units='cc'))
#print('-')

#print(abc.V50)
#print(abc.statistic('V50'))
#print(abc.volume_constraint(50.3))
#print('-')

#print(abc.V50Gy)
#print(abc.statistic('V50Gy'))
#print(abc.volume_constraint(50, dose_units='Gy'))
#print('-')


# In[ ]:

#Dx = dose_constraint(x)
#Dxcc = dose_constraint(x, volume_units='cc')
#Vx = volume_constraint(x)
#VxGy = volume_constraint(x, dose_units='Gy')


# In[447]:

#abc.statistic('V50Gy')


# In[ ]:



