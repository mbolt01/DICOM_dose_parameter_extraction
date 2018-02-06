# -*- coding: utf-8 -*-
"""
Created on Thu Jun 15 14:51:57 2017

@author: mb22
"""

import dvh_analysis
#%%
import importlib
importlib.reload(dvh_analysis)

#%%

my_res = dvh_analysis.dicom_dvh_stats_multi(stats='import',include_body=False,
                                            struct_labels=True,save_df=True,
                                            user_structures=None)

