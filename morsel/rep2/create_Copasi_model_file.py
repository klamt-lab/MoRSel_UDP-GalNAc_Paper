#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys

sys.path.append("..")
from func_lib import *

# load parameter ensemble of the best model variant that was selected from an ImpExtSearch result
reading_file = open('ImpExtSearch_v23b_rep2_ModelVar5_50runs_evaluated_vars_log', 'rb')
selected_model_var = pickle.load(reading_file)
reading_file.close()

# export the model variant with the best median and mean rank as .cps file to use in the batch optimization
output_file_name = "ImpExtSearch_v23b_rep2_MV5.cps"
build_copasi_model_file_from_log_entry(selected_model_var[0], output_file_name)
