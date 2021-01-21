#!/usr/bin/env python

"""
Script that allows to run selected tasks as parallel subprocesses. Useful for
using multiple GPUs to run each one a different task.
"""

import libcommon
from importlib import import_module

# Loads the information about the refinement job and gets the tasks to run.
job, module, task_s = libcommon.setup_external_mode()
# Select the module to use.
selected_module = import_module(module)
# Actually runs the tasks.
selected_module.run_tasks(job, task_s)
