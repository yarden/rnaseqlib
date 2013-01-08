import os
import sys
import time

import ConfigParser
from collections import defaultdict
import json


def set_default_misowrap_settings(settings_info,
                                  chunk_jobs=600):
    """
    Set default misowrap settings.
    """
    # Set default MISO binary dir to empty (so that 
    # by default we look for MISO scripts already in path)
    if "miso_bin_dir" not in settings_info["settings"]:
        settings_info["settings"]["miso_bin_dir"] = ""
    if "overhanglen" not in settings_info["settings"]:
        settings_info["settings"]["overhanglen"] = 1
    # If given a cluster type, then assume we run on cluster
    # and assign a default chunk_jobs parameter
    if "cluster_type" in settings_info["settings"]:
        if not "chunk_jobs" in settings_info["settings"]:
            settings_info["settings"]["chunk_jobs"] = chunk_jobs
    # If sample labels is not given, set default labels
    if "sample_labels" not in settings_info["data"]:
        # Set default sample labels
        settings_info["data"]["sample_labels"] = []
        for bam_label, bam_file in settings_info["data"]["bam_files"]:
            settings_info["data"]["sample_labels"].append([bam_label,
                                                           bam_label])
    # If no comparison groups are given, treat the entire
    # sample set as a group (i.e. compute all pairwise
    # comparisons between samples)
    if "comparison_groups" not in settings_info["data"]:
        sample_labels = [label[0] for label in \
                         settings_info["data"]["sample_labels"]]
        settings_info["data"]["comparison_groups"] = sample_labels
    return settings_info


def check_misowrap_settings(settings_info):
    """
    Check that we were given the proper settings.
    """
    print "Checking the validity of misowrap settings..."
    # Required parameters, indexed by section
    required_params = {"settings": ["readlen",
                                    "miso_events_dir",
                                    "miso_settings_filename",
                                    "miso_output_dir",
                                    "cluster_type"],
                       "pipeline-files": ["init_dir"],
                       "data": ["bam_files"]}
    for sect, params in required_params.iteritems():
        if sect not in settings_info:
            print "Error: Section %s not found in settings." %(sect)
            sys.exit(1)
        for param in params:
            if param not in settings_info[sect]:
                print "Error: %s is not set in section %s in " \
                      "settings file." %(param, sect)
                sys.exit(1)
    print "Settings appear valid."


def load_misowrap_settings(config_filename,
                           # Integer parameters
                           INT_PARAMS=["readlen",
                                       "overhanglen",
                                       "chunk_jobs",
                                       # Filters for events
                                       "atleast_inc",
                                       "atleast_exc",
                                       "atleast_sum",
                                       "atleast_const"],
                           BOOL_PARAMS = [],
                           # Parameters to be interpreted as Python lists or
                           # data structures,
                           STR_PARAMS=["cluster_type",
                                       "init_dir",
                                       "miso_bin_dir",
                                       "miso_events_dir",
                                       "miso_settings_filename",
                                       "miso_output_dir",
                                       "events_to_genes_dir",
                                       "insert_lens_dir"],
                           DATA_PARAMS=["bam_files",
                                        "sample_labels",
                                        "comparison_groups"]):
    config = ConfigParser.ConfigParser()
    print "Loading settings from: %s" %(config_filename)
    if not os.path.isfile(config_filename):
        print "Error: settings file %s not found." %(config_filename)
        sys.exit(1)
    parsed_settings = config.read(config_filename)
    settings_info = defaultdict(dict)
    for section in config.sections():
        for option in config.options(section):
            if option in INT_PARAMS:
                settings_info[section][option] = config.getint(section, option)
            elif option in BOOL_PARAMS:
                settings_info[section][option] = config.getboolean(section, option)
            elif option in STR_PARAMS:
                settings_info[section][option] = str(config.get(section, option))
            elif option in DATA_PARAMS:
                print "Loading %s" %(option)
                settings_info[section][option] = json.loads(config.get(section, option))
            else:
                settings_info[section][option] = config.get(section, option)
    # Error-check the existing settings
    check_misowrap_settings(settings_info)
    # Set default values for settings 
    settings_info = set_default_misowrap_settings(settings_info)
    return settings_info, parsed_settings
