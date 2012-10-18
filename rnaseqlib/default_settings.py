##
## Setting of default settings for pipeline
##
import os
import sys
import time

def set_settings_value(settings_info, section,
                       param, value):
    """
    Set settings value if there isn't a value already.
    """
    if section not in settings_info:
        return settings_info
    # If parameter does not have a value already, set it
    if param not in settings_info[section]:
        settings_info[section][param] = value
    return settings_info


def set_default_rnaseq_settings(settings_info):
    """
    Default settings that are RNA-Seq specific.
    """
    if settings_info["mapping"]["paired"]:
        # Compute mate inner dist based on read length
        mate_inner_dist = \
            settings_info["mapping"]["paired_end_frag"] - (2 * settings_info["mapping"]["readlen"])
        settings_info = set_settings_value(settings_info,
                                           "mapping",
                                           "mate_inner_dist",
                                           mate_inner_dist)
    return settings_info


def set_default_riboseq_settings(settings_info):
    """
    Default settings that are Ribo-Seq specific.
    """
    return settings_info


def set_default_clip_settings(settings_info):
    """
    Default settings that are CLIP-Seq specific.
    """
    return settings_info


def section_error(section):
    print "Error in settings: cannot find %s section." \
        %(section)
    sys.exit(1)


def param_error(param):
    print "Error in settings: cannot find parameter %s." \
        %(param)
    sys.exit(1)


def check_settings(settings_info):
    """
    Error-check the settings.
    """
    required_sections = ["pipeline", "pipeline-files", "mapping", "data"]
    # Check that the major sections are in place
    for section in required_sections:
        if section not in settings_info:
            section_error(section)
    # Check that the major parameters are in place
    mapping_params = ["readlen"]
    for param in mapping_params:
        if param not in settings_info["mapping"]:
            param_error(param)
    # Check that paired-end specific parameters are correct
    if ("paired" in settings_info["mapping"]) and settings_info["mapping"]["paired"]:
        if "paired_end_frag" not in settings_info["mapping"]:
            print "Error: Need \'paired_end_frag\' to be set for paired " \
                  "samples. Is your data paired-end?"
            sys.exit(1)
            

def set_default_settings(settings_info):
    """
    Set default global settings and data type
    specific settings.
    """
    data_type = settings_info["pipeline"]["data_type"]
    if data_type == "rnaseq":
        settings_info = set_default_rnaseq_settings(settings_info)
    elif data_type == "riboseq":
        settings_info = set_default_riboseq_settings(settings_info)
    elif data_type == "clip":
        raise Exception, "Not implemented."
    elif data_type == "selex":
        raise Exception, "Not implemented."
    else:
        print "Unknown data type %s. Are you on crack?" \
            %(data_type)
    return settings_info
