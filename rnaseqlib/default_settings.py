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
    if (settings_info["mapping"]["paired"]) and \
       "paired_end_frag":
        # If paired-end, choose a default paired-end fragment length
        # if there isn't one already
        paired_end_frag = 300
        settings_info = set_settings_value(settings_info,
                                           "mapping",
                                           "paired_end_frag",
                                           paired_end_frag)
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


def check_settings(settings_info):
    """
    Error-check the settings.
    """
    required_sections = ["pipeline", "pipeline-files", "mapping", "data"]
    # Check that the major sections are in place
    for section in required_sections:
        if section not in settings_info:
            section_error(section)


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
