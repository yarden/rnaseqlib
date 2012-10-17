##
## Setting of default settings for pipeline
##
import os
import sys
import time


def set_default_rnaseq_settings(settings_info):
    """
    Default settings that are RNA-Seq specific.
    """
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


def section_error(section)
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
            settings_error(section)


def set_default_settings(settings_info):
    """
    Set default global settings and data type
    specific settings.
    """
    check_settings()
    data_type = settings_info["pipeline"]["data_type"]
    if data_type == "rnaseq":
        settings_info = set_default_rnaseq_settings(settings)
    elif data_type == "riboseq":
        settings_info = set_default_riboseq_settings(settings)
    elif data_type == "clip":
        raise Exception, "Not implemented."
    elif data_type == "selex":
        raise Exception, "Not implemented."
    else:
        print "Unknown data type %s. Are you on crack?" \
            %(data_type)
    return settings_info
