import ConfigParser
from collections import defaultdict
import json


def set_default_misowrap_settings(settings_info):
    """
    Set default misowrap settings.
    """
    # Set default MISO binary dir to empty (so that 
    # by default we look for MISO scripts already in path)
    settings_info["settings"]["miso_bin_dir"] = ""
    settings_info["settings"]["overhanglen"] = 1
    return settings_info


def load_misowrap_settings(config_filename,
                           # Integer parameters
                           INT_PARAMS=["readlen",
                                       "overhanglen"],
                           # Boolean parameters
                           BOOL_PARAMS=["paired"],
                           # Parameters to be interpreted as Python lists or
                           # data structures,
                           STR_PARAMS=["indir",
                                       "outdir",
                                       "stranded",
                                       "mapper"],
                           DATA_PARAMS=["sequence_files",
                                        "sample_groups"]):
    config = ConfigParser.ConfigParser()
    print "Loading settings from: %s" %(config_filename)
    parsed_settings = config.read(config_filename)
    settings_info = defaultdict(dict)
    for section in config.sections():
        for option in config.options(section):
            if option in FLOAT_PARAMS:
                settings_info[section][option] = config.getfloat(section, option)
            elif option in INT_PARAMS:
                settings_info[section][option] = config.getint(section, option)
            elif option in BOOL_PARAMS:
                settings_info[section][option] = config.getboolean(section, option)
            elif option in STR_PARAMS:
                settings_info[section][option] = str(config.get(section, option))
            elif option in DATA_PARAMS:
                settings_info[section][option] = json.loads(config.get(section, option))
            else:
                settings_info[section][option] = config.get(section, option)
    # Error-check the existing settings
    default_settings.check_settings(settings_info)
    # Set default values for settings 
    settings_info = set_default_misowrap_settings(settings_info)
    return settings_info, parsed_settings
