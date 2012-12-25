import ConfigParser
from collections import defaultdict
import json

import default_settings

def load_settings(config_filename,
                  # Float parameters
                  FLOAT_PARAMS=[],
                  # Integer parameters
                  INT_PARAMS=["readlen",
                              "overhanglen",
                              "num_processors",
                              "paired_end_frag"],
                  # Boolean parameters
                  BOOL_PARAMS=["paired",
                               "prefilter_miso"],
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
                settings_info[section][option] = \
                    config.getfloat(section, option)
            elif option in INT_PARAMS:
                settings_info[section][option] = \
                    config.getint(section, option)
            elif option in BOOL_PARAMS:
                settings_info[section][option] = \
                    config.getboolean(section, option)
            elif option in STR_PARAMS:
                settings_info[section][option] = \
                    str(config.get(section, option))
            elif option in DATA_PARAMS:
                settings_info[section][option] = \
                    json.loads(config.get(section, option))
            else:
                settings_info[section][option] = \
                    config.get(section, option)
    # Error-check the existing settings
    default_settings.check_settings(settings_info)
    # Set default values for settings 
    settings_info = default_settings.set_default_settings(settings_info)
    return settings_info, parsed_settings
