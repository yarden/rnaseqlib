import ConfigParser
from collections import defaultdict
import json

# def tryEval(s):
#   try:
#     return eval(s, {}, {})
#   except:
#     return s

# def evalDict(d):
#     for k, v in d.iteritems():
# 	d[k] = tryEval(v)
#     return d

def load_settings(config_filename,
                  # Float parameters
                  FLOAT_PARAMS=[],
                  # Integer parameters
                  INT_PARAMS=["num_processors"],
                  # Boolean parameters
                  BOOL_PARAMS=["paired",
                               "compressed"],
                  # Parameters to be interpreted as Python lists or
                  # data structures,
                  STR_PARAMS=["indir",
                              "outdir",
                              "stranded"],
                  DATA_PARAMS=["sequence_files", 
                               "sample_groups"]):
    config = ConfigParser.ConfigParser()
    print "Loading settings from: %s" %(config_filename)
    parsed_settings = config.read(config_filename)
    settings_info = defaultdict(dict)
    
    for section in config.sections():
        print "Parsing section: ", section
        for option in config.options(section):
            print "option -> ", option
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
    return settings_info, parsed_settings

        
