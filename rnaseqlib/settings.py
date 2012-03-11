import ConfigParser
from collections import defaultdict

def tryEval(s):
  try:
    return eval(s, {}, {})
  except:
    return s

def evalDict(d):
    for k, v in d.iteritems():
	d[k] = tryEval(v)
    return d

def load_settings(config_filename,
                  # Float parameters
                  FLOAT_PARAMS=[],
                  # Integer parameters
                  INT_PARAMS=["num_processors"],
                  # Boolean parameters
                  BOOL_PARAMS=["paired_end"],
                  # Parameters to be interpreted as Python lists or
                  # data structures
                  DATA_PARAMS=["sequence_files", "indir", "outdir",
                               "sample_groups", "stranded"]):
    config = ConfigParser.ConfigParser()

    print "Loading settings from: %s" %(config_filename)
    parsed_settings = config.read(config_filename)

    settings_info = defaultdict(dict)
    
    for section in config.sections():
        for option in config.options(section):
            if option in FLOAT_PARAMS:
                settings[section][option] = config.getfloat(section, option)
            elif option in INT_PARAMS:
                settings[section][option] = config.getint(section, option)
            elif option in BOOL_PARAMS:
                settings[section][option] = config.getboolean(section, option)
            elif option in DATA_PARAMS:
                settings[section][option] = json.loads(config.get(section, option))
            else:
                settings[section][option] = config.get(section, option)
    
    return settings_info, parsed_settings


def parse_plot_settings(settings_filename, event=None, chrom=None,
                        # Float parameters
                        no_posteriors=False):
    """
    Populate a settings dictionary with the plotting parameters, parsed
    as the right datatype.
    """
    settings = get_default_settings()
    
    config = ConfigParser.ConfigParser()

    print "Reading settings from: %s" %(settings_filename)
    config.read(settings_filename)
    
        
