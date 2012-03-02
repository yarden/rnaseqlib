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

def load_settings(config_filename):
    config = ConfigParser.ConfigParser()

    print "Loading settings from: %s" %(config_filename)
    parsed_settings = config.read(config_filename)

    settings_info = defaultdict(dict)

    for section in config.sections():
        for option in config.options(section):
            settings_info[section][option] = tryEval(config.get(section, option))
    
    return settings_info, parsed_settings

