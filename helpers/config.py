import ConfigParser as cfg
import yaml
import os

def Config(path):
    config = cfg.RawConfigParser(allow_no_value=True)
    config.readfp(open(path))
    parameters = {section: {key: value for key, value in config.items(section)} 
                                 for section in config.sections()}
    return parameters

class Config2(object):
    '''Flexible Config class that
    reads in dictionaries, cfg- and yaml-files.
    Saves in the yaml format.
    '''
    def _copy_data(self, data):
        '''copies all the data into a single
        dictionary.
        '''
        for key,val in data.items():
            self.parameters[key] = val
    
    def read(self, source=None):
        if not source:
            source = self.config_path
        if isinstance(source, basestring) and source.endswith('.yml'):
            self._read_yaml(source)
        elif isinstance(source, basestring) and source.endswith('.cfg'):
            self._read_cfg(source)
        elif isinstance(source, dict):
            self._copy_data(source)
    
    def _read_yaml(self, path):
        # read in parameters
        with open(path) as stream:
            config = yaml.load(stream)
        self._copy_data(config)
        
    def _read_cfg(self, path):
        config = cfg.RawConfigParser(allow_no_value=True)
        config.readfp(open(path))
        parameters = {section: {key: value for key, value in config.items(section)} 
                                 for section in config.sections()}
        self._copy_data(parameters)
        
    def save(self, path=None):
        '''Only saves to yaml format.
        '''
        if not path:
            root, ext = os.path.splitext(self.config_path)
            path = root + '_save.yml'
        with open(path, 'w') as f:
            yaml.dump(self.parameters, f)

