import os

class Plugin(object):
    IS_DAHU_PLUGIN = True
    DEFAULT_SET_UP = "setup"      # name of the method used to set-up the plugin (close connection, files)
    DEFAULT_PROCESS = "process"   # specify how to run the default processing
    DEFAULT_TEAR_DOWN = "teardown"# name of the method used to tear-down the plugin (close connection, files)
    def __init__(self):
        """         
        
        """
        self.init_param=None
        self.output = None
        self._logging = [] #

    def setup(self, kargs=None):
        """
        This is the second constructor to setup 
        input variables and possibly initialize
        some objects 
        """
        self.init_param = kargs


    def process(self, kargs=None):
        """
        main processing of the plugin
        """
        pass

    def teardown(self):
        """
        method used to tear-down the plugin (close connection, files)
        """
        pass

    def get_info(self):
        """
        """
        return os.path.linesep.join(self._logging)

if __name__ == "__main__":
    #here I should explain how to run the plugin as stand alone:
    p = Plugin()
    p.setup()
    p.process()
    p.teardown()
