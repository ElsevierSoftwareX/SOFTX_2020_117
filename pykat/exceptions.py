import exceptions

class BasePyKatException:
    def __init__(self, msg):
        self.__msg = msg
        
    def __str__(self):
        return "PyKat Exception message: ", self.__msg

class MissingFinesseEnvVar(BasePyKatException) :    
    def __init__(self):
        BasePyKatExeception.__init__("The environment variable FINESSE_DIR was not defined")

class MissingFinesse(BasePyKatException) :    
    def __init__(self):
        BasePyKatExeception.__init__("Could not find the finesse executable 'kat' in '{0}'," \
                                     "or you do not have the permissions to run it." \
                                      .format(os.environ.get('FINESSE_DIR')))
    
class FinesseRunError(BasePyKatException) :
    def __init__(self, err, kat):
        self.__err = err
        self.__kat = kat
        
        BasePyKatExeception.__init__("Finesse returned an error running {1}: {0}".format(self.__err, self.__kat))
        