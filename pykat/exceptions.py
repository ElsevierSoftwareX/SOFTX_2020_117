import exceptions

class BasePyKatException(Exception):
    def __init__(self, msg):
        self.__msg = msg
        
    def __str__(self):
        return self.__msg

class FinesseParse(BasePyKatException) :    
    def __init__(self, msg):
        BasePyKatException.__init__(self, "Error parsing Finesse input\n{0}".format(msg))
    
class MissingFinesseEnvVar(BasePyKatException) :    
    def __init__(self):
        BasePyKatException.__init__(self, "The environment variable FINESSE_DIR was not defined")

class MissingFinesse(BasePyKatException) :    
    def __init__(self):
        BasePyKatException.__init__(self, "Could not find the finesse executable 'kat' in '{0}'," \
                                     "or you do not have the permissions to run it." \
                                      .format(os.environ.get('FINESSE_DIR')))
    
class FinesseRunError(BasePyKatException) :
    def __init__(self, err, kat):
        self.__err = err
        self.__kat = kat
        
        BasePyKatException.__init__(self, "Finesse error (file {1}):\n{0}".format(self.__err, self.__kat))
        #BasePyKatException.__init__(self, "{0}".format(self.__err))
        
