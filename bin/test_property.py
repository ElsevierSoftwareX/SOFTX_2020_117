class Param(float):
    def __new__(self,name,value):
        return float.__new__(self,value)
         
    def __init__(self,name,value):
        self.__name = name
        
    name = property(lambda self: self.__name)
    
    
class Beer(object):
    def __init__(self, temp):
        self.__T = temp
        
    @property
    def temp(self):
        return Param('Beer Temperature', self.__T)
    
    @temp.setter
    def temp(self,value):
        self.__T = float(value) 
        

b = Beer(100)

print b.temp.name, b.temp.__class__

print b.temp * 4.5

print b.temp.name, b.temp.__class__

b.temp = 101

print b.temp.name, b.temp.__class__

print b.temp 