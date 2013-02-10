# -*- coding: utf-8 -*-
"""
Created on Sun Jan 27 10:02:41 2013

@author: Daniel
"""
import exceptions
from pykat.components import Component
from pykat.detectors import Detector

class NodeNetwork(object):
    def __init__(self, kat):
        self._nodes = {}
        self.__kat = kat
        
    def createNode(self, node_name):
        if node_name == 'dump':
            return DumpNode()
            
        if node_name in self._nodes:
            n = self._nodes[node_name]
            
            return n
        else:
            n = Node(node_name)
            self._nodes[node_name] = n
            return n
            
    def hasNode(self, name):
        return (name in self._nodes)
      
    def dumpInfo(self):
        
        for name in self._nodes:
            
            n = self._nodes[name]
            items = n.getComponents()
            comp = items[0][:]
            det = items[1]
            
            if comp[0] == None:
                comp1 = 'dump'
            else:
                comp1 = comp[0].name
            
            if comp[1] == None:
                comp2 = 'dump'
            else:
                comp2 = comp[1].name    
            
            detectors = ""
            
            if len(det) > 0:
                detectors = "Detectors: "
                
                for d in det:
                    detectors = detectors + d.name + " "
                
            print "node: {0} connected:{1} {2}->{3} {4}".format(
                    n.name,n.isConnected(),comp1, comp2, detectors)
        
class Node(object):
    
    def __init__(self, name):
        self._comp1 = None
        self._comp2 = None
        self._detectors = []
        self.__name = name
    
    def isConnected(self):
        return (self._comp1!=None) and (self._comp2!=None)
    
    def connect(self, obj):

        if not (isinstance(obj,Component) or isinstance(obj,Detector)):
            raise exceptions.ValueError("Object is not a component or detector")
        
        if isinstance(obj, Component):
            
            if self._comp1 == None:
                self._comp1 = obj
            elif self._comp2 == None:
                self._comp2 = obj
            else:
                raise exceptions.RuntimeError("Cannot attach {0} to node {1} as it is already connected: {2} and {3}".format(
                                    obj.name, self.__name,self._comp1.name,self._comp2.name))            
        else:
            # we must have a detector as we check above            
            self._detectors.append(obj)
    
    def getComponents(self):
        ''' Returns a tuple with the first being the 2 components that 
        connect to the node. The second item is a list of the detectors at the
        node'''
        return [(self._comp1, self._comp2),self._detectors]
        
    def amIConnected(self, obj):
        """
        Checks if obj is connected oto the node. Returns true or false in tuple
        with None or the other object and the node index which it is attached to
        """ 
        if obj == self._comp1:
            if self._comp2 == None:
                ix = -1
            else:
                ix = self._comp2.getNodes().index(self)
                
            return [True, self._comp2, ix]
        elif obj == self._comp2:
            if self._comp1 == None:
                ix = -1
            else:
                ix = self._comp1.getNodes().index(self)
                
            return [True, self._comp1, ix]
        else:
            return [False, None]
        
    def __getname(self):
        return self.__name      
        
    name = property(__getname)
    
class DumpNode(Node):
    def __init__(self):
        Node.__init__(self,'dump')
        
        