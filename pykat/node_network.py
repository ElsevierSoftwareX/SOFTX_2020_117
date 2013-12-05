# -*- coding: utf-8 -*-
"""
Created on Sun Jan 27 10:02:41 2013

@author: Daniel
"""
import exceptions
import pykat.gui.graphics
import pykat.exceptions as pkex
from pykat.components import Component
from pykat.detectors import Detector

class NodeNetwork(object):
    def __init__(self, kat):
        self.__nodes = {}
        self.__kat = kat
        self.__nodeComponents = {} # dictionary of tuples containing which components are connected to a node
        self.__componentNodes = {} # dictionary of tuples containing which nodes are connected to a given component
        self.__componentCallback = {}
        self.__node_id = 1
        
    def registerComponentNodes(self, comp, node_names, change_callback):
        """
        For a given component we create some nodes or get existing ones and 
        attach them to this component. Also specify a callback function that
        is called whenever the nodes attached to this component are changed
        , e.g. connected, disconnected, name change, etc.
        """
        if not isinstance(comp, Component):
            raise exceptions.ValueError("comp argument is not of type Component")
        
        if comp.id in self.__componentNodes:
            raise pkex.BasePyKatException("Component has already been registered")
        
        list = []
        
        for name in node_names:
            n = self.createNode(name)
            self.connectNodeToComp(n, comp, do_callback=False)
            list.append(n)
        
        self.__componentNodes[comp.id] = tuple(list)
        self.__componentCallback[comp.id] = change_callback
        
        change_callback()
    
    def connectNodeToComp(self, node, comp, do_callback=True):
        if node.id in self.__nodeComponents:
            comps = self.__nodeComponents[node.id]
        else:
            comps = ()
        
        if len(comps) >= 2:
            raise pkex.BasePyKatException("Node is already connected to 2 components")
        
        l = list(comps)
        l.append(comp)
        
        self.__nodeComponents[node.id] = tuple(l)
        
        if do_callback: self.__componentCallback[comp.id]()
        
    def createNode(self, node_name):
        if node_name == 'dump':
            return DumpNode()
            
        if node_name in self.__nodes:
            # then this node already exists
            return self.__nodes[node_name]
        else:
            n = Node(node_name, self, self.__node_id)
            
            self.__node_id += 1
            self.__add_node_attr(n) # add node as a member of this object, e.g. kat.nodes.n
            self.__nodes[node_name] = n
            return n
        
    def removeNode(self, node):
        if not isinstance(node,Node):
            raise exceptions.ValueError("node argument is not of type Node")
        
        if node.name not in self.__nodes:
            raise exceptions.RuntimeError("Trying to remove node {0} when it has not been added".format(node.name))
        
        C = self.getNodeComponents(node)
        
        if C[0] is not None or C[1] is not None:
            raise exceptions.RuntimeError("Cannot remove a node which is attached to components")
            
        if len(node.getDetectors()) > 0:
            raise exceptions.RuntimeError("Cannot remove a node which is attached to detectors still")
        
        self.__remove_node_attr(node)
        del self.__nodes[node.name] 
        
    def hasNode(self, name):
        return (name in self.__nodes)
    
    def getNodes(self):
        return self.__nodes.copy()
    
    def dumpInfo(self):
        
        for name in self.__nodes:
            
            n = self.__nodes[name]
            
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
    
    def getComponentNodes(self, comp): return self.__componentNodes[comp.id]
    def getNodeComponents(self, node): return self.__nodeComponents[node.id]
    
    def __add_node_attr(self, node):

        if not isinstance(node, Node):
            raise exceptions.ValueError("Argument is not of type Node")
        
        name = node.name
        fget = lambda self: self.__get_node_attr(name)
        
        setattr(self, name, property(fget))
        setattr(self, '__node_' + name, node)                   
    
    def __remove_node_attr(self, node):
        if not isinstance(node, Node):
            raise exceptions.ValueError("Argument is not of type Node")
        
        name = node.name
        detattr(self, '__node_' + name)
        delattr(self, name)
        
    def __get_node_attr(self, name):
        return getattr(self, '__node_' + name)        
        
class Node(object):

    def __init__(self, name, network, id):
        self._detectors = []
        self.__name = name
        self._item = None
        self._network = network
        self.__q_x = None
        self.__q_y = None
        self.__q_comp = None
        self.__id = id
        
    @property
    def id(self): return self.__id
    
    @property
    def network(self): return self._network
    
    @property
    def components(self): return self._network.getNodeComponents(self)
    
    @property
    def q(self):
        if self.__q_x == self.__q_y:
            return self.__q_x
        else:
            return (self.__q_x, self.__q_y)
            
    @property
    def qx(self): return self.__q_x
    @property
    def qy(self): return self.__q_y
    
    def removeGauss():
        self.__q_x = None
        self.__q_y = None
        self.__q_comp = None
    
    def setGauss(self, component, *args):
        self.__q_comp = component
        
        if len(args) == 1:
            self.__q_x = args[0]
            self.__q_y = args[0]
        elif len(args) == 2:
            self.__q_x = args[0]
            self.__q_y = args[1]
        else:
            raise pkex.BasePyKatException("Must specify either 1 Gaussian beam parameter or 2 for astigmatic beams")
        
    def getFinesseText(self):    
        if self.__q_x is None or self.__q_y is None or self.__q_comp is None:
            return []
            
        rtn = []
        
        if self.__q_x == self.__q_y:
            rtn.append("gauss* g_{node} {comp} {node} {z} {zr}".format(node=self.name, comp=self.__q_comp.name, z=self.__q_x.real, zr=self.__q_x.imag))
        else:
            rtn.append("gauss* g_{node} {comp} {node} {zx} {zrx} {zy} {zry}".format(node=self.name, comp=self.__q_comp.name, zx=self.__q_x.real, zrx=self.__q_x.imag, zy=self.__q_y.real, zry=self.__q_y.imag))
            
        return rtn
        
    def isConnected(self):
        if (self.components[0] is not None) and (self.self.components[1] is not None):
            return True
        else:
            return False
      
    def remove(self):
        self._network.removeNode(self)
        
        if self._item != None:
            self._item.scene().removeItem(self._item)
    
    def getQGraphicsItem(self,dx=0,dy=0,nsize=8,parent=None):
        if self._item == None:
            self._item = pykat.gui.graphics.NodeQGraphicItem(self,
                                                             dx,dy,
                                                             -nsize/2,-nsize/2,
                                                             nsize, nsize, parent)
            
        return self._item
    
    def getDetectors(self):
        return self._detectors[:]
        
    def amIConnected(self, obj):
        """
        Checks if obj is connected to the node. Returns true or false in tuple
        with None or the other object and the node index which it is attached to
        """ 
        comps = self.components
        
        if obj == comps[0]:
            if comps[1] == None:
                ix = -1
            else:
                ix = comps[1].getNodes().index(self)
                
            return [True, comps[1], ix]
            
        elif obj == comps[1]:
            if comps[0] == None:
                ix = -1
            else:
                ix = comps[0].getNodes().index(self)
                
            return [True, comps[0], ix]
        else:
            return [False, None]
        
    @property
    def name(self): return self.__name      
    
    
class DumpNode(Node):
    def __init__(self):
        Node.__init__(self, 'dump', None, -1)
        
        