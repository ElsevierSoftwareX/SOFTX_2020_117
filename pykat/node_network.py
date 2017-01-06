# -*- coding: utf-8 -*-
"""
Created on Sun Jan 27 10:02:41 2013

@author: Daniel
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import pykat
from pykat import USE_GUI, NoGUIException

if USE_GUI:
    import pykat.gui.graphics
    
import pykat.exceptions as pkex
import pykat.external.six as six

from pykat.components import Component, NodeGaussSetter
from pykat.detectors import BaseDetector as Detector
from pykat.optics.gaussian_beams import BeamParam
from copy import deepcopy

id___ = 0
    
class NodeNetwork(object):
    
    def __new__(cls, *args, **kwargs):
        # This creates an instance specific class for the component
        # this enables us to add properties to instances rather than
        # all classes
        global id___
        id___ += 1
        cnew_name = str("%s.%s_%i" % (cls.__module__, cls.__name__, id___))
        
        cnew = type(cnew_name, (cls,), {})
        
        #return object.__new__(cnew, *args, **kwargs)
        return object.__new__(cnew)

    def __init__(self, kat):
        self.__nodes = {}
        self.__kat = kat
        self.__nodeComponents = {} # dictionary of tuples containing which components are connected to a node
        self.__componentNodes = {} # dictionary of tuples containing which nodes are connected to a given component
        self.__componentCallback = {}
        self.__node_id = 1
    
    def __deepcopy__(self, memo):
        """
        When deep copying a kat object we need to take into account
        the instance specific properties.
        """
        
        # Here we create a copy of this object based of the base class
        # of this one, otherwise we're making a copy of a copy of a copy...
        result = self.__class__.__new__(self.__class__.__base__)
        result.__dict__ = deepcopy(self.__dict__, memo)
        
        return result
                
    @property
    def kat(self): return self.__kat
        
    def registerComponentNodes(self, comp, node_names, change_callback):
        """
        For a given component we create some nodes or get existing ones and 
        attach them to this component. Also specify a callback function that
        is called whenever the nodes attached to this component are changed
        , e.g. connected, disconnected, name change, etc.
        """
        if not isinstance(comp, Component):
            raise pkex.BasePyKatException("comp argument is not of type Component")
        
        if comp.id in self.__componentNodes:
            raise pkex.BasePyKatException("Component has already been registered")
        
        list = []
        
        for name in node_names:
            n = self.createNode(name)
            self.__connectNodeToComp(n, comp, do_callback=False)
            list.append(n)
        
        self.__componentNodes[comp.id] = tuple(list)
        self.__componentCallback[comp.id] = change_callback
        
        change_callback()
    
    def replaceNode(self, comp, node_old, node_new):
        """
        For a particular pykat component this will replace a node that is currently
        connected to it with another. This can be used to dynamically change layouts
        once components have been read into the pykat.finesse.kat object.
        
        node_old is the node that is attached to the component. This will accept
             str - name of a node
             pykat.node_network.Node - The direct node object
             NodeGaussSetter - the node object that is used to set gaussian parameters
             
        This will call a components __on_node_change callback function to let it know that the nodes
        connected to it have changed.
        """
        
        if isinstance(node_old, six.string_types):
            node_old = self.__kat.nodes[node_old]
        
        if isinstance(node_new, six.string_types):
            node_new = self.__kat.nodes[node_new]
            
        if isinstance(node_old, NodeGaussSetter):
            node_old = node_old.node
        
        if isinstance(node_new, NodeGaussSetter):
            node_new = node_new.node
            
        if not node_new.isDump and node_new.components.count(None) == 0:
            raise pkex.BasePyKatException("New node already connected to two components")
            
        if comp not in node_old.components:
            raise pkex.BasePyKatException("Old node not attached to component")
        
        if not node_new.isDump and comp in node_new.components:
            raise pkex.BasePyKatException("New node already attached to component")
        
        # add component to new node component list
        new_node_comps = list(node_new.components)
        new_node_comps[new_node_comps.index(None)] = comp
        self.__nodeComponents[node_new.id] = tuple(new_node_comps)
        del new_node_comps
        
        # remove component from old node list
        old_node_comps = list(node_old.components)
        old_node_comps[old_node_comps.index(comp)] = None
        self.__nodeComponents[node_old.id] = tuple(old_node_comps)
        del old_node_comps
        
        comp_nodes = list(comp.nodes)
        comp_nodes[comp_nodes.index(node_old)] = node_new
        self.__componentNodes[comp.id] = tuple(comp_nodes)
        del comp_nodes
        
        # if old node is no longer connected to anything then delete it
        if node_old.components.count(None) == 2:
            self.removeNode(node_old)
        
        # Call component callback to let it know that we have changed the 
        # nodes attached to it
        self.__componentCallback[comp.id]()
            
    def __connectNodeToComp(self, node, comp, do_callback=True):
        """
        This is an internal function used to create connections between nodes
        """
        if node.id in self.__nodeComponents:
            comps = self.__nodeComponents[node.id]
        else:
            comps = (None,) * 2
        
        if len(comps) >= 2 and comps[0] != None and comps[1] != None:
            raise pkex.BasePyKatException("Node '{0}' is already connected to 2 components ({1}, {2})".format(node.name, comps[0], comps[1]))
        
        l = list(comps)
        
        if l[0] is None:
            l[0] = comp
        elif l[1] is None:
            l[1] = comp
        else:
            raise pkex.BasePyKatException("Connected to two coponents already")
        
        self.__nodeComponents[node.id] = tuple(l)
        
        if do_callback: self.__componentCallback[comp.id]()
    
    def __update_nodes_properties(self):
        # check if any node setters have already been added. If so we
        # need to remove them. This function should get called if the nodes
        # are updated, either by some function call or the GUI
        key_rm = [k for k in self.__dict__ if k.startswith("__node_", 0, 7)]

        # now we have a list of which to remove
        for key in key_rm:
            ns = self.__dict__[key]
            name = str(ns.name)
            
            if '__node_' + name in self.__dict__:
                delattr(self, '__node_' + name)
            
            if name in self.__class__.__dict__:
                delattr(self.__class__, name)
        
        # Now re-add them pointing to the recent nodes
        for node in self.__nodes:
            if not self.__nodes[node].isDump:
                self.__add_node_attr(self.__nodes[node])
         
    def createNode(self, node_name):
        """
        This creates a new node object. It won't be connected to anything or added to a
        pykat.finesse.kat object until it is specifically attached to a particular 
        component. This should be used in conjunction with kat.nodes.replaceNode to 
        add a new node into a system, as every component will already have the nodes
        setup, including dump nodes.
        
        This will return a dump node if the name of the node is "dump" (case senstive)
        """
            
        if node_name != 'dump' and node_name in self.__nodes:
            # then this node already exists
            return self.__nodes[node_name]
        else:
            if node_name == 'dump':
                n = DumpNode(self)
            else:
                n = Node(node_name, self, self.__node_id)
            
            self.__node_id += 1
            self.__nodeComponents[n.id] = (None, None)
            
            if not n.isDump:
                self.__add_node_attr(n) # add node as a member of this object, e.g. kat.nodes.n
                self.__nodes[node_name] = n
                
            
            return n
    
    def _removeComponent(self, comp):
        """
        This is an internal function that shouldn't be used directly. This removes
        a particular component from the node network. For this to work it has to be 
        detached from all other connections first.
        """
        C = self.__componentNodes[comp.id]
        
        for n in C:
           if comp in self.__nodeComponents[n.id]:
               l = list(self.__nodeComponents[n.id])
               l[l.index(comp)] = None
               self.__nodeComponents[n.id] = tuple(l)
               
               if l.count(None) == 2:
                   self.removeNode(n) 
               
               del l
               
        del self.__componentCallback[comp.id]
        del self.__componentNodes[comp.id]
        
    def removeNode(self, node):
        """
        This will remove a particular node object from the network. The node in question
        must be fully detached from all components and connections first. This function is 
        called by replaceNode directly so a replaced node, that is no longer connected to 
        anything, is removed automatically.
        
        node_old is the node that is attached to the component. This will accept
             str - name of a node
             pykat.node_network.Node - The direct node object
             NodeGaussSetter - the node object that is used to set gaussian parameters
              
        """
        
        if isinstance(node, six.string_types):
            node = self.__kat.nodes[node]
            
        if isinstance(node, NodeGaussSetter):
            node = node.node
            
        if not isinstance(node, Node):
            raise pkex.BasePyKatException("node argument is not of type Node")
        
        if not isinstance(node, DumpNode) and node.name not in self.__nodes:
            raise pkex.BasePyKatException("Trying to remove node {0} when it has not been added".format(node.name))
        
        C = self.getNodeComponents(node)
        
        if C[0] is not None or C[1] is not None:
            raise pkex.BasePyKatException("Cannot remove a node which is attached to components still")
            
        if len(node.getDetectors()) > 0:
            raise pkex.BasePyKatException("Cannot remove a node which is attached to detectors still")
        
        if not isinstance(node, DumpNode):
            self.__remove_node_attr(node)
            del self.__nodes[node.name] 
            
        del self.__nodeComponents[node.id]
        
    def hasNode(self, name):
        ""
        return (name in self.__nodes)
    
    def getNodes(self):
        """
        Returns a copy of the node dictionary, this is for infomration purposes any edits won't make
        any changes to the node network.
        """
        return self.__nodes.copy()
    
    def getComponentNodes(self, comp):
        """
        This function returns a tuple of the nodes connected to the component specified.
        For information only, you cannot edit the connections using this function.
        """
        return self.__componentNodes[comp.id]
    
    def getNodeComponents(self, node):
        """
        This function returns a tuple of the components connected to the node specified.
        For information only, you cannot edit the connections using this function.
        """
        return self.__nodeComponents[node.id]
    
    def __add_node_attr(self, node):

        if not isinstance(node, Node):
            raise pkex.BasePyKatException("Argument is not of type Node")
        
        name = node.name
        fget = lambda self: self.__get_node_attr(name)
        
        setattr(self.__class__, name, property(fget))
        setattr(self, '__node_' + name, node)                   
    
    def __remove_node_attr(self, node):
        if not isinstance(node, Node):
            kat.nodes.replaceNode(kat.bs1, "n1", kat.nodes.createNode("test1"))
        
        name = node.name
        
        delattr(self, '__node_' + name)
        delattr(self.__class__, name)
        
    def __get_node_attr(self, name):
        return getattr(self, '__node_' + name)        
        
    def __getitem__(self, value):
        if str(value) in self.__nodes:
            return self.__nodes[str(value)]
        else:
            raise pkex.BasePyKatException("The node '%s' could not be found in the network." % str(value))
        
    def __contains__(self, value):
        return value in self.__nodes
    
    def __nodeSearch(self, fnode, currcomp, branches, tnode):
        
        if fnode == tnode:
            branches[-1][0] = True
            branches[-1][1] = True
            return True # Hurrah, we have found a path to the node
        elif fnode.isDump:
            branches[-1][0] = True
            return False # if the current node is a dump, we need to check another branch

        nextnode = None
        
        if isinstance(currcomp, pykat.components.beamSplitter):
            # if we get to a beamsplitter then we need to 
            # create new branches to search: the rfelected
            # and transmitted
            
            # set this branch as searched
            branches[-1][0] = True
            # add two new ones
            
            if fnode == currcomp.nodes[0]:
                rn = currcomp.nodes[1]
                tn = currcomp.nodes[2]
            elif fnode == currcomp.nodes[1]:
                rn = currcomp.nodes[0]
                tn = currcomp.nodes[3]
            elif fnode == currcomp.nodes[2]:
                rn = currcomp.nodes[3]
                tn = currcomp.nodes[0]
            elif fnode == currcomp.nodes[3]:
                rn = currcomp.nodes[2]
                tn = currcomp.nodes[1]
            else:
                raise pkex.BasePyKatException("Node not attached in path find to BS")
            
            nextcomp = None
            
            if tn.components[0] == currcomp:
                nextcomp = tn.components[1]
            elif tn.components[1] == currcomp:
                nextcomp = tn.components[0]
            
            if nextcomp != None:
                branches.append([False, False, tn, nextcomp, []])
            
            if rn.components[0] == currcomp:
                nextcomp = rn.components[1]
            elif rn.components[1] == currcomp:
                nextcomp = rn.components[0]
            
            if nextcomp != None:
                branches.append([False, False, rn, nextcomp, []])
            
            branches[-1][-1].append(currcomp)
            
            return False
            
        elif isinstance(currcomp, pykat.components.isolator):
            print ("isol")
        elif isinstance(currcomp, pykat.components.laser):
            # if we are at a laser then we can't go any further
            # and it isn;t this node as we checked before
            branches[-1][0] = True
            return False
        elif len(currcomp.nodes) == 2:
            if currcomp.nodes[0] == fnode:
                nextnode = currcomp.nodes[1]
            elif currcomp.nodes[1] == fnode:
                nextnode = currcomp.nodes[0]
            else:
                raise pkex.BasePyKatException("Unexpeceted condition")
        else:
            raise pkex.BasePyKatException("Did not handle component {0} correctly, has more or less than 2 nodes.".format(currcomp))
        
        if nextnode is None:
            branches[-1][0] = True
            return False
        elif nextnode == tnode:
            branches[-1][0] = True
            branches[-1][1] = True
            branches[-1][-1].append(currcomp)
            return True
        else:
            # Now we have the next node, we need the next component
            if nextnode.components[0] == currcomp:
                nextcomp = nextnode.components[1]
            elif nextnode.components[1] == currcomp:
                nextcomp = nextnode.components[0]
            else:
                raise pkex.BasePyKatException("Unexpeceted condition")

            if nextcomp is None:
                branches[-1][0] = True
                return False
            
            branches[-1][-1].append(currcomp)
            
            return self.__nodeSearch(nextnode, nextcomp, branches, tnode)
            
    def getComponentsBetween(self, from_node, to_node):
        """
        This function will trace the path between the two nodes specified and return a list
        of the components it finds between them.
        """
        
        if isinstance(from_node, six.string_types):
            from_node = self.__kat.nodes[from_node]
            
        if isinstance(from_node, NodeGaussSetter):
            from_node = from_node.node
            
        if isinstance(to_node, six.string_types):
            to_node = self.__kat.nodes[to_node]
            
        if isinstance(to_node, NodeGaussSetter):
            to_node = to_node.node
            
        if to_node == from_node:
            return []
    
        if from_node.name not in self.__nodes:
            raise pkex.BasePyKatException("Node {0} cannot be found in this kat object".format(from_node))

        if to_node.name not in self.__nodes:
            raise pkex.BasePyKatException("Node {0} cannot be found in this kat object".format(to_node))
        
        branches = []
        
        fn = self.__nodes[from_node.name]
        tn = self.__nodes[to_node.name]
        
        branches.append([False, False, fn, fn.components[1], []])
        branches.append([False, False, fn, fn.components[0], []])
        
        while len(branches) > 0 and branches[-1][1] != True:
            while branches[-1][0] == False:
                branch = branches[-1]
            
                if not self.__nodeSearch(branch[2], branch[3], branches, tn):
                    if len(branches) > 0 and branches[-1][0] != False:
                        branches.pop()
            
            if branches[-1][1] != True:
                while len(branches) > 0 and branches[-1][0] == True:
                    branches.pop()
            
            
        comps = []
        
        if len(branches) > 0  and branches[-1][0] == True and branches[-1][1] == True:
            # select the branches that form the path from node to node
            br = [b for b in branches if b[0] == True]
        
            for b in br:
                comps.extend(b[-1])
        
        return comps
    
    
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
        self._isDump = False
        
    def __str__(self): return self.__name

    @property
    def isDump(self): return self._isDump
        
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
    
    def removeGauss(self):
        self.__q_x = None
        self.__q_y = None
        self.__q_comp = None
    
    def setGauss(self, component, *args):
        self.__q_comp = component
        
        if len(args) == 1:  
            self.__q_x = BeamParam(self._network.kat.lambda0, q=args[0])
            self.__q_y = BeamParam(self._network.kat.lambda0, q=args[0])
        elif len(args) == 2:
            self.__q_x = BeamParam(self._network.kat.lambda0, q=args[0])
            self.__q_y = BeamParam(self._network.kat.lambda0, q=args[1])
        else:
            raise pkex.BasePyKatException("Must specify either 1 Gaussian beam parameter or 2 for astigmatic beams")
                
    def getFinesseText(self):    
        if self.__q_x is None or self.__q_y is None or self.__q_comp is None:
            return []
            
        rtn = []
        
        # to get the name of the gauss parameter is a bit convoluted...
        # firstly the name is set in the NodeGaussSetter object which is
        # connected to each component, so this has to be retrieved and 
        # then applied.
        if hasattr(self.__q_comp, self.name):
            ns = getattr(self.__q_comp, self.name)
            
            # if no name is present give it a default one
            if ns.name != None:
                name = ns.name
            else:
                name = "g_%s" % self.name
        else:
            raise pkex.BasePyKatException("Node {0} is not connected to {1}".format(self.name, self.__q_comp.name))
        
        
        if self.__q_x == self.__q_y:
            rtn.append("gauss {name} {comp} {node} {w0:.15g} {z:.15g}".format(name=name, node=self.name, comp=self.__q_comp.name, w0=self.__q_x.w0, z=self.__q_x.z))
        else:
            rtn.append("gauss {name} {comp} {node} {w0x:.15g} {zx:.15g} {w0y:.15g} {zy:.15g}".format(name=name, node=self.name, comp=self.__q_comp.name, w0x=self.__q_x.w0, zx=self.__q_x.z, w0y=self.__q_y.w0, zy=self.__q_y.z))
            
        return rtn
        
    def isConnected(self):
        if (self.components[0] is not None) and (self.components[1] is not None):
            return True
        else:
            return False
      
    def remove(self):
        self._network.removeNode(self)
        
        if self._item != None:
            self._item.scene().removeItem(self._item)
    
    def getQGraphicsItem(self,dx=0,dy=0,nsize=8,parent=None):
        if not USE_GUI:
            raise NoGUIException
            
        if self._item is None:
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
            if comps[1] is None:
                ix = -1
            else:
                ix = comps[1].nodes.index(self)
                
            return [True, comps[1], ix]
            
        elif obj == comps[1]:
            if comps[0] is None:
                ix = -1
            else:
                ix = comps[0].nodes.index(self)
                
            return [True, comps[0], ix]
        else:
            return [False, None]
        
    @property
    def name(self): return self.__name      
        
        
class DumpNode(Node):
    __total_dump_node_id = 0
    
    def __init__(self, network):
        DumpNode.__total_dump_node_id -= 1
        Node.__init__(self, 'dump', network, DumpNode.__total_dump_node_id)
        self._isDump = True
        
        
