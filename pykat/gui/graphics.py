# -*- coding: utf-8 -*-
"""
Created on Fri Feb 01 09:13:03 2013

@author: Daniel
"""

from __future__ import print_function
import pykat.external.six as six
if six.PY2:
	import exceptions
from PyQt4.QtGui import *
from PyQt4.Qt import *
from PyQt4 import QtSvg
from PyQt4.QtSvg import QGraphicsSvgItem
import pykat.components
import weakref

nsize = 10
    
class NodeQGraphicItem(QGraphicsRectItem):
    
    def __init__(self, node, x,y, *args, **kwargs):
        QGraphicsRectItem.__init__(self, *args, **kwargs)
        self.__node = node
        self.setPos(x,y)
        
        item = QGraphicsTextItem(node.name, self)
        rect = item.boundingRect()       
        item.setPos(-0.5*rect.width(), 0)
        
        self.setAcceptHoverEvents(True)
        
        self.marked = False
        
    @property
    def node(self): return self.__node
        
    def refresh(self):
        if not self.marked:
            if self.__node.isConnected():
                self.setBrush(QBrush(Qt.red))
            else:
                self.setBrush(QBrush(Qt.green))
        else:
            self.setBrush(QBrush(Qt.yellow))
            
class SpaceQGraphicsItem(QGraphicsLineItem):
    def __init__(self, spaceComponent):
        QGraphicsLineItem.__init__(self)
        self.__n1 = None
        self.__n2 = None
        self.__space = spaceComponent
    
        item = QGraphicsTextItem(self.__space.name, self)
        rect = item.boundingRect()       
        item.setPos(-0.5*rect.width(),0*rect.height())
        self.refresh()
        
    @property
    def space(self): return self.__space
    
    def refresh(self):
        nodes = self.__space.nodes
                
        conn = nodes[0].amIConnected(self.__space)
        
        x1 = 0
        y1 = 0
        x2 = 0
        y2 = 0
                
        if conn[0]:
            if conn[1] != None:
                if self.__n1 is not None:
                    # i.e. we have a node graphic item but now it is connected to something
                    self.__n1.scene().removeItem(self.__n1)
                    self.__n1 = None
                    
                # now check if a connected component was returned too
                if conn[1] != None:
                    # so this node should be attached to something
                    # in this case we get the position of their node 
                    # and draw the the line from their
                    itm=conn[1].getQGraphicsItem()
                    x1 = itm.x() + itm.nodedx[conn[2]][0]
                    y1 = itm.y() + itm.nodedx[conn[2]][1]
            else:
                if self.__n1 == None:
                    self.__n1 = NodeQGraphicItem(nodes[0],0,0,-nsize/2,-nsize/2,nsize,nsize,self)
                    self.__n1.setPen(QPen(Qt.black,1))
                    
                self.__n1.setVisible(True)
                self.__n1.setBrush(QBrush(Qt.green))
                p = self.__n1.pos()
                x1 = self.x()+p.x()
                y1 = self.y()+p.y()
                
        conn = nodes[1].amIConnected(self.__space)
        
        if conn[0]:
            if conn[1] != None:
                
                if self.__n2 is not None:
                    # i.e. we have a node graphic item but now it is connected to something
                    self.__n2.scene().removeItem(self.__n2)
                    self.__n2 = None
                    
                # now check if a connected component was returned too
                if conn[1] != None:
                    # so this node should be attached to something
                    # in this case we get the position of their node 
                    # and draw the the line from their
                    itm=conn[1].getQGraphicsItem()
                    x2 = itm.x() + itm.nodedx[conn[2]][0]
                    y2 = itm.y() + itm.nodedx[conn[2]][1]
            else:
                if self.__n2 == None:
                    self.__n2 = NodeQGraphicItem(nodes[1],0,0,-nsize/2,-nsize/2,nsize,nsize,self)
                    self.__n2.setPen(QPen(Qt.black,1))
                    
                self.__n2.setVisible(True)
                self.__n2.setBrush(QBrush(Qt.green))
                p = self.__n2.pos()
                x2 = self.x()+p.x()
                y2 = self.y()+p.y()
        
        # convert x1,y1,x2 and y2 into the local coordinates of the 
        # space object
        p = QPointF((x1-x2)*0.5,(y1-y2)*0.5)
        self.setPos(x1 - p.x(), y1 - p.y())
        
        # if the nodes are visible then reposition them in the 
        # component reference frame
        if self.__n1 is not None and self.__n1.isVisible():
            self.__n1.setPos(QPointF(p.x(),p.y()))
            self.__n1.refresh()
        
        if self.__n2 is not None and self.__n2.isVisible():
            self.__n2.setPos(QPointF(p.x()+x2-x1, p.y()+y2-y1))
            self.__n2.refresh()
            
        self.setLine(p.x(), p.y(), p.x()+x2-x1, p.y()+y2-y1)
        self.setPen(QPen(Qt.red, 3))
        
class ComponentQGraphicsItem(QGraphicsObject): #(QtSvg.QGraphicsSvgItem):
        
    def __init__(self, svgfile, component, nodes):
        #QGraphicsSvgItem.__init__(self,svgfile)
        QGraphicsObject.__init__(self)
        self.__svggraphic = QGraphicsSvgItem(svgfile) 
        
        rect = self.__svggraphic.boundingRect()
        
        self.__nodeGraphics = []
        self.__component = weakref.ref(component)
        
        # this signals the itemChange() method when this item is moved
        # used for refreshing the spaces between components
        self.setFlags(QGraphicsItem.ItemSendsGeometryChanges)
        self.nodedx = [] # stores the node square offsets
                
        item = QGraphicsTextItem(component.name,self)
        rect = item.boundingRect()       
        item.setPos(-0.5*rect.width(),40-0.5*rect.height())
        
        self.setAcceptsHoverEvents(True)
        
        for n in nodes:
            self.nodedx.append([n[0],n[1]])
            node = n[2].getQGraphicsItem(n[0],n[1],nsize,self)
            node.setPen(QPen(Qt.black))
            node.refresh()
            self.__nodeGraphics.append(node)
            
        self.refresh()
        self.installEventFilter(self)
        self.setHandlesChildEvents(True)

    def boundingRect(self):
        return self.__svggraphic.boundingRect()

    def paint(self, arg1, arg2, arg3):
        self.__svggraphic.rotate(45)
        self.__svggraphic.paint( arg1, arg2, arg3)
        
    @property
    def component(self): return self.__component()
    
    def refresh(self):
        for n in self.__nodeGraphics:
            n.refresh()
        
    def itemChange(self, change, value):
        # if the item is moved then update any spaces attached to it
        if change == QGraphicsItem.ItemPositionHasChanged:
            nodes = self.component.nodes
            
            for n in nodes:
                conn = n.amIConnected(self.component)
                
                if conn[0] and isinstance(conn[1],  pykat.components.space):
                    conn[1].getQGraphicsItem().refresh()
                   
        return QGraphicsSvgItem.itemChange(self, change, value)
            
