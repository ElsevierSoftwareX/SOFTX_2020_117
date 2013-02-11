# -*- coding: utf-8 -*-
"""
Created on Fri Feb 01 09:13:03 2013

@author: Daniel
"""

from PyQt4.QtGui import *
from PyQt4.Qt import *
import pykat.components

nsize = 8

class NodeQGraphicItem(QGraphicsRectItem):
    
    def __init__(self, node, x,y, *args, **kwargs):
        QGraphicsRectItem.__init__(self, *args, **kwargs)
        self.__node = node
        
        self.setPos(x,y)
        
        item = QGraphicsTextItem(node.name, self)
        rect = item.boundingRect()       
        item.setPos(-0.5*rect.width(), 0)
        
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
        
    def refresh(self):    
        nodes = self.__space.getNodes()
        
        if self.__n1 == None:
            self.__n1 = NodeQGraphicItem(nodes[0],0,0,-nsize/2,-nsize/2,nsize,nsize,self)
            self.__n1.setPen(QPen(Qt.black))
            
        conn = nodes[0].amIConnected(self.__space)
        
        if conn[0]:
            if conn[1] != None:
                self.__n1.setVisible(False)
                # now check if a connected component was returned too
                if conn[1] != None:
                    # so this node should be attached to something
                    # in this case we get the position of their node 
                    # and draw the the line from their
                    itm=conn[1].getQGraphicsItem()
                    x1 = itm.x() + itm.nodedx[conn[2]][0]
                    y1 = itm.y() + itm.nodedx[conn[2]][1]
            else:
                self.__n1.setVisible(True)
                self.__n1.setBrush(QBrush(Qt.red))
                x1 = 0
                y1 = 0
                
        if self.__n2 == None:
            self.__n2 = NodeQGraphicItem(nodes[1],0,0,-nsize/2,-nsize/2,nsize,nsize,self)
            self.__n2.setPen(QPen(Qt.black))
                
        conn = nodes[1].amIConnected(self.__space)
        
        if conn[0]:
            if conn[1] != None:
                self.__n2.setVisible(False)
                # now check if a connected component was returned too
                if conn[1] != None:
                    # so this node should be attached to something
                    # in this case we get the position of their node 
                    # and draw the the line from their
                    itm=conn[1].getQGraphicsItem()
                    x2 = itm.x() + itm.nodedx[conn[2]][0]
                    y2 = itm.y() + itm.nodedx[conn[2]][1]
            else:
                self.__n2.setVisible(True)
                self.__n2.setBrush(QBrush(Qt.red))
                p = self.__n2.pos()
                x2 = -p.x()
                y2 = -p.y()
            
        p = QPointF((x1-x2)*0.5,(y1-y2)*0.5)
        self.setPos(x1 - p.x(), y1 - p.y())
        
        # if the nodes are visible then reposition them in the 
        # component reference frame
        if self.__n1.isVisible():
            self.__n1.setPos(QPointF(p.x(),p.y()))
        
        if self.__n2.isVisible():
            self.__n2.setPos(QPointF(p.x()+x2-x1, p.y()+y2-y1))
            
        self.setLine(p.x(), p.y(), p.x()+x2-x1, p.y()+y2-y1)
        self.setPen(QPen(Qt.red, 3))
        
        
class ComponentQGraphicsItem(QGraphicsSvgItem):
    
    def __init__(self, svgfile, component, nodes):
        QGraphicsSvgItem.__init__(self,svgfile)
        
        self.__component = component
        # this signals the itemChange() method when this item is moved
        # used for refreshing the spaces between components
        self.setFlags(QGraphicsItem.ItemSendsGeometryChanges)
        self.nodedx = [] # stores the node square offsets
                
        item = QGraphicsTextItem(component.name,self)
        rect = item.boundingRect()       
        item.setPos(-0.5*rect.width(),40-0.5*rect.height())
        
        for n in nodes:
            self.nodedx.append([n[0],n[1]])
            node = NodeQGraphicItem(n[2],n[0],n[1],-nsize/2,-nsize/2,nsize,nsize,self)
            node.setBrush(QBrush(Qt.red))
            node.setPen(QPen(Qt.black))
        
    def itemChange(self, change, value):
        # if the item move then update any spaces
        if change == QGraphicsItem.ItemPositionHasChanged:
            nodes = self.__component.getNodes()
            
            for n in nodes:
                conn = n.amIConnected(self.__component)
                
                if conn[0] and isinstance(conn[1],  pykat.components.space):
                    conn[1].getQGraphicsItem().refresh()
                   
        return QGraphicsSvgItem.itemChange(self, change, value)
            