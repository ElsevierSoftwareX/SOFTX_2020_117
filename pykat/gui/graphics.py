# -*- coding: utf-8 -*-
"""
Created on Fri Feb 01 09:13:03 2013

@author: Daniel
"""

from PyQt4.QtGui import *
from PyQt4.Qt import *

class NodeQGraphicItem(QGraphicsRectItem):
    pass

class ComponentQGraphicsItem(QGraphicsSvgItem):
    
    def __init__(self, svgfile, component, nodes):
        QGraphicsSvgItem.__init__(self, svgfile)
        self.__component = component
                
        item = QGraphicsTextItem(component.name,self)
        rect = item.boundingRect()       
        item.setPos(-0.5*rect.width(),40-0.5*rect.height())
        
        for n in nodes:
            node = NodeQGraphicItem(n[0],n[1],8,8,self)
            node.setBrush(QBrush(Qt.red))
            node.setPen(QPen(Qt.black))