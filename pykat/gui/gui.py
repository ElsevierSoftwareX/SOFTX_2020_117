# -*- coding: utf-8 -*-
"""
Created on Tue Jan 29 11:35:48 2013

@author: Daniel
"""

from PyQt4 import QtGui, QtCore
from PyQt4.Qt import *
from PyQt4.QtGui import QCursor
from pykat.gui.graphics import *
import qt_gui
     
def openGUI(kat):
    app = QtGui.QApplication([""])
    pykatgui = pyKatGUI(kat)
    pykatgui.main()
    app.exec_()
        
class pyKatGUI(QtGui.QMainWindow, qt_gui.Ui_MainWindow):
    def __init__(self, kat,parent=None):
        super(pyKatGUI, self).__init__(parent)
        self.setupUi(self)
        
        # create a new scene
        self.__scene = QGraphicsScene()  

        brush = QBrush()
        brush.setStyle(Qt.CrossPattern)
        brush.setColor(QColor(230,230,230))
        self.__scene.setBackgroundBrush(brush)
        
        # add scene to the graphics view
        self.graphicsView.setScene(self.__scene)
                
        self._kat = kat        
        
    def main(self):
        self.show()
        
        self.addComponentsToScene()
        
    def scene(self):
        return self.__scene

    def addComponentsToScene(self):
        
        for c in self._kat.getComponents():
            itm = c.getQGraphicsItem()
            
            if itm != None:
                itm.setPos(0,0)
                self.__scene.addItem(itm)
                
class pyKatGraphicsView(QGraphicsView):
    def __init__(self,val):
        QGraphicsView.__init__(self,val)
        
        self.__selected_item = None
        self.__prev_pt = None
            
    def contextMenuEvent(self, ev):  
        pt = self.mapToScene(ev.pos())
        
        menu = QMenu(self)
        addmenu = menu.addMenu("Add...")
        addmenu.addAction("Mirror")
        addmenu.addAction("Laser")
        addmenu.addAction("Beamsplitter")
        addmenu.addAction("Photodiode")
        
        item = self.itemAt(pt.x(),pt.y())
        
        if item != None :
            if isinstance(item,Component):           
                menu.addSeparator()
                menu.addAction("Edit")
                menu.addAction("Delete")
            if isinstance(item,NodeQGraphicItem):
                menu.addSeparator()
                menu.addAction("Disconnect")

        menu.popup(ev.globalPos());        
        
    def mousePressEvent(self, ev):
        
        if ev.button()==Qt.LeftButton:
            
            pt = self.mapToScene(ev.pos())
            
            item = self.scene().itemAt(pt)
            
            if isinstance(item, ComponentQGraphicsItem):
                if item == None:
                    self.__selected_item = None
                    self.__prev_pt = None
                else:
                    item.setFocus(Qt.MouseFocusReason)
        
                    self.__selected_item = item
                    self.__prev_pt = pt
            
    def mouseReleaseEvent(self, ev):
        self.__selected_item = None
        self.setCursor(QCursor(Qt.ArrowCursor))
        pass
        
    def mouseMoveEvent(self, ev):
        if self.__selected_item != None:
            self.setCursor(QCursor(Qt.ClosedHandCursor))
            
            item = self.__selected_item
            pt_ = self.__prev_pt
            pt = self.mapToScene(ev.pos())
            
            item.moveBy(pt.x()-pt_.x(), pt.y()-pt_.y())
            self.__prev_pt = pt