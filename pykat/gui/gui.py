# -*- coding: utf-8 -*-
"""
Created on Tue Jan 29 11:35:48 2013

@author: Daniel
"""

from pykat.components import Component
from pykat.detectors import Detector

from PyQt4 import QtGui, QtCore
from PyQt4.Qt import *
from PyQt4.QtGui import QCursor, QGraphicsItem
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
        self.graphicsView = pyKatGraphicsView(self.centralwidget)
        self.graphicsView.setObjectName("graphicsView")
        self.graphicsView.setViewportUpdateMode(QGraphicsView.FullViewportUpdate)
        self.gridLayout.addWidget(self.graphicsView, 0, 0, 1, 1)
        
        # create a new scene
        if kat.scene == None: 
            kat.scene = pyKatGraphicsScene()  
        
        self.__scene = kat.scene
        
        # add scene to the graphics view
        self.graphicsView.setScene(self.__scene)
        self.graphicsView.setRenderHint(QtGui.QPainter.Antialiasing)
                
        self.actionExport_to_SVG.triggered.connect(lambda: self.exportToSVG())
        self.actionClose.triggered.connect(lambda: self.close)

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
                # uncomment this line to stop background bitmap caching of
                # svg rendering. Important to make sure when rendering to 
                # svg file that it is in a vector format. Gradients however
                # don't work...
                itm.setCacheMode(QGraphicsItem.NoCache)
                self.__scene.addItem(itm)
                
    def exportToSVG(self):
        self.statusbar.showMessage("Saving to 'output.svg'...")
        
        svg = QSvgGenerator()
        svg.setFileName("./output.svg")
        svg.setSize(QSize(self.__scene.width(), self.__scene.height()))
        svg.setViewBox(QRect(0,0,self.__scene.width(), self.__scene.height()))
        svg.setTitle("pyKat output of example.kat")
        
        pntr = QPainter(svg)
        self.__scene.render(pntr)
        pntr.end()
        
        self.statusbar.showMessage("Complete: Saved to 'output.svg'")
                
class pyKatGraphicsScene(QGraphicsScene):
    def drawBackground(self, painter, rect):
        size = 10
        painter.setPen(QPen(QColor(200,200,255,255),0.5))
        
        start = round(rect.top(), size)
        
        if start > rect.top():
            start =- size
        
        y = start - size
        
        while y < rect.bottom():
            y += size
            painter.drawLine(rect.left(),y, rect.right(), y)
        
        start = round(rect.left(), size)
        
        if start > rect.left():
            start =- size
        
        y = start - size
        
        while y < rect.right():
            y += size
            painter.drawLine(y, rect.top(), y, rect.bottom())
                    
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
            if isinstance(item, Component):           
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
            elif isinstance(item, NodeQGraphicItem):
                if item == None:
                    self.__selected_item = None
                    self.__prev_pt = None
                else:
                    if isinstance(item.parentItem(),SpaceQGraphicsItem):        
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
            #pt_ = self.__prev_pt
            pt = self.mapToScene(ev.pos())
            
            # smooth moving of item depending on where you click
            #item.moveBy(pt.x()-pt_.x(), pt.y()-pt_.y())
            # then snap to some integer value
            snap = 10.0
            
            
            item.setPos(int(round(pt.x()/snap)*snap),int(round(pt.y()/snap)*snap))
            self.__prev_pt = pt
            
            
            