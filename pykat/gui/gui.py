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
import functools

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
        self.graphicsView.viewport().setMouseTracking(True)
        self.setMouseTracking(True)
        
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

        self.kat = kat        
        
    def main(self):
        self.show()
        
        self.addComponentsToScene()
        
    def scene(self):
        return self.__scene

    def addComponentsToScene(self):
        
        for c in self.kat.getComponents():
            self.addComponentToScene(c)
                
    def addComponentToScene(self,c,x=0,y=0):
        itm = c.getQGraphicsItem()
            
        if itm != None:
            itm.setPos(x,y)
            # uncomment this line to stop background bitmap caching of
            # svg rendering. Important to make sure when rendering to 
            # svg file that it is in a vector format. Gradients however
            # don't work...
            itm.refresh()
            itm.setCacheMode(QGraphicsItem.NoCache)
            self.__scene.addItem(itm)
            
    def exportToSVG(self):
        self.statusbar.showMessage("Saving to 'output.svg'...")
        
        svg = QSvgGenerator()
        filename = QtGui.QFileDialog.getSaveFileNameAndFilter(self,'Save SVG File',filter=".svg")
        
        if filename == None:
            return None
        
        svg.setFileName(filename)
        svg.setSize(QSize(self.__scene.width(), self.__scene.height()))
        svg.setViewBox(QRect(0,0,self.__scene.width(), self.__scene.height()))
        svg.setTitle("pyKat output of example.kat")
        
        pntr = QPainter(svg)
        self.__scene.render(pntr)
        pntr.end()
        
        self.statusbar.showMessage("Complete: Saved to 'output.svg'")
    
    def addMirror(self, x,y):
        name = self.kat.getNewComponentName('m')
        n = self.kat.getNewNodeNames('n',2)
        m = pykat.components.mirror(self.kat,name,n[0],n[1])
               
        self.addComponentToScene(m,x,y)
                
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
        self.setMouseTracking(True)
        self.__itemHover = None
        self.__marked = None
        
    def contextMenuEvent(self, ev):  
        pt = self.mapToScene(ev.pos())
        
        gui = self.parentWidget().parent() # get the main gui window
        
        menu = QMenu(self)
        addmenu = menu.addMenu("Add...")
        
        action = addmenu.addAction("Mirror")
        action.triggered.connect(functools.partial(gui.addMirror,pt.x(),pt.y()))
                
        addmenu.addAction("Laser")
        addmenu.addAction("Beamsplitter")
        addmenu.addAction("Photodiode")
        
        item = self.scene().itemAt(pt.x(),pt.y())
        
        print pt.x(),pt.y(),item
        
        if item is not None :
            if isinstance(item, ComponentQGraphicsItem):           
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
                
                 
                if isinstance(item.parentItem(),SpaceQGraphicsItem):        
                    self.__selected_item = item
                    self.__prev_pt = pt
                    
                elif isinstance(item.parentItem(),ComponentQGraphicsItem):
                    self.__selected_item = item.parentItem()
                    self.__prev_pt = pt
            
            if self.__selected_item is not None:
                self.setCursor(QCursor(Qt.ClosedHandCursor))
                        
    def mouseReleaseEvent(self, ev):
        # if we have dragged a node and marked another to connect it too
        if self.__selected_item is not None and isinstance(self.__selected_item, NodeQGraphicItem) and self.__marked is not None:
            # node attached to space which needs to be removed
            node_s = self.__selected_item.node
            
            # get the selected node which must be attached to a space, which is the nodes parent
            # and then get the space component from it
            qspace = self.__selected_item.parentItem()
            space = qspace.space
            
            # marked node, then get the parent object and component
            node_c = self.__marked.node
            qcomp = self.__marked.parentItem()
            
            # connect space of node dragged to the component node
            space.changeNode(node_s,node_c)
            node_s.remove()
            
            # then refresh the graphical items
            qspace.refresh()
            qcomp.refresh()
            
        if self.__marked is not None:
            self.__marked.marked = False
            self.__marked.refresh()
            self.__marked = None
            
        self.__selected_item = None
        self.setCursor(QCursor(Qt.ArrowCursor))
            
    def mouseMoveEvent(self, ev):
        if self.__selected_item != None:
            
            item = self.__selected_item
            #pt_ = self.__prev_pt
            pt = self.mapToScene(ev.pos())
            
            # smooth moving of item depending on where you click
            #item.moveBy(pt.x()-pt_.x(), pt.y()-pt_.y())
            # then snap to some integer value
            snap = 10.0
            
            # if we are moving a node, it must be attached to a space
            # component otherwise we shouldn't be here
            if isinstance(item, NodeQGraphicItem) and isinstance(item.parentItem(), SpaceQGraphicsItem):
                space = item.parentItem()
                
                item.setPos(pt.x()-space.x(),pt.y()-space.y())
                space.refresh()
                
                # now check to see if any other connectable nodes are within reach
                # and if so hightlight them
                select_size = 20
                rect = QRectF(pt.x()-select_size/2,pt.y()-select_size/2,select_size,select_size)
                itms = item.scene().items(rect)
                
                # remove the node we are dragging
                if item in itms:
                    itms.remove(item)
                
                if self.__marked is not None:
                    self.__marked.marked = False
                    self.__marked.refresh()
                    self.__marked = None
                    
                if len(itms) > 0:
                    for i in itms:
                        if isinstance(i,NodeQGraphicItem) and i != item:
                            i.marked = True
                            i.refresh()
                            self.__marked = i
                            break
            else:
                item.setPos(int(round(pt.x()/snap)*snap),int(round(pt.y()/snap)*snap))
                
            self.__prev_pt = pt
            
            return
        
        else:
            item = self.itemAt(ev.pos())
            
            if isinstance(item, (NodeQGraphicItem, ComponentQGraphicsItem)) or item is None:
                #if item is not None or self.__itemHover is not None:
                    
                if isinstance(item, ComponentQGraphicsItem):
                    self.__itemHover = item
                elif isinstance(item,NodeQGraphicItem) and (not item.node.isConnected()):
                    self.__itemHover = item
                else:
                    self.__itemHover = None
            
            if self.__itemHover is not None:                
                self.setCursor(QCursor(Qt.OpenHandCursor))
            else:
                self.setCursor(QCursor(Qt.ArrowCursor))
                        
            
            