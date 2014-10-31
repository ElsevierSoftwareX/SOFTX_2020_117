# -*- coding: utf-8 -*-
"""
Created on Tue Jan 29 11:35:48 2013

@author: Daniel
"""

from pykat.components import Component, space
from pykat.detectors import BaseDetector as Detector

from PyQt4 import QtGui, QtCore
from PyQt4.Qt import *
from PyQt4.QtGui import QCursor, QGraphicsItem
from pykat.gui.graphics import *
import qt_gui
import functools

class pyKatGUI(QtGui.QMainWindow, qt_gui.Ui_MainWindow):
    def __init__(self, kat, parent=None):
        super(pyKatGUI, self).__init__(parent)
        
        self.setupUi(self)
        self.graphicsView = pyKatGraphicsView(self.centralwidget, kat)
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

        self._kat = kat        
        
    @property
    def kat(self): return self._kat
    
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
    
    def _onComponentRemoved(self, comp, nodes):
        """
        When a component has been removed from the kat object this function should update
        all gui objects. 
            comp - object that is removed
            nodes - nodes that this comp was attached too, as that information may no longer be accessible
        """
        itm = comp.getQGraphicsItem()
            
        if itm != None:

            itm.refresh()
            self.__scene.removeItem(itm)
            
            for n in nodes:
                for cc in self._kat.nodes.getNodeComponents(n):
                    if cc != None:
                        ccitm = cc.getQGraphicsItem()
                        if ccitm != None:
                            ccitm.refresh()
            
            
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
        m = pykat.components.mirror(name,n[0],n[1],R=0.5,T=0.5)
        
        self.kat.add(m)
        self.addComponentToScene(m,x,y)
    
    def addBeamsplitter(self, x, y):
        name = self.kat.getNewComponentName('bs')
        n = self.kat.getNewNodeNames('n', 4)
        m = pykat.components.beamSplitter(name,n[0],n[1],n[2],n[3],R=0.5,T=0.5)
        
        self.kat.add(m)
        self.addComponentToScene(m,x,y)
        
    def addSpace(self, x,y):
        name = self.kat.getNewComponentName('s')
        n = self.kat.getNewNodeNames('n',2)
        s = pykat.components.space(name, n[0], n[1])

        self.kat.add(s)
        self.addComponentToScene(s,x,y) 
     
    def addLaser(self, x,y):
        name = self.kat.getNewComponentName('l')
        n = self.kat.getNewNodeNames('n',1)
        l = pykat.components.laser(name, n[0])

        self.kat.add(l)
        self.addComponentToScene(l,x,y)   
    
    def addModulator(self, x,y):
        name = self.kat.getNewComponentName('mod')
        n = self.kat.getNewNodeNames('n',2)
        l = pykat.components.modulator(name, n[0], n[1], 1e6, 0.1, 1)

        self.kat.add(l)
        self.addComponentToScene(l,x,y)
    
    def addLens(self, x,y):
        name = self.kat.getNewComponentName('lens')
        n = self.kat.getNewNodeNames('n',2)
        l = pykat.components.lens(name, n[0], n[1])

        self.kat.add(l)
        self.addComponentToScene(l,x,y)
        
    def addIsolator(self, x,y):
        name = self.kat.getNewComponentName('isol')
        n = self.kat.getNewNodeNames('n',3)
        l = pykat.components.isolator(name, n[0], n[1], node3=n[2])

        self.kat.add(l)
        self.addComponentToScene(l,x,y)
            
    def addPhotodiode(self, x, y):
        name = self.kat.getNewDetectorName('pd')
        n = self.kat.getNewNodeNames('n',1)
        l = pykat.detectors.photodiode(name, n[0], [])

        self.kat.add(l)
        self.addComponentToScene(l,x,y)   
    
    def deleteComponent(self, comp):
        comp.component.remove()
    
    def disconnect(self, node):
        comps = self.kat.nodes.getNodeComponents(node)
        
        spaces = [c for c in comps if isinstance(c, space)]
        
        if len(spaces) > 0:
            dis_comp = spaces[0]
        else:
            dis_comp = comps[0]
        
        new_node_name = self.kat.getNewNodeNames("n", 1)
        new_node = self.kat.nodes.createNode(new_node_name[0])
        
        self.kat.nodes.replaceNode(dis_comp, node, new_node)
        
        # refresh all the graphics that might be affected
        for c in node.components + new_node.components:
            if c != None:
                c.getQGraphicsItem().refresh()
    
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
    def __init__(self, val, kat):
        QGraphicsView.__init__(self, val)
        self._kat = kat
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
        
        action = addmenu.addAction("Space")
        action.triggered.connect(functools.partial(gui.addSpace, pt.x(), pt.y()))
        
        action = addmenu.addAction("Mirror")
        action.triggered.connect(functools.partial(gui.addMirror, pt.x(), pt.y()))
        
        action = addmenu.addAction("Beamsplitter")
        action.triggered.connect(functools.partial(gui.addBeamsplitter, pt.x(), pt.y()))
        
        action = addmenu.addAction("Laser")
        action.triggered.connect(functools.partial(gui.addLaser, pt.x(), pt.y()))
        
        action = addmenu.addAction("Lens")
        action.triggered.connect(functools.partial(gui.addLens, pt.x(), pt.y()))
        
        action = addmenu.addAction("Isolator")
        action.triggered.connect(functools.partial(gui.addIsolator, pt.x(), pt.y()))
        
        action = addmenu.addAction("Modulator")
        action.triggered.connect(functools.partial(gui.addModulator, pt.x(), pt.y()))
        
        action = addmenu.addAction("Photodiode")
        action.triggered.connect(functools.partial(gui.addPhotodiode, pt.x(), pt.y()))
        
        item = self.scene().itemAt(pt.x(),pt.y())
        
        print pt.x(),pt.y(),item
        
        if item is not None :
            if isinstance(item, ComponentQGraphicsItem):           
                menu.addSeparator()
                menu.addAction("Edit")
                action = menu.addAction("Delete")
                action.triggered.connect(functools.partial(gui.deleteComponent, item))
            if isinstance(item,NodeQGraphicItem):
                menu.addSeparator()
                comps = self._kat.nodes.getNodeComponents(item.node)
                
                if(comps.count(None) == 0):
                    action = menu.addAction("Disconnect")
                    action.triggered.connect(functools.partial(gui.disconnect, item.node))

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
            # the space node that has been dragged gets deleted and we
            # replace it with the components
            self._kat.nodes.replaceNode(space, node_s, node_c)
            
            # then refresh the graphical items
            qspace.refresh()
            qcomp.refresh()

            self.setCursor(QCursor(Qt.ArrowCursor))
                        
        if self.__marked is not None:
            self.__marked.marked = False
            self.__marked.refresh()
            self.__marked = None
            
        self.__selected_item = None
        
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
                        
            
            