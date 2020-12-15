import numpy as np
import random
import tfim
import itertools as it
import networkx as nx
from line_profiler import LineProfiler

import threading
import sys, os
from itertools import combinations
from PyQt5 import QtCore, QtWidgets
from PyQt5.QtWidgets import QMainWindow, QWidget, QLabel, QLineEdit
from PyQt5.QtWidgets import QPushButton
from PyQt5.QtCore import QSize
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QPainter, QBrush, QPen, QFont, QColor
from PyQt5.QtWidgets import QApplication



class Window(QMainWindow):
    def __init__(self):
        QMainWindow.__init__(self)

        
        
        
        #Set these to True to see periodic boundary conditions or Energy info, respectively
        ######################################################################
        self.PBC = True
        self.draw_lattice = False
        self.spins = False
        ######################################################################

        #Window size and label
        self.title = "SGViz_InfD"
        self.top= 50
        self.left= 100
        self.width = 1000
        self.height = 750

        self.yheight = 4
        self.xwidth = 4
        self.spacing = 35
        self.L = [self.yheight, self.xwidth]
        self.lattice = tfim.Lattice(self.L, self.PBC)
        self.N = self.lattice.N
        self.seed = 150
        if self.draw_lattice == False:
            self.bonds = bond_list(self)
            self.Jij = make_Jij(self)
        self.cnfg = 0
        self.configuration = list(map(int,list(bin(self.cnfg)[2:].zfill(self.N))))
        self.center = (500,375)
        self.coordList = spinCoords(self)
        self.plaq = make_plaquettes(self)
        self.f_plaq = frustrated(self)
        
        self.node_pairs = plaq_pairing(self)
        init_ground = initial_ground(self)
        self.p_pairings = init_ground[0]
        self.ground_distance = init_ground[1]
        
        v_edges = viable_edges(self)
        self.plaq_pairs = v_edges[0]
        self.str_length = v_edges[1]
        self.edges = v_edges[2]
        self.matchings = plaq_groups(self)
        self.string_groups = add_all_strings(self)
        self.string_index = 0
        self.b_bonds = broken_bonds(self)
        true_ground = make_config(self)
        self.ground_config = true_ground[0]
        self.true_ground_strings = true_ground[1]
        self.number_ground_states = len(self.true_ground_strings)
        
        if self.number_ground_states == 0:
            found_ground = False
        else:
            found_ground = True
            
        while found_ground == False:
            self.ground_distance += 1
            v_edges = viable_edges(self)
            self.plaq_pairs = v_edges[0]
            self.str_length = v_edges[1]
            self.edges = v_edges[2]
            if len(self.edges) != 0:
                self.matchings = plaq_groups(self)
                self.string_groups = add_all_strings(self)
                self.string_index = 0
                self.b_bonds = broken_bonds(self)
                true_ground = make_config(self)
                self.ground_config = true_ground[0]
                self.true_ground_strings = true_ground[1]
                self.number_ground_states = len(self.true_ground_strings)
                if self.number_ground_states != 0:
                    found_ground = True
                

        self.setStyleSheet("background-color: white;")
        
        self.Clabel = QLabel(self)
        self.Clabel.setFont(QFont('Decorative', 13))
        self.Clabel.setStyleSheet("color: black")
        self.Clabel.resize(400,30)
        self.Clabel.move(500,20)
        self.Clabel.setText('Configuration: ' + str(self.configuration))
        
        self.countLabel = QLabel(self)
        self.countLabel.setStyleSheet("color: black")
        self.countLabel.move(20,20)
        self.countLabel.resize(400,30)
        self.countLabel.setText("Number of ground states: " + str(self.number_ground_states))
        
        self.cEdit = QLineEdit(self)
        self.cEdit.setStyleSheet("color: black")
        self.cEdit.move(80, 650)
        self.cEdit.resize(100,32)
        self.cLabel = QLabel(self)
        self.cLabel.setStyleSheet("color: black")
        self.cLabel.setText("Config")
        self.cLabel.move(20,650)
        self.cLabel.resize(50,32)

        self.sEdit = QLineEdit(self)
        self.sEdit.setStyleSheet("color: black")
        self.sEdit.move(80, 600)
        self.sEdit.resize(100,32)
        self.sLabel = QLabel(self)
        self.sLabel.setStyleSheet("color: black")
        self.sLabel.setText("Seed")
        self.sLabel.move(20,600)
        self.sLabel.resize(50,32)

        
        self.pEdit = QLineEdit(self)
        self.pEdit.setStyleSheet("color: black")
        self.pEdit.move(850, 650)
        self.pEdit.resize(100,32)
        self.pLabel = QLabel(self)
        self.pLabel.setStyleSheet("color: black")
        self.pLabel.setText("PBC")
        self.pLabel.move(740,650)
        self.pLabel.resize(75,32)
        
        self.dlEdit = QLineEdit(self)
        self.dlEdit.setStyleSheet("color: black")
        self.dlEdit.move(850, 600)
        self.dlEdit.resize(100,32)
        self.dlLabel = QLabel(self)
        self.dlLabel.setStyleSheet("color: black")
        self.dlLabel.setText("Draw Lattice")
        self.dlLabel.move(740,600)
        self.dlLabel.resize(75,32)
        
        
        self.spEdit = QLineEdit(self)
        self.spEdit.setStyleSheet("color: black")
        self.spEdit.move(850, 550)
        self.spEdit.resize(100,32)
        self.spLabel = QLabel(self)
        self.spLabel.setStyleSheet("color: black")
        self.spLabel.setText("Draw Spins")
        self.spLabel.move(740,550)
        self.spLabel.resize(75,32)
        
        self.siEdit = QLineEdit(self)
        self.siEdit.setStyleSheet("color: black")
        self.siEdit.move(850, 500)
        self.siEdit.resize(100,32)
        self.siLabel = QLabel(self)
        self.siLabel.setStyleSheet("color: black")
        self.siLabel.setText("String Index")
        self.siLabel.move(740,500)
        self.siLabel.resize(75,32)
        

        self.xwidthEdit = QLineEdit(self)
        self.xwidthEdit.setStyleSheet("color: black")
        self.xwidthEdit.move(80, 550)
        self.xwidthEdit.resize(100,32)
        self.xwidthLabel = QLabel(self)
        self.xwidthLabel.setStyleSheet("color: black")
        self.xwidthLabel.setText("Width")
        self.xwidthLabel.move(20,550)
        self.xwidthLabel.resize(50,32)
        
        self.yheightEdit = QLineEdit(self)
        self.yheightEdit.setStyleSheet("color: black")
        self.yheightEdit.move(80, 500)
        self.yheightEdit.resize(100,32)
        self.yheightLabel = QLabel(self)
        self.yheightLabel.setStyleSheet("color: black")
        self.yheightLabel.setText("Height")
        self.yheightLabel.move(20,500)
        self.yheightLabel.resize(50,32)
        

        self.scountLabel = QLabel(self)
        self.scountLabel.setStyleSheet("color: black")
        self.scountLabel.setText("Unsatisfied: ")
        self.scountLabel.move(20,50)
        self.scountLabel.resize(100,32)
        self.scountLabel.setText("Unsatisfied: " + str(self.str_length))

        self.pybutton = QPushButton('Enter', self)
        self.pybutton.setStyleSheet("QPushButton {color: white; background-color: black}")
        self.pybutton.clicked.connect(self.clickMethod)
        self.pybutton.resize(100,32)
        self.pybutton.move(300, 650)

        
        self.InitWindow()

    def keyPressEvent(self, qKeyEvent):
        if qKeyEvent.key() == QtCore.Qt.Key_Return: 
            self.clickMethod()


    def clickMethod(self):
        thread = threading.Thread(target=self.DynMethod())
        thread.start()


    def InitWindow(self):
        self.setWindowTitle(self.title)
        self.setGeometry(self.left,self.top, self.width, self.height)
        self.show()

    def clickMethod(self):
        thread = threading.Thread(target=self.DynMethod())
        thread.start()

    def clearMethod(self):
        if not self.clear:
            self.Clabel.move(1000,20)
            self.cEdit.move(1000, 700)
            self.cLabel.move(1000,700)
            self.sEdit.move(1000, 650)
            self.sLabel.move(1000,650)
            self.NEdit.move(1000, 600)
            self.NLabel.move(1000,600)
            self.pybutton.move(1000, 700)
            self.pLabel.move(1000,800)
            self.pEdit.move(1000,900)
            self.countLabel.move(1000,200)
            self.scountLabel.move(1000,250)
            self.clear = True
        else:
            self.Clabel.move(700,20)
            self.cEdit.move(80, 700)
            self.cLabel.move(20,700)
            self.sEdit.move(80, 650)
            self.sLabel.move(20,650)
            self.NEdit.move(80, 600)
            self.NLabel.move(20,600)
            self.pybutton.move(300, 700)
            self.pLabel.move(740,650)
            self.pEdit.move(850, 650)
            self.countLabel.move(20,200)
            self.scountLabel.move(20,250)
            self.clear = False




    def DynMethod(self):

        change = False

        ws = self.xwidthEdit.text()
        hs = self.yheightEdit.text()
        

        if len(ws) > 0 or len(hs) > 0:
            if len(ws)>0:
                ws = int(ws)
            else:
                ws = self.xwidth
            if len(hs) > 0:
                hs = int(hs)
            else:
                hs = self.yheight
            if ws != self.xwidth or hs!= self.yheight:
                change = True
                self.xwidth = ws
                self.yheight = hs
                self.L = [self.yheight, self.xwidth]
                self.coordList = spinCoords(self)
                
                
        ss = self.sEdit.text()
        if len(ss) > 0:
            s = int(ss)
            if s != self.seed:
                change = True
                self.seed = s
                ss = self.sEdit.text()

        cs = self.cEdit.text()
        if len(cs) > 0:
            c = int(cs)
            if c != self.cnfg:
                change = True
                self.cnfg = c
                
        ps = self.pEdit.text()
        if len(ps) > 0:
            if str(ps) != str(self.PBC):
                change = True
                if self.PBC == False:
                    self.PBC = True
                else:
                    self.PBC = False
                
        dl = self.dlEdit.text()
        if len(dl) > 0:
            if str(dl) != str(self.draw_lattice):
                change = True
                if self.draw_lattice == True:
                    self.draw_lattice = False
                else:
                    self.draw_lattice = True
                    
        spin = self.spEdit.text()
        if len(spin) > 0:
            if str(spin) != str(self.spins):
                change = True
                if self.spins == True:
                    self.spins = False
                else:
                    self.spins = True
                    
        string_set = False
        string = self.siEdit.text()
        if len(string) > 0:
            string_set = True
            if string != str(self.string_index):
                self.string_index = int(string)
                self.repaint()
                
        if change:
            self.lattice = tfim.Lattice(self.L, self.PBC)
            self.N = self.lattice.N
            if self.draw_lattice == False:
                if not string_set:
                    self.string_index = 0
                self.bonds = bond_list(self)
                self.Jij = make_Jij(self)
                self.plaq = make_plaquettes(self)
                self.f_plaq = frustrated(self)
                self.node_pairs = plaq_pairing(self)
                
                init_ground = initial_ground(self)
                self.p_pairings = init_ground[0]
                self.ground_distance = init_ground[1]
                v_edges = viable_edges(self)
                self.plaq_pairs = v_edges[0]
                self.str_length = v_edges[1]
                self.edges = v_edges[2]
                self.matchings = plaq_groups(self)
                
                self.string_groups = add_all_strings(self)
                self.b_bonds = broken_bonds(self)
                true_ground = make_config(self)
                self.ground_config = true_ground[0]
                self.true_ground_strings = true_ground[1]
                self.number_ground_states = len(self.true_ground_strings)
                
                if self.number_ground_states == 0:
                    found_ground = False
                else:
                    found_ground = True
                    
                incremented = 0
                while found_ground == False:
                    incremented += 1
                    self.ground_distance += 1
                    v_edges = viable_edges(self)
                    self.plaq_pairs = v_edges[0]
                    self.str_length = v_edges[1]
                    self.edges = v_edges[2]
                    if len(self.edges) != 0:
                        self.matchings = plaq_groups(self)
                        self.string_groups = add_all_strings(self)
                        self.string_index = 0
                        self.b_bonds = broken_bonds(self)
                        true_ground = make_config(self)
                        self.ground_config = true_ground[0]
                        self.true_ground_strings = true_ground[1]
                        self.number_ground_states = len(self.true_ground_strings)
                        if self.number_ground_states != 0:
                            found_ground = True
                
                print('incremented: ', incremented)
                self.countLabel.setText("Number of ground states: " + str(self.number_ground_states))
                
            if self.spins == True:
                self.configuration = list(map(int,list(bin(self.cnfg)[2:].zfill(self.xwidth*self.yheight))))
                self.Clabel.setText('Configuration: ' + str(self.configuration))
        
        self.repaint()



    def paintEvent(self,event):
        qp = QPainter()
        qp.begin(self)
        qp.setPen(QPen(Qt.blue, 3))
        self.drawBonds(qp,self.coordList)
        self.drawSpins(qp,self.coordList)
        if self.spins == True:
            self.drawConfiguration(qp,self.configuration,self.coordList)
            self.drawBonds_spins(qp, self.coordList)
        self.draw_frustrated(qp)
        self.drawStrings(qp)
        qp.end()


    
    def drawSpins(self, qp, coordList):
        qp.setBrush(QBrush(QColor(0,0,0), Qt.SolidPattern))
        for i in range(len(coordList)):
            qp.setPen(QPen(Qt.black, 5))
            c = coordList[i]
            qp.drawEllipse(c[0]-6,c[1]-6,12,12)
            qp.setFont(QFont('Decorative', 9.5))
            qp.drawText(c[0]+5, c[1]+15, str(i))

        
    def drawBonds(self, qp, coordList):
        for i, coord in enumerate(coordList):
            NNs = self.lattice.NN(i)
            for j in NNs:
                if self.draw_lattice == False:
                    qp.setPen(QPen(bondColor(self.Jij,i,j), 5))
                else:
                    qp.setPen(QPen(Qt.black, 5))
                drawn = False
                coord2 = coordList[j]
                if j == (i + (self.yheight-1)*self.xwidth) or i == (j + (self.yheight-1)*self.xwidth):
                    if (self.yheight != 2 or self.PBC == True) and i < j:
                        drawn = True
                        qp.drawLine(coord2[0], coord2[1], coord2[0], coord2[1] + (self.spacing)/2)
                        if self.yheight == 2:
                            qp.drawLine(coord[0], coord[1], coord2[0], coord2[1])
                if j == i + self.xwidth - 1 or i == j + self.xwidth - 1:
                    if (self.xwidth != 2 or self.PBC == True) and i < j:
                        drawn = True
                        qp.drawLine(coord2[0], coord2[1], coord2[0] + (self.spacing)/2, coord2[1])
                        if self.xwidth == 2:
                            qp.drawLine(coord[0], coord[1], coord2[0], coord2[1])
                if drawn == False:
                    if i<j:
                        qp.drawLine(coord[0], coord[1], coord2[0], coord2[1])
                        
    def drawBonds_spins(self, qp, coordList):
        for i in range(self.N):
            for j in range(i+1,self.N):
                if int(self.Jij[i][j]*(2*self.configuration[-i-1]-1)*(2*self.configuration[-j-1]-1)) == -1: #bond broken
                    qp.setPen(QPen(Qt.white, 2, Qt.DashLine))
                    diff = False
                    #ci = coordList[i]
                    cj = coordList[j]
                    if self.PBC:
                        if j == i+self.xwidth-1:
                            diff = True
                            ci = [cj[0]+23.5, cj[1]] #to the right
                        elif j == i + (self.yheight-1)*self.xwidth: #down from bottom
                            diff = True
                            ci = [cj[0], cj[1]+23.5]
                        else:
                            ci = coordList[i]
                    else:
                        ci = coordList[i]
                    if ci[0] == cj[0]:
                        if not diff:
                            qp.drawLine(ci[0],ci[1]+8.5,cj[0],cj[1]-8.5)
                        else:
                            qp.drawLine(ci[0],ci[1]+8.5,cj[0],cj[1]+8.5) #drawing line down
                    else:
                        if not diff:
                            qp.drawLine(ci[0]+8.5, ci[1], cj[0]-8.5, cj[1])
                        else:
                            qp.drawLine(ci[0]+8.5,ci[1],cj[0]+8.5,cj[1])
            

    def drawConfiguration(self,qp,config,coordList):
        i=0
        for c in coordList:
            if config[-i - 1] == 1:
                qp.setPen(QPen(Qt.magenta, 2))
                self.drawArrow(qp,1,c)
            else:
                qp.setPen(QPen(Qt.cyan, 2))
                self.drawArrow(qp,-1,c)
            i+=1
            
    def draw_frustrated(self,qp):
        if self.draw_lattice == False:
            for plaq in self.f_plaq:
                qp.setPen(QPen(Qt.magenta, 4))
                coord = self.coordList[plaq]
                qp.drawPoint(coord[0] + (self.spacing)/2, coord[1] + (self.spacing)/2)
                
            if self.PBC:
                num_of_bonds = 2*self.N
            else:
                num_of_bonds = (self.xwidth - 1)*(self.yheight) + (self.xwidth)*(self.yheight - 1)
                
            #self.scountLabel.setText("Unsatisfied: " + str(self.str_length))
        
    def drawArrow(self, qp, dir, coords):
        qp.drawLine(coords[0],coords[1]+(dir*6),coords[0],coords[1]-(dir*6))
        qp.drawLine(coords[0]-6,coords[1]-(dir*2.2),coords[0],coords[1]-(dir*6))
        qp.drawLine(coords[0]+6,coords[1]-(dir*2.2),coords[0],coords[1]-(dir*6))
        
    def drawStrings(self,qp):
        strings = self.true_ground_strings[self.string_index]
        qp.setPen(QPen(Qt.magenta, 4))
        for string in strings:
            for index in range(0, len(string)-1):
                switched = False
                pp1 = string[index]
                pp2 = string[index + 1]
                if pp2 < pp1:
                    p1 = pp2
                    p2 = pp1
                else:
                    p1 = pp1
                    p2 = pp2
                coord1x = self.coordList[p1][0] + (self.spacing)/2
                coord1y = self.coordList[p1][1] + (self.spacing)/2
                coord2x = self.coordList[p2][0] + (self.spacing)/2
                coord2y = self.coordList[p2][1] + (self.spacing)/2
                xweight = abs(coord2x - coord1x)/self.spacing
                yweight = abs(coord1y - coord2y)/self.spacing
                
                qp.setPen(QPen(Qt.magenta, 4))
                if xweight > 1:
                    qp.drawLine(coord1x, coord1y, coord1x - self.spacing, coord1y)
                    qp.drawPoint(coord2x, coord2y)
                elif yweight > 1:
                    qp.drawLine(coord1x, coord1y, coord1x, coord1y - self.spacing)
                    qp.drawPoint(coord2x, coord2y)
                else:
                    qp.drawLine(coord1x, coord1y, coord2x, coord2y)
                
            qp.setPen(QPen(Qt.black, 2.5))
            qp.drawPoint(self.coordList[string[0]][0] + (self.spacing)/2, self.coordList[string[0]][1] + (self.spacing)/2)
            qp.drawPoint(self.coordList[string[-1]][0] + (self.spacing)/2, self.coordList[string[-1]][1] + (self.spacing)/2)
                
            '''if index == 0:
                    if switched:
                        qp.drawPoint(coord2x, coord2y)
                    else:
                        qp.drawPoint(coord1x, coord1y)
                if index == len(string) - 2:
                    if switched:
                        qp.drawPoint(coord1x, coord1y)
                    else:
                        qp.drawPoint(coord2x, coord2y)'''
        return
        
        

def spinCoords(self):
    coords = []
    y_init = self.center[1] - ((self.yheight*self.spacing)/2)
    for j in range(0, self.yheight):
        y = y_init + (j * self.spacing)
        x_init = self.center[0] - ((self.xwidth*self.spacing)/2)
        for i in range(0, self.xwidth):
            x = x_init + (i * self.spacing)
            c = (x,y)
            coords.append(c)
    return coords


def bond_list(self):
    np.random.seed(self.seed)
    # Generates a random list of bonds with equal numbers of ferromagnetic and antiferromagnetic bonds
    N = self.N
    if self.PBC == True:
        num_of_bonds = 2*N
    else:
        num_of_bonds = (self.xwidth - 1)*(self.yheight) + (self.xwidth)*(self.yheight - 1)
    if num_of_bonds%2 == 0:
        a1 = [-1 for i in range(num_of_bonds//2)]
    else:
        a1 = [-1 for i in range((num_of_bonds//2) + 1)]
    a2 = [1 for i in range(num_of_bonds//2)]
    a = list(np.random.permutation(a1+a2))
    return a
    


def make_Jij(self):
    bond_index = 0
    N = self.N
    b_list = self.bonds
    Jij = np.zeros((N,N))
    for i in range(0,N):
        lattice = self.lattice
        NNs = self.lattice.NN(i)
        for j in NNs:
            if Jij[i][j] == 0:
                Jij[i][j] = b_list[bond_index]
                Jij[j][i] = b_list[bond_index]
                bond_index += 1
    return Jij
    
    
def make_plaquettes(self):
    p_list = []
    if self.PBC:
        for i in range(0, self.N):
            NNs = self.lattice.NN(i)
            plaq = [i]
            plaq.append(NNs[3])
            NNs2 = self.lattice.NN(NNs[3])
            plaq.append(NNs2[1])
            plaq.append(NNs[1])
            p_list.append(plaq)
    else:
        for y in range(0,self.yheight):
            for x in range(0, self.xwidth):
                if y == self.yheight-1 or x == self.xwidth-1: #This part adds empty plaquettes so the first number is also the index of each plaquette
                    #if i want to take it out, just need to subtract 1 from x and y range
                    p_list.append([])
                else:
                    plaq = []
                    i = y*self.xwidth + x
                    plaq.append(i)
                    plaq.append(i+1)
                    plaq.append(i+self.xwidth+1)
                    plaq.append(i+self.xwidth)
                    p_list.append(plaq)
    return p_list
    
    
def frustrated(self):
    Jij = self.Jij
    f_plaq = []
    for plaq in self.plaq:
        count = 0
        if len(plaq)!=0:
            if Jij[plaq[0]][plaq[1]] == -1:
                count += 1
            if Jij[plaq[1]][plaq[2]] == -1:
                count += 1
            if Jij[plaq[2]][plaq[3]] == -1:
                count += 1
            if Jij[plaq[0]][plaq[3]] == -1:
                count += 1
            if count == 1 or count == 3:
                f_plaq.append(plaq[0])
    return f_plaq
        
        
def plaq_pairing(self):
    '''Function returns a list of all possible pairs between frustrated plaquettes with the distances between them–The distance is stored as the maximum possible distance between two plaquettes minus the actual distance'''
    pair_list = []
    for index, p1 in enumerate(self.f_plaq):
        coord1 = self.coordList[p1]
        for p2 in self.f_plaq[index+1:]:
            coord2 = self.coordList[p2]
            x1 = coord1[0] + (self.spacing)/2
            x2 = coord2[0] + (self.spacing)/2
            
            y1 = coord1[1] + (self.spacing)/2
            y2 = coord2[1] + (self.spacing)/2
            
            xdiff = abs((x1 - x2)/self.spacing)
            ydiff = abs((y2-y1) / self.spacing)
            if self.PBC:
                if xdiff > (self.xwidth)//2:
                    xdiff = (self.xwidth) - xdiff
                if ydiff > (self.yheight)//2:
                    ydiff = (self.yheight) - ydiff
            else:
                if xdiff > (self.xwidth-1)//2:
                    xdiff = (self.xwidth - 1) - xdiff
                if ydiff > (self.yheight-1)//2:
                    ydiff = (self.yheight) - ydiff
            tot_dist = int(xdiff + ydiff)
            #Here we build a list of pairs with the distance between them
            max = self.xwidth//2 + self.yheight//2
            op_dist = (max - tot_dist)
            if p1 > p2:
                pair_list.append((p2, p1, op_dist))
            else:
                pair_list.append((p1, p2, op_dist)) #Pair_list does not have true distance between pairs, the true distance is subtracted from the max possible distance to help with node matching later
    return pair_list
    

def initial_ground(self):

    pair_list = self.node_pairs
    G = nx.Graph()
    G.add_weighted_edges_from(pair_list) #makes graph of all node pairs
    matching = nx.max_weight_matching(G) #gives one solution
    ground_dist = 0
    p_pairs = []
    for pair in matching:
        edge = G.get_edge_data(pair[0], pair[1])
        pair_dist = (self.xwidth//2 + self.yheight//2)-edge['weight']
        ground_dist += pair_dist #total string length for a ground state soln
        if pair[0] > pair[1]:
            p0 = pair[1]
            p1 = pair[0]
            pair = (p0, p1)
        p_pairs.append([pair, pair_dist]) #adds solution from above to list with pairs and pair distance
    return p_pairs, ground_dist


  
def viable_edges(self):
    '''Function takes the list of all possible pairings of nodes and returns a list of all combinations of those nodes that would result in a ground state energy'''
    pair_list = self.node_pairs # The list of all pairs of nodes
    p_pairs = self.p_pairings
    ground_dist = self.ground_distance
    edge_lst = []
    plaq_dict = {}
    
    #Make the list for the edges grouped by plaquette
    for index, plaq in enumerate(self.f_plaq):
        edge_lst.append([])
        plaq_dict[plaq] = index
    
    G = nx.Graph()
    G.add_weighted_edges_from(pair_list)
   
    first = False
    p_dist = 0
    for pair in p_pairs:
        dist = pair[1]
        p_dist += dist
    if p_dist == ground_dist:
        first = True #So we only remove certain edges if p_pairs is relevant
    
    if first:
        for pair in p_pairs:
            plaq = pair[0][0]
            ind = plaq_dict.get(plaq)
            edge_lst[ind].append(pair)
        
    for plaq in self.f_plaq:
        G2 = G.copy()
        if first:
            for pair in p_pairs:
                if pair[0][1] == plaq or pair[0][0] == plaq:
                    G2.remove_edge(*pair[0])
                    break
                
        ground_energy = True
        while ground_energy == True:
            matching = nx.max_weight_matching(G2)
            if len(matching) != len(self.f_plaq)/2: #This would happen if we have taken out all edges for a particular plaquette
                ground_energy = False
                break
            new_length = 0
            new_group = []
            for pair in matching:
                edge = G2.get_edge_data(pair[0], pair[1])
                if pair[0] == plaq or pair[1] == plaq:
                    rem_edge = (pair[0], pair[1])
                pair_dist = (self.xwidth//2 + self.yheight//2)-edge['weight']
                new_length += pair_dist
                if pair[0] > pair[1]:
                    p0 = pair[1]
                    p1 = pair[0]
                    pair = (p0, p1)
                new_group.append([pair, pair_dist])
            if new_length == ground_dist:
                G2.remove_edge(*rem_edge)
                for pair in new_group:
                    plaq = pair [0][0]
                    ind = plaq_dict.get(plaq)
                    
                    if pair not in edge_lst[ind]:
                        edge_lst[ind].append(pair)
            elif new_length < ground_dist:
                G2.remove_edge(*rem_edge)
            else:
                ground_energy = False
                ind = plaq_dict[plaq]
    zeroes = True
    for plaq in edge_lst:
        if len(plaq) != 0:
            zeroes = False
            break
    if zeroes:
        edge_lst = []
    return p_pairs, ground_dist, edge_lst



def plaq_groups(self):
    edges = self.edges
    
    group = []
    used_plaquettes = []
    all_groups = []
    current_plaq = 0
    p_ind = 0
    index = 0
    loop_count = 0
    
    new = False
    running = True
    
    plaq_dict = {}
    for index, plaq in enumerate(self.f_plaq): #Allows me to find the index of the frustrated plaquettes
        plaq_dict[plaq] = index
    
    #The main piece of the function
    while running:
        for group_index, p_edges in enumerate(edges[current_plaq:]):
            if new:
                new = False #new allows me to restart the for loop when i change current_plaq
                break
                
            if group_index + current_plaq == len(edges) - 1: #if we get to the last group without having a full ground state, we need to use a different combination of edges
                try_new = False
                for_loop = False
                for e_ind, edge in enumerate(group[::-1]): #going through the ground state group backwards to see if other edge choices will work
                    loop_count +=1
                    if loop_count > 1000000:
                        running = False
                        new = True
                    for_loop = True
                    if try_new == True:
                        break
                    else:
                        plaq_ind = plaq_dict.get(edge[0][0])
                        for e_index, edge2 in enumerate(edges[plaq_ind]):
                            if edge2 == edge:
                                if e_index == len(edges[plaq_ind])-1 and plaq_ind == 0: #end of program, we've reached the last entry of the first plaquette
                                    running = False
                                    new = True
                                    try_new = True
                                    break
                                elif e_index == len(edges[plaq_ind])-1: #Move to the previous plaquette list to find a viable edge
                                    break
                                else:
                                    current_plaq = plaq_ind #move to the next edge in the list for the plaquette
                                    p_ind = e_index + 1
                                    try_new = True
                                    new = True
                                    break
                if for_loop:
                    group = group[:-e_ind]
                    used_plaquettes = used_plaquettes[:-e_ind]
                    if len(group) == 1 and group[0][0][0] == plaq_dict.get(current_plaq):
                        group = []
                        used_plaquettes = []
            for pair in p_edges[p_ind:]: #This is a single edge from edge list
                p_ind = 0 #Resets p_ind for the next loop through
                if pair[0][0] in used_plaquettes:
                    break #Can move to next plaquette
                elif pair[0][1] in used_plaquettes:
                    continue #Need to go through to the next pair
                else:
                    group.append(pair)
                    used_plaquettes.append(pair[0][1])#maybe add the first element of pair if this is the first plaquette we visit
                    if len(group) == len(self.f_plaq)//2: #Group is full
                    
                        all_groups.append(group)
                        last_pair = group[-2] #This is the pair that we remove and replace before cycling through other options
                        ind = plaq_dict.get(last_pair[0][0]) #The plaquette index
                        group = group[:-2]
                        used_plaquettes = used_plaquettes[:-2]
                        found = False
                        while found == False:
                            loop_count += 1
                            if loop_count > 1000000:
                                running = False
                                new = True
                                found = True
                                break
                            for index, pairing in enumerate(edges[ind]):
                                if pairing == last_pair and index == len(edges[ind])-1: #This happens if we are at the last pair of a particular plaquette
                                    if len(group) == 0: #This happens if we have gotten through the last edge of the first plaquette, function is done
                                        running = False
                                        found = True
                                        break
                                    last_pair = group[-1] #Take off the last pair and go to that plaquette to see if there are further pairs to use
                                    ind = plaq_dict.get(last_pair[0][0])
                                    group = group[:-1]
                                    used_plaquettes = used_plaquettes[:-1]
                                elif pairing == last_pair: #This means there are more pairs for the plaquette in question, so we adjust current_plaq and p_ind, and go through the for loops again from there
                                    current_plaq = ind
                                    p_ind = index + 1
                                    found = True
                                    new = True
                                    break
                    break
    return all_groups
    

def add_all_strings(self):
    edges = []
    groups = self.matchings
    for i in range(len(self.coordList)):
        NNs = self.lattice.NN(i)
        for j in NNs:
            if i < j:
                edges.append((i,j))
    G = nx.Graph()
    G.add_edges_from(edges) #G has edges connecting all points in a lattice with PBC
    
    all_groups = []
    index = 0
    for group in groups:
    
        if index > 20000:
            print('Not all ground states found')
            break
    
        single_pairing = []
        index += 1
        for pairing in group:
            paths = nx.all_shortest_paths(G, pairing[0][0], pairing[0][1]) #finds all possible paths between two points
            paths_list = []
            for path in paths:
                paths_list.append(path)
            single_pairing.append(paths_list)
        path_combos = it.product(*single_pairing)
        for combo in path_combos:
            all_groups.append(combo)
    return all_groups
    

def broken_bonds(self):
    '''Returns an NxN matrix with 1's where there are broken bonds between two spins'''
    config_Jij_list = []
    for str_index, state in enumerate(self.string_groups):
        config_Jij = np.zeros((self.N,self.N))
        for string in state:
            for index in range(0, len(string)-1):
                p1 = string[index]
                p2 = string[index + 1]
                if p1>p2:
                    hold = p1
                    p1 = p2
                    p2 = hold
                c1x = self.coordList[p1][0]
                c2x = self.coordList[p2][0]
                if c1x == c2x:
                    if p2 + self.xwidth > self.N - 1 and p1 < self.xwidth:
                        sp1 = p1
                        if (p1+1) % self.xwidth == 0:
                            sp2 = p1 - self.xwidth + 1
                        else:
                            sp2 = p1 + 1
                    else:
                        sp1 = p2
                        if (p1+1) % self.xwidth == 0: #on the far right
                            sp2 = p2 - self.xwidth + 1
                        else:
                            sp2 = p2 + 1
                else:
                    if p2 + self.xwidth > self.N - 1:
                        if (p2+1) % self.xwidth == 0:
                            if p1 % self.xwidth == 0:
                                sp1 = p1
                        else:
                            sp1 = p2
                        sp2 = sp1 - (self.xwidth * (self.yheight - 1))
                    elif (p2+1) % self.xwidth == 0:
                        if p1 % self.xwidth == 0:
                            sp1 = p1
                        else:
                            sp1 = p2
                        sp2 = sp1 + self.xwidth
                    else:
                        sp1 = p2
                        sp2 = p2 + self.xwidth
                bond = (sp1, sp2)
                config_Jij[sp1][sp2] = 1
                config_Jij[sp2][sp1] = 1
        config_Jij_list.append([config_Jij, str_index])
    return config_Jij_list
    
def make_config(self):
    ground_states = []
    true_strings = []
    s_index = 0
    #print('initial number potential ground states: ', len(self.string_groups))
    for Jij in self.b_bonds:
        broken = Jij[0]
        bonds = self.Jij
        spin_list = []
        spin_list.append(0) #Set the first spin as down
        valid = True
        
        #Loop through all other spins
        for sp1 in range(1, self.N):
            if valid == False:
                #print('String group: ', self.string_groups[Jij[1]])
                break
            if sp1 % self.xwidth == 0:
                sp2 = sp1 - self.xwidth
            else:
                sp2 = sp1 - 1
            spin2 = spin_list[sp2]
            bond = bonds[sp1][sp2]
            status = broken[sp1][sp2]
            
            #Set spin
            if bond == 1:
                #Spins want to be the same
                if status == 1: #broken
                    spin1 = abs(spin2 - 1)
                else:
                    spin1 = spin2
            else:
                #Spins want to be opposite
                if status == 1: #broken
                    spin1 = spin2
                else:
                    spin1 = abs(spin2 - 1)
            spin_list.append(spin1)
            
            #Check bonds to all lower spins
            NNs = self.lattice.NN(sp1)
            for i in NNs:
                if i < sp1:
                    spini = spin_list[i]
                    bond = bonds[sp1][i]
                    status = broken[sp1][i]
                    if bond == 1: #Spins want to be same
                        if status == 1: #Spins should be opposite
                            if spin1 == spini:
                                valid = False
                                #print("1")
                                break
                        else: #Spins should be same
                            if spin1 != spini:
                                valid = False
                                #print("2")
                                break
                    else: #Spins want to be opposite
                        if status == 1: #spins should be same
                            if spini != spin1:
                                valid = False
                                #print("3")
                                break
                        else: #Spins should be opposite
                            if spini == spin1:
                                valid = False
                                #print("4")
                                break
        if valid:
            index = 0
            for i in range(0, self.N):
                if spin_list[i] == 1:
                    index += 2**i
            #spin_list.reverse()
            ground_states.append([index, s_index])
            s_index +=1
            true_strings.append(self.string_groups[Jij[1]])
            
    #print('ground states: ', ground_states)
    return ground_states, true_strings
    
def bondColor(Jij,i,j):
    if Jij[i][j] == 1:
        return Qt.red
    else:
        return Qt.blue


if __name__ == "__main__":
    App = QApplication(sys.argv)
    window = Window()
    sys.exit(App.exec())
    
