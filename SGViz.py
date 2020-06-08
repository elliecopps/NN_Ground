import numpy as np
import random
import tfim
import itertools as it

import threading
import sys, os
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

		#Set these to True to see bipartitions or Energy info, respectively
		######################################################################
		self.show_bps = False
		self.show_Ens = True
		######################################################################

		#Window size and label
		self.title = "SGViz_InfD"
		self.top= 50
		self.left= 100
		self.width = 1000
		self.height = 750


		# Set initial parameters
		self.N = 8
		self.seed = 0
		self.G = tfim.Jij_instance(self.N,1,"bimodal",self.seed,True)
		self.Jij = makeJij(self.G, self.N)
		self.cnfg = 0
		self.configuration = list(map(int,list(bin(self.cnfg)[2:].zfill(self.N))))
		self.center = (500,375)
		self.coordList = spinCoords(self.N, self.center,290)
		combos = list(it.combinations([i for i in range(int(self.N))],int(self.N / 2)))
		self.permlist = []
		for i in range(int(len(combos)/2)):
			self.permlist.append(combos[i] + combos[int(len(combos))-1-i])
		self.pi = 0
		self.p = self.permlist[self.pi]

		if self.show_Ens:
			self.lattice = tfim.Lattice([self.N])
			self.basis = tfim.IsingBasis(self.lattice)
			self.JZZ = tfim.JZZ_SK(self.basis,self.G)
			self.ea = -(self.JZZ.diagonal())
			self.gs = np.where(self.ea == self.ea.min())[0]



		self.setStyleSheet("background-color: white;")

		if self.show_Ens:
			self.Ealabel = QLabel(self)
			self.Ealabel.setFont(QFont('Decorative', 13))
			self.Ealabel.setStyleSheet("color: black")
			self.Ealabel.resize(300,30)
			self.Ealabel.move(20,50)
			self.Ealabel.setText('GS Configurations: ' + str(self.gs))

			self.Enlabel = QLabel(self)
			self.Enlabel.setFont(QFont('Decorative', 13))
			self.Enlabel.setStyleSheet("color: black")
			self.Enlabel.resize(200,30)
			self.Enlabel.move(20,80)
			self.Enlabel.setText('Current Energy: ' + str(self.ea[self.cnfg]))

		self.Clabel = QLabel(self)
		self.Clabel.setFont(QFont('Decorative', 13))
		self.Clabel.setStyleSheet("color: black")
		self.Clabel.resize(400,30)
		self.Clabel.move(700,20)
		self.Clabel.setText('Configuration: ' + str(self.configuration))
		
		self.cEdit = QLineEdit(self)
		self.cEdit.setStyleSheet("color: black")
		self.cEdit.move(80, 700)
		self.cEdit.resize(100,32)
		self.cLabel = QLabel(self)
		self.cLabel.setStyleSheet("color: black")
		self.cLabel.setText("Config")
		self.cLabel.move(20,700)
		self.cLabel.resize(50,32)

		self.countLabel = QLabel(self)
		self.countLabel.setStyleSheet("color: black")
		self.countLabel.setText("unsatisfied: ")
		self.countLabel.move(20,200)
		self.countLabel.resize(100,32)


		self.scountLabel = QLabel(self)
		self.scountLabel.setStyleSheet("color: black")
		self.scountLabel.setText("satisfied: ")
		self.scountLabel.move(20,250)
		self.scountLabel.resize(100,32)

		self.sEdit = QLineEdit(self)
		self.sEdit.setStyleSheet("color: black")
		self.sEdit.move(80, 650)
		self.sEdit.resize(100,32)
		self.sLabel = QLabel(self)
		self.sLabel.setStyleSheet("color: black")
		self.sLabel.setText("Seed")
		self.sLabel.move(20,650)
		self.sLabel.resize(50,32)

		self.pEdit = QLineEdit(self)
		self.pEdit.setStyleSheet("color: black")
		self.pEdit.move(850, 650)
		self.pEdit.resize(100,32)
		self.pLabel = QLabel(self)
		self.pLabel.setStyleSheet("color: black")
		self.pLabel.setText("Bipartition")
		self.pLabel.move(740,650)
		self.pLabel.resize(75,32)


		self.NEdit = QLineEdit(self)
		self.NEdit.setStyleSheet("color: black")
		self.NEdit.move(80, 600)
		self.NEdit.resize(100,32)
		self.NLabel = QLabel(self)
		self.NLabel.setStyleSheet("color: black")
		self.NLabel.setText("N")
		self.NLabel.move(20,600)
		self.NLabel.resize(50,32)

		self.pybutton = QPushButton('Enter', self)
		self.pybutton.setStyleSheet("QPushButton {color: white; background-color: black}")
		self.pybutton.clicked.connect(self.clickMethod)
		self.pybutton.resize(100,32)
		self.pybutton.move(300, 700)

		self.clear = False
		self.Cpybutton = QPushButton('Clear', self)
		self.Cpybutton.setStyleSheet("QPushButton {color: white; background-color: black}")
		self.Cpybutton.clicked.connect(self.clearMethod)
		self.Cpybutton.resize(100,32)
		self.Cpybutton.move(10, 10)



		
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
			if self.show_Ens:
				self.Ealabel.move(1000,50)
				self.Enlabel.move(1000,80)
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
			if self.show_Ens:
				self.Ealabel.move(20,50)
				self.Enlabel.move(20,80)
			self.clear = False




	def DynMethod(self):

		change = False

		ns = self.NEdit.text()	
		if len(ns) > 0:
			change = True
			n = int(ns)
			if n != self.N:
				self.N = n
				self.coordList = spinCoords(self.N, (500,375),300)
				lattice = tfim.Lattice([self.N])
				basis = tfim.IsingBasis(lattice)
				combos = list(it.combinations([i for i in range(int(self.N))],int(self.N / 2)))
				self.permlist = []
				for i in range(int(len(combos)/2)):
					self.permlist.append(combos[i] + combos[int(len(combos))-1-i])
				self.pi = 0
				self.p = self.permlist[self.pi]

		ss = self.sEdit.text()
		if len(ss) > 0:
			change = True
			s = int(ss)
			if s != self.seed:
				self.seed = s
				ss = self.sEdit.text()

		cs = self.cEdit.text()
		if len(cs) > 0:
			change = True
			c = int(cs)
			if c != self.cnfg:
				self.cnfg = c

		if change:
			self.G = tfim.Jij_instance(self.N,1,"bimodal",self.seed,True)
			self.Jij = makeJij(self.G, self.N)
			self.configuration = list(map(int,list(bin(self.cnfg)[2:].zfill(self.N))))
			self.Clabel.setText('Configuration: ' + str(self.configuration))


			if self.show_Ens:
				self.lattice = tfim.Lattice([self.N])
				self.basis = tfim.IsingBasis(self.lattice)
				self.JZZ = tfim.JZZ_SK(self.basis,self.G)
				self.ea = -(self.JZZ.diagonal())
				self.gs = np.where(self.ea == self.ea.min())[0]
				self.Ealabel.setText('GS Configurations: ' + str(self.gs))
				self.Enlabel.setText('Current Energy: ' + str(self.ea[self.cnfg]))






		ps = self.pEdit.text()
		if len(ps) > 0:
			bp = int(ps)
			if bp != self.pi:
				self.pi = bp
				self.p = self.permlist[self.pi]

		self.repaint()




	def paintEvent(self,event):
		qp = QPainter()
		qp.begin(self)
		qp.setPen(QPen(Qt.blue, 3))
		self.drawBonds(qp,self.coordList,self.configuration)
		self.drawSpins(qp,self.coordList)
		self.labelSpins(qp,self.coordList)
		self.drawConfiguration(qp,self.configuration,self.coordList)
		qp.end()



	def drawSpins(self, qp, coordList):
		qp.setBrush(QBrush(QColor(0,0,0), Qt.SolidPattern))
		for i in range(len(coordList)-1,-1,-1):
			qp.setPen(QPen(spinColor(self.p,i,self.show_bps), 5))

			c = coordList[i]
			qp.drawEllipse(c[0]-25,c[1]-25,50,50)

	def drawBonds(self, qp, coordList,config):
		self.ubondCountArray = [0 for i in range(self.N)]
		self.sbondCountArray = [0 for i in range(self.N)]
		for i in range(self.N):
			for j in range(i+1,self.N):
				if int(self.Jij[i][j]*(2*config[i]-1)*(2*config[j]-1)) == -1:
					self.ubondCountArray[i]+=1
					self.ubondCountArray[j]+=1
					qp.setPen(QPen(bondColor(self.Jij,i,j), 5, Qt.DashLine))
				else:
					self.sbondCountArray[i]+=1
					self.sbondCountArray[j]+=1
					qp.setPen(QPen(bondColor(self.Jij,i,j), 5, Qt.SolidLine))
				ci = coordList[i]
				cj = coordList[j]
				qp.drawLine(ci[0],ci[1],cj[0],cj[1])
		self.countLabel.setText("Unsatisfied: " + str(sum(self.ubondCountArray)))
		self.scountLabel.setText("Satisfied: " + str(sum(self.sbondCountArray)))

	def labelSpins(self, qp, coordList):
		qp.setPen(Qt.black)
		
		i=0
		for c in coordList:
			angle = (2*i*(np.pi))/self.N 
			qp.setFont(QFont('Decorative', 15))
			#uncomment following line for labels
			qp.drawText(self.center[0] -5- 350*np.cos(angle),self.center[1] +5- 350*np.sin(angle), str(i))

			i+=1 

	def drawConfiguration(self,qp,config,coordList):
		i=0
		for c in coordList:
			if config[i] == 1:
				qp.setPen(QPen(Qt.magenta, 4))
				self.drawArrow(qp,1,c)
			else:
				qp.setPen(QPen(Qt.cyan, 4))
				self.drawArrow(qp,-1,c)
			i+=1
			

	def drawArrow(self, qp, dir, coords):
		qp.drawLine(coords[0],coords[1]+(dir*15),coords[0],coords[1]-(dir*15))
		qp.drawLine(coords[0]-10,coords[1]-(dir*5),coords[0],coords[1]-(dir*15))
		qp.drawLine(coords[0]+10,coords[1]-(dir*5),coords[0],coords[1]-(dir*15))





        


def spinCoords(N,center,r):
	coords = []
	for i in range(N):
		angle = (2*i*(np.pi))/N 
		x = center[0] - r*np.cos(angle)
		y = center[1] - r*np.sin(angle)
		c = (x,y)
		coords.append(c)
	return coords


def makeJij(G,N):
	"""
	Turns Jij matrix from form built in tfim.py to standard Jij where J[i][j] is the bond between spins i and j
	"""
	Jij = np.zeros((N,N))
	for j in range(N//2):
		for i in range(N):
			Jij[i][(i-j+N-1) % N] = Jij[(i-j+N-1) % N][i] = G[j][i]
	return Jij


def bondColor(Jij,i,j):
	if Jij[i][j] == 1:
		return Qt.red
	else:
		return Qt.blue

#BIPARTITION COLOR

def spinColor(p,i,show_bps):
	if show_bps:
		A = p[0:len(p)//2]
		if i in A:
			return Qt.yellow
		else:
			return Qt.green
	else:
		return Qt.black


if __name__ == "__main__":
	App = QApplication(sys.argv)
	window = Window()
	sys.exit(App.exec())