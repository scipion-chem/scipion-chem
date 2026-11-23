#Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# # -*- coding: utf-8 -*-
# # # **************************************************************************
# # # *
# # # * Authors: Daniel Del Hoyo Gómez (ddelhoyo@cnb.csic.es)
# # # *
# # # *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************
import os

import math, sys
import numpy as np

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Patch
import matplotlib.image as mpimg
from matplotlib.lines import Line2D
from matplotlib.figure import Figure

from rdkit import Chem
from rdkit.Chem import rdDepictor, AllChem
from rdkit.Geometry import Point3D
from rdkit.Chem.Draw import rdMolDraw2D

from utils import parseParams, getMolFilesDic

RES_COLORS = {'Apolar': 'green', 'Negative': 'blue', 'Positive': 'red', 'Polar': 'cyan', 'Aromatic': 'purple'}
RES_GROUPS = {'Apolar': ['GLY', 'ALA', 'VAL', 'LEU', 'ILE', 'MET'], 'Aromatic': ['PHE', 'TYR', 'TRP'],
							'Polar': ['SER', 'THR', 'CYS', 'PRO', 'ASN', 'GLN'], 'Negative': ['ASP', 'GLU'],
							'Positive': ['LYS', 'ARG', 'HIS']}

COLOR_DIC = {res: RES_COLORS[group] for group, resGroup in RES_GROUPS.items() for res in resGroup}

# Hbond and waterbridge are both blue, told apart by how dark they are rather than by their hue: the
# lines are drawn half transparent, which washes a light colour out
COLOR_INT_DIC = {'halogen': 'red', 'hbond': 'royalblue', 'hydrophobic': 'green', 'pication': 'orange',
								 'pistacking': 'brown', 'saltbridge': 'gray', 'waterbridge': 'navy'}

# Unit of the number written on each interaction line: a distance for a single pose, the
# fraction of frames it is present in for a trajectory
DEF_VALUE_UNIT = ''
# Free space left around what is drawn, in the units of the drawing, so there is room to move the residues around
VIEW_MARGIN = 300


class ResidueBox(FancyBboxPatch):
	"""Rounded rectangle wide enough to hold the name of a residue, drawn where a circle used to be.
	It keeps the centre as the circles had it, which is what the dragging and the lines work with"""
	def __init__(self, center, width, height, **kwargs):
		self.size = (width, height)
		FancyBboxPatch.__init__(self, (center[0] - width / 2, center[1] - height / 2), width, height,
														boxstyle='round,pad=0,rounding_size={}'.format(height / 2), **kwargs)

	def setSize(self, width, height):
		center = self.get_center()
		self.size = (width, height)
		self.set_bounds(center[0] - width / 2, center[1] - height / 2, width, height)
		self.set_boxstyle('round', pad=0, rounding_size=height / 2)

	def get_center(self):
		return self.get_x() + self.size[0] / 2, self.get_y() + self.size[1] / 2

	@property
	def center(self):
		return self.get_center()

	@center.setter
	def center(self, coords):
		self.set_x(coords[0] - self.size[0] / 2)
		self.set_y(coords[1] - self.size[1] / 2)

	def getBorderDistance(self, dx, dy):
		"""How far its border is from its centre in a direction, so a line can stop there"""
		length = np.hypot(dx, dy)
		if not length:
			return max(self.size) / 2

		width, height = self.size
		distances = [size / 2 / abs(d / length) for size, d in ((width, dx), (height, dy)) if d]
		return min(distances)


def getExternalPosition(atomPoint, molCenter, distance=80):
	"""Returns a position afar from the molecule to draw a residue circle"""
	dx = atomPoint[0] - molCenter[0]
	dy = atomPoint[1] - molCenter[1]

	length = math.sqrt(dx ** 2 + dy ** 2)
	if length > 0:
		dx /= length
		dy /= length

	newX = atomPoint[0] + dx * distance
	newY = atomPoint[1] + dy * distance

	return newX, newY


def orientHorizontally(mol):
	"""Rotate a 2D depiction so the molecule lies along the horizontal axis, which makes the most of
	the width of the figure and leaves the room above and below it for the residues"""
	conf = mol.GetConformer()
	coords = np.array([[conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y]
										 for i in range(mol.GetNumAtoms())])
	center = coords.mean(axis=0)
	# First principal component of the atom positions: the direction the molecule is longest in
	components = np.linalg.svd(coords - center, full_matrices=False)[2]
	angle = -math.atan2(components[0][1], components[0][0])

	rotation = np.array([[math.cos(angle), math.sin(angle)], [-math.sin(angle), math.cos(angle)]])
	for i, (x, y) in enumerate((coords - center) @ rotation + center):
		conf.SetAtomPosition(i, Point3D(float(x), float(y), 0.0))


def getMoleculeCenter(drawer, mol, canvasSize):
	"""Get the center of the drawn molecules"""
	coordsX, coordsY = [], []

	for i in range(mol.GetNumAtoms()):
		try:
			atomPoint = getDrawCoords(drawer, canvasSize, i)
			coordsX.append(atomPoint.x)
			coordsY.append(atomPoint.y)
		except RuntimeError:
			continue

	if coordsX and coordsY:
		centroX, centroY = np.mean(coordsX), np.mean(coordsY)
		return centroX, centroY
	else:
		return None

def getDrawCoords(drawer, canvasSize, i):
	atomPoint = drawer.GetDrawCoords(i)
	atomPoint.y = int(canvasSize * 0.5) - atomPoint.y
	return atomPoint

def setConnectionLineCoords(c1, c2, displac, margin=50, spacing=45):
	'''Draw the connection lines based on the extreme coordinates, with some margin.
	Also, the lines are displaced perpendicularly to themselves, by a fixed spacing, so several
	interactions of the same residue stay apart whatever direction they run in'''
	x1, x2 = c1[0], c2[0]
	y1, y2 = c1[1], c2[1]
	dx, dy = x2 - x1, y2 - y1
	length = np.hypot(dx, dy)

	ux, uy = dx / length, dy / length
	offX, offY = -uy * displac * spacing, ux * displac * spacing
	x1New, x2New = x1 + ux * margin/2 + offX, x2 - ux * margin + offX
	y1New, y2New = y1 + uy * margin/2 + offY, y2 - uy * margin + offY
	return (x1New, y1New), (x2New, y2New)


def getBoundingBox(allPoints, margin=VIEW_MARGIN):
	'''Get the bounding box for the figure given all the points drawn, with a given margin to accommodate them'''
	allX = [p[0] for p in allPoints]
	allY = [p[1] for p in allPoints]

	minX, maxX = min(allX), max(allX)
	minY, maxY = min(allY), max(allY)

	bboxX = [minX - margin, maxX + margin]
	bboxY = [minY - margin, maxY + margin]
	return bboxX, bboxY

def parsePDBAtomNumbers(ligFile):
		'''Parse the PDB to transform PDB atom numbers in RDKit atom numbers'''
		d = {}
		with open(ligFile) as f:
				for i, line in enumerate(f):
						d[int(line[6:11])] = i
		return d

def removeNumbers(s):
	'''Remove digits form string'''
	return ''.join(c for c in s if not c.isdigit())


def getCirclesRepulsion(circleI, allCircles, interactDistance, molCenter, kRep):
		circle, i = circleI
		repulsionForce = np.zeros(2)
		for j, otherCircle in enumerate(allCircles):
			if i != j:
				direction = circle - otherCircle
				distance = np.linalg.norm(direction)
				if distance < interactDistance:
					if distance == 0:
						direction = circle - molCenter
						direction[0] = -direction[0]
						distance = 100

					direction /= distance
					repF = direction * kRep / (distance)
					repulsionForce += repF
		return repulsionForce

def getAtomsRepulsion(circle, atomPositions, atomInteractRadius, kRep):
	repulsionForce = np.zeros(2)
	for atomPos in atomPositions:
		direction = circle - atomPos
		distance = np.linalg.norm(direction)
		if distance < atomInteractRadius:
			direction /= distance
			repForce = direction * kRep / (distance)
			repulsionForce += repForce
	return repulsionForce

def optimizeCirclePositions(circlePositions, atomPositions, molCenter,
														interactDistance, atomInteractRadius, iterations=100):
	"""Optimize the initial positions of the resdiue circles to avoid overlaps using repulsion forces between circles and 
  other circles and atoms. 
  Originally, the circle positions are next to the interactions atoms.
  """
	circles = np.array(circlePositions)
	kRep = 1000

	for _ in range(iterations):
		repForces = []
		for i, circle in enumerate(circles):
			# Repulsion against other circles
			repulsionForce = getCirclesRepulsion((circle, i), circles, interactDistance, molCenter, kRep)
			# Repulsion against closeby atoms
			repulsionForce += getAtomsRepulsion(circle, atomPositions, atomInteractRadius, kRep)

			# Update position
			circles[i] += repulsionForce
			repForces.append(sum(repulsionForce))

		# If the positions are not updated much anymore, exit
		if max(repForces) < 0.1:
			break
	return circles

def getDisplacements(n):
	a = list(range(n))
	meanA = sum(a) / len(a)
	return [i-meanA for i in a]


def normalizeAngle(angle, ofx, ofy):
	"""
  Normalize the angle in range[-90, 90] degrees so the text is always upside down
  """
	angle = angle % 360
	if angle > 180:
		angle -= 360

	if angle > 90:
		angle -= 180
		ofx, ofy = -ofx, -ofy
	elif angle < -90:
		angle += 180
		ofx, ofy = -ofx, -ofy

	return angle, ofx, ofy


def addParallelLabel(ax, x, y, texto, color, distancia=5, offSign=1, **kwargs):
	"""
  Add a label parallel to the interaction line, with the desired color at a desired perpedicular distance to the line.
  Offsign can be used to control whether the label will be placed in one or the other side of the line

  """
	dx = x[-1] - x[0]
	dy = y[-1] - y[0]
	angle = np.degrees(np.arctan2(dy, dx))

	length = np.sqrt(dx ** 2 + dy ** 2)
	if length > 0:
		ux = -dy / length  
		uy = dx / length  
	else:
		ux, uy = 0, 1

	midX = np.mean(x)
	midY = np.mean(y)
	offsetX = distancia * ux * offSign
	offsetY = distancia * uy * offSign
	angle, offsetX, offsetY = normalizeAngle(angle, offsetX, offsetY)

	label = ax.annotate(texto, xy=(midX, midY), xytext=(offsetX, offsetY), color=color, textcoords='offset points',
							rotation=angle, ha='center', va='center', annotation_clip=False, **kwargs)
	return label


class MoleculeInteractions:
	'''Main object to draw a molecule and its interacting residues, obtained from OpenMMDL.
	The residues are placed surrounding the moelcule, next to their interacting atoms, but their positions can be 
	dynamically changed.
	'''
	def __init__(self, molFile, inDic, poseFile=None, valueUnit=DEF_VALUE_UNIT):
		self.canvasSize = 2000
		self.circleRadius = 40
		self.interactDistance = 200
		self.valueUnit = valueUnit
		
		self.circles, self.resCircles = {}, {}
		self.cirLabels = {}
		self.selectedCircle = None
		
		self.mol = self.loadMolecule(molFile, poseFile)
		self.atomNumDic = parsePDBAtomNumbers(molFile)

		self.nDic = self.processInputDic(inDic)
		self.drawMolecule()
		self.cirDic, allPoints = self.getDrawnPoints(self.nDic)

		self.getBoundingBoxes(allPoints)
		self.addLegends()
		self.drawAminoacids()
		self.drawInteractions()

	########## Actions tools ########################
	def onPress(self, event):
		if event.inaxes != self.ax:
			return

		for circle in self.circles:
			contains, _ = circle.contains(event)
			if contains:
				self.selectedCircle = circle
				break

	def onRelease(self, event):
		self.selectedCircle = None

	def onMotion(self, event):
		if self.selectedCircle is None:
			return

		xData, yData = self.ax.transData.inverted().transform((event.x, event.y))

		#  Moving residue Circle and label
		self.selectedCircle.center = (xData, yData)
		self.cirLabels[self.selectedCircle].set_position((xData, yData))

		# Redoing interactions in new position
		self.removeCircleInteractions(self.selectedCircle)
		self.drawCircleInteractions(self.selectedCircle)

		self.fig.canvas.draw()

		if hasattr(self, 'plot_element'):
			self.plot_element.set_cursor_data(False)

	########## Molecule tools ########################
	def loadMolecule(self, molFile, poseFile):
		"""Molecule to draw.

		The PDB fixes the atom numbering the interactions refer to, but carries no bond orders, so the
		pose file (sdf, mol2...) is used as a template to recover them and draw the right molecule.
		"""
		molsDict, _ = getMolFilesDic([molFile], True)
		mol = next(iter(molsDict))

		if poseFile:
			try:
				tmplDict, _ = getMolFilesDic([poseFile], True)
				mol = AllChem.AssignBondOrdersFromTemplate(next(iter(tmplDict)), mol)
			except Exception as e:
				print(f'Bond orders could not be taken from {poseFile} ({type(e).__name__}: {e})')

		return mol

	########## Main tools ########################
	def drawMolecule(self):
			'''Draw the moelcule using RDKit'''
			rdDepictor.Compute2DCoords(self.mol)
			orientHorizontally(self.mol)
			Chem.RemoveStereochemistry(self.mol)
			self.mol = rdMolDraw2D.PrepareMolForDrawing(self.mol)

			self.drawer = rdMolDraw2D.MolDraw2DCairo(int(self.canvasSize * 0.5), int(self.canvasSize * 0.5))
			self.drawer.DrawMolecule(self.mol)
			self.drawer.FinishDrawing()

	def getDrawnPoints(self, nDic):
		'''Places the residues given the molecule atoms positions and returns that information
		'''
		molPoints = self.getMolPoints()
		cirDic = self.getCirclesDic(nDic)

		centerX, centerY = getMoleculeCenter(self.drawer, self.mol, self.canvasSize)
		newCircles = optimizeCirclePositions(list(cirDic.values()), molPoints, [centerX, centerY],
											 self.interactDistance, self.interactDistance)

		for i, resId in enumerate(cirDic):
			newCircle = newCircles[i]
			molPoints.append(newCircle)
			cirDic[resId] = newCircle

		return cirDic, molPoints

	def getBoundingBoxes(self, allPoints):
			'''Sets the bounding boxes based on the drawn positions'''
			bboxX, bboxY = getBoundingBox(allPoints)
			bboxWidth, bboxHeight = bboxX[1] - bboxX[0], bboxY[1] - bboxY[0]

			self.drawer.WriteDrawingText('mol_temp_full.png')
			imgFull = mpimg.imread('mol_temp_full.png')
			imgFull = np.flipud(imgFull)

			matplotlib.use('tkagg')
			self.fig, self.ax = plt.subplots(figsize=(10, 10 * bboxHeight / bboxWidth))
			# The plot takes the whole window and not the part of it matplotlib leaves free around the axes
			self.fig.subplots_adjust(left=0.01, right=0.99, bottom=0.01, top=0.99)
			self.ax.imshow(imgFull)
			self.ax.set_xlim(bboxX[0], bboxX[1])
			self.bboxY = bboxY
			self.ax.set_ylim(bboxY[0], bboxY[1])
			self.ax.set_axis_off()

			self.fig.canvas.mpl_connect('button_press_event', self.onPress)
			self.fig.canvas.mpl_connect('button_release_event', self.onRelease)
			self.fig.canvas.mpl_connect('motion_notify_event', self.onMotion)

	def drawAminoacids(self):
		'''Draw all the interacting resiudes and their labels'''
		for resId, (circleX, circleY) in self.cirDic.items():
			circleCoords = (circleX, circleY)
			# Wide enough for the name it holds, so "GLY101" is not written outside of its own box
			width = max(len(resId) * self.circleRadius * 0.55, self.circleRadius * 2)
			circle = ResidueBox(circleCoords, width, self.circleRadius * 1.2,
											facecolor=COLOR_DIC[removeNumbers(resId.upper())],
											linewidth=2, alpha=0.2, transform=self.ax.transData, clip_on=False)
			self.ax.add_patch(circle)
			cirLabel = self.ax.text(circleX, circleY, resId, color='black', weight='bold', fontsize=8,
															ha='center', va='center', transform=self.ax.transData)

			self.resCircles[resId] = circle
			self.circles[circle] = resId
			self.cirLabels[circle] = cirLabel

		self.fitBoxesToLabels()

	def fitBoxesToLabels(self):
		'''Each box is given in the coordinates of the drawing and its name in points, so the one
		only holds the other at a certain scale: the boxes are fitted to the names as they are
		drawn, and the view widened for those ending up over its border. Widening it changes the
		scale, so it is all done again'''
		for _ in range(2):
			self.fig.canvas.draw()
			for circle, label in self.cirLabels.items():
				labelBox = label.get_window_extent().transformed(self.ax.transData.inverted())
				circle.setSize(labelBox.width + labelBox.height * 0.9, labelBox.height * 1.9)

			pad = max([circle.size[1] for circle in self.circles])
			corners = [(x, y) for circle in self.circles
						 for x in (circle.get_x(), circle.get_x() + circle.size[0])
						 for y in (circle.get_y(), circle.get_y() + circle.size[1])]
			xs, ys = [c[0] for c in corners], [c[1] for c in corners]
			xLims, yLims = self.ax.get_xlim(), self.ax.get_ylim()
			self.ax.set_xlim(min(xLims[0], min(xs) - pad), max(xLims[1], max(xs) + pad))
			self.ax.set_ylim(min(yLims[0], min(ys) - pad), max(yLims[1], max(ys) + pad))

	def drawInteractions(self):
		'''Draw all the interactions between residues and moelcule atoms'''
		self.circleInteractions = {}
		for resId in self.cirDic:
			circle = self.resCircles[resId]
			self.drawCircleInteractions(circle)

	def addLegends(self):
		'''Add residues and interaction legends above the graph, where no residue can fall on them'''
		residuePatches = []
		for resGroup, color in RES_COLORS.items():
			element = Patch(facecolor=color, alpha=0.4, edgecolor='black', label=resGroup, linewidth=1)
			residuePatches.append(element)


		interactionPatches = []
		for intType, color in COLOR_INT_DIC.items():
			element = Line2D([0], [0], color=color, linewidth=2.5,
												linestyle='--', label=intType.capitalize())
			interactionPatches.append(element)

		# Both in a single block, each one a column of its own under its heading
		heading = lambda title: Patch(facecolor='none', edgecolor='none', label=title)
		nRows = max(len(residuePatches), len(interactionPatches)) + 1
		columns = ([heading('Residues type')] + residuePatches, [heading('Interactions type')] + interactionPatches)
		handles = [h for column in columns for h in column + [heading(' ')] * (nRows - len(column))]

		legend = self.ax.legend(handles=handles, ncol=len(columns), frameon=True, loc='upper left',
											handlelength=1.8, columnspacing=1.5)

		for row in (0, nRows):
			legend.get_texts()[row].set_fontweight('bold')

		for _ in range(3):
			self.fig.canvas.draw()
			legendBox = legend.get_window_extent().transformed(self.ax.transData.inverted())
			self.ax.set_ylim(self.bboxY[0], self.bboxY[1] + legendBox.height * 1.05)


	def display(self):
		'''Main function to display the generated graph'''
		plt.show()


	########## HELP tools ########################
	def processInputDic(self, intDic):
		'''Processes the input dictionary containing the interacting information'''
		nDic = {}
		for pairId, iDic in intDic.items():
			resId, atomIds = pairId[0].replace('_', ''), pairId[1]
			if '_' in atomIds:
				nAtomIds = '_'.join([str(self.atomNumDic[int(atomId)]) for atomId in atomIds.split('_')])
			else:
				nAtomIds = self.atomNumDic[int(eval(atomIds))]

			if resId not in nDic:
				nDic[resId] = {}
			nDic[resId][nAtomIds] = {}
			for intType, n in iDic.items():
				nDic[resId][nAtomIds][intType] = n
		return nDic

	def getCirclesDic(self, nDic):
		'''Gets the original positions of the circles, close to their interacting atoms'''
		cirDic = {}
		centerX, centerY = getMoleculeCenter(self.drawer, self.mol, self.canvasSize)
		for resId, atomDic in nDic.items():
			atomIdx = next(iter(atomDic))
			if isinstance(atomIdx, str):
				atomIdx = int(atomIdx.split('_')[0])
				iDist = self.interactDistance * 2
			else:
				iDist = self.interactDistance

			atomPoint = getDrawCoords(self.drawer, self.canvasSize, atomIdx)
			atomPoint = [atomPoint.x, atomPoint.y]
			circleX, circleY = getExternalPosition(atomPoint, (centerX, centerY), iDist)
			cirDic[resId] = [circleX, circleY]
		return cirDic
	
	def getMolPoints(self):
		'''Returns the molecule atom coordinates'''
		allPoints = []
		for i in range(self.mol.GetNumAtoms()):
			atomPoint = getDrawCoords(self.drawer, self.canvasSize, i)
			allPoints.append((atomPoint.x, atomPoint.y))
		return allPoints

	def addInteractLine(self, atomCoords, circleCoords, disp, intColor, n, offSign, margin=50):
		'''Adds a interaction line and its label given the circle and atom coordinates
		'''
		c1, c2 = setConnectionLineCoords(atomCoords, circleCoords, disp, margin=margin)
		xs, ys = [c1[0], c2[0]], [c1[1], c2[1]]
		line = Line2D(xs, ys, color=intColor, alpha=0.6, linewidth=1.5, linestyle='--', clip_on=False)
		self.ax.add_line(line)
		n = '{} {}'.format(round(float(n), 2), self.valueUnit).strip()
		label = addParallelLabel(self.ax, xs, ys, n, intColor, offSign=offSign)
		return line, label

	def drawCircleInteractions(self, circle):
		'''Functions that draws all the interactions of a given circle'''
		resId = self.circles[circle]
		circleCoords = circle.get_center()

		# Displaced over all the interactions of the residue and not over those of each atom, so lines
		# reaching different atoms do not lie on top of each other when they come out nearly parallel
		displacs = getDisplacements(sum([len(intsDic) for intsDic in self.nDic[resId].values()]))

		circleLines, di = [], 0
		for atomIdx, intsDic in self.nDic[resId].items():
			if isinstance(atomIdx, str):
				piCoords = []
				for atomAro in atomIdx.split('_'):
					atomPoint = getDrawCoords(self.drawer, self.canvasSize, int(atomAro))
					piCoords.append([atomPoint.x, atomPoint.y])

				piPoint = np.mean(np.array(piCoords), axis=0)
				atomCoords = (piPoint[0], piPoint[1])

			else:
				atomPoint = getDrawCoords(self.drawer, self.canvasSize, atomIdx)
				atomCoords = (atomPoint.x, atomPoint.y)

			for intType, n in intsDic.items():
				offSign = 1 if di % 2 == 0 else -1
				intColor = COLOR_INT_DIC[intType.lower()]
				border = circle.getBorderDistance(circleCoords[0] - atomCoords[0],
																		  circleCoords[1] - atomCoords[1])
				linea, label = self.addInteractLine(atomCoords, circleCoords, displacs[di], intColor, n,
																	 offSign, margin=border + 10)
				circleLines.append([linea, label])
				di += 1

		self.circleInteractions[circle] = circleLines

	def removeCircleInteractions(self, circle):
		'''Remove the existing interactions of a circle'''
		for line, label in self.circleInteractions[circle]:
			line.remove(), label.remove()


if __name__ == "__main__":
		pDic = parseParams(sys.argv[1], sep='::')
		molFile, intDic = pDic['molFile'], eval(pDic['intDic'])
		# Optional: molecule with the bond orders a PDB lacks
		poseFile = pDic.get('poseFile', None)
		valueUnit = pDic.get('valueUnit', DEF_VALUE_UNIT)

		app = MoleculeInteractions(molFile, intDic, poseFile=poseFile, valueUnit=valueUnit)
		app.display()

