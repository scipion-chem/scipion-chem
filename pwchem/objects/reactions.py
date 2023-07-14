# -*- coding: utf-8 -*-
#  **************************************************************************
# *
# * Authors:     Martín Salinas Antón (martin.salinas@cnb.csic.es)
# *
# *
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

# General imports
from typing import Tuple

# Scipion em imports
import pwem.objects.data as data
import pyworkflow.object as pwobj

# Scipion chem imports
from pwchem.objects import SmallMolecule, SetOfSmallMolecules

class ChemicalReaction(data.EMObject):
	"""
	This class represents a chemical reaction, where one or more
	chemical compounds (SmallMolecules) produce other compounds.
	"""
	def __init__(self, inputMolecules: SetOfSmallMolecules=None, outputMolecules: SetOfSmallMolecules=None,
			name: pwobj.String = '', energyRequired: pwobj.Float = 0.0, energyUnits: pwobj.String = None, **kwargs):
		"""
		### Constructor for the ChemicalReaction class.

		#### Parameters:
		inputMolecules (SetOfSmallMolecules): Optional. Input molecules for the chemical reaction.
		outputMolecules (SetOfSmallMolecules): Optional. Output molecules produced by the chemical reaction.
		name (String): Optional. Name of the chemical reaction.
		energyRequired (Float): Optional. Required energy for the reaction to take place.
		energyUnits (String): Optional. Units for the required energy.
		"""
		# Calling parent constructor to set its variables
		super().__init__(self, **kwargs)

		# Setting reaction name
		self.__name = name

		# Setting the required energy for the reaction to take place
		self.__energyRequired = energyRequired
		self.__energyUnits = energyUnits

		# Setting input and output molecules
		self.__inputMolecules = inputMolecules
		self.__outputMolecules = outputMolecules

	#--------------------------------------- PRIVATE FUNCTIONS ---------------------------------------#
	def __getMoleculeFromSet(self, moleculeSet: SetOfSmallMolecules, name: pwobj.String = '', fileName: pwobj.String = ''):
		"""
		### This function returns the molecule in the given set that matches all of the provided attributes (name and/or filename).

		#### Parameters:
		moleculeSet (SetOfSmallMolecules): Set to look the molecules from.
		name (String): Optional. Name of the molecule to look for.
		fileName (String): Optional. File name of the molecule to look for.

		#### Returns:
		(SmallMolecule): The molecule matching the given parameters. None if it was not found or if neither a name or a filename are provided.
		"""
		# At least one of the search params needs to be introduced
		if not name and not fileName:
			return None
		
		# Search in the molecule set for the first one that matches the required attribute/s
		for molecule in moleculeSet:
			# Search by name if filename is missing
			# Search by filename if name is missing
			# Match both attributes if both are present
			if (
				(not fileName and name == molecule.getMolName()) or
				(not name and fileName == molecule.getFileName()) or
				(name == molecule.getMolName() and fileName == molecule.getFileName())
			):
				return molecule

	#--------------------------------------- PUBLIC FUNCTIONS ---------------------------------------#
	################## NAME ATTRIBUTE FUNCTIONS ##################
	def getName(self) -> pwobj.String:
		"""
		### This function returns the name of the chemical reaction.

		#### Returns:
		(String): The name of the chemical reaction.
		"""
		return self.__name

	def setName(self, name: pwobj.String):
		"""
		### This function sets the name of the chemical reaction.

		#### Parameters:
		name (String): The name of the chemical reaction.
		"""
		self.__name = name
	
	################## INPUT MOLECULES ATTRIBUTE FUNCTIONS ##################
	def getInputMolecules(self) -> SetOfSmallMolecules:
		"""
		### This function returns the input molecules of the chemical reaction.

		#### Returns:
		(SetOfSmallMolecules): The set of input molecules.
		"""
		return self.__inputMolecules
	
	def setInputMolecules(self, inputMolecules: SetOfSmallMolecules):
		"""
		### This function sets the input molecules of the chemical reaction.

		#### Parameters:
		inputMolecules (SetOfSmallMolecules): The set of input molecules.
		"""
		self.__inputMolecules = inputMolecules

	def addInputMolecule(self, inputMolecule: SmallMolecule):
		"""
		### This function adds the given SmallMolecule to the input set.

		#### Parameters:
		inputMolecule (SmallMolecule): Molecule to be added to the input set.
		"""
		# If input set has not been yet initialized, raise exception
		if not self.__inputMolecules:
			raise AttributeError("Set of input molecules has not been initialized yet. You need to call setInputMolecules first.")
		
		# Adding molecule to set
		self.__inputMolecules.append(inputMolecule)
	
	def getInputMolecule(self, name: pwobj.String = '', fileName: pwobj.String = ''):
		"""
		### This function returns the molecule in the input set that matches all of the provided attributes (name and/or filename).

		#### Parameters:
		name (String): Optional. Name of the input molecule to look for.
		fileName (String): Optional. File name of the input molecule to look for.

		#### Returns:
		(SmallMolecule): The input molecule matching the given parameters. None if it was not found or if neither a name or a filename are provided.
		"""
		return self.__getMoleculeFromSet(self.getInputMolecules(), name=name, fileName=fileName)

	################## OUTPUT MOLECULES ATTRIBUTE FUNCTIONS ##################
	def getOutputMolecules(self) -> SetOfSmallMolecules:
		"""
		### This function returns the output molecules of the chemical reaction.

		#### Returns:
		(SetOfSmallMolecules): The set of output molecules.
		"""
		return self.__outputMolecules
	
	def setOutputMolecules(self, outputMolecules: SetOfSmallMolecules):
		"""
		### This function sets the output molecules of the chemical reaction.

		#### Parameters:
		outputMolecules (SetOfSmallMolecules): The set of output molecules.
		"""
		self.__outputMolecules = outputMolecules
	
	def addOutputMolecule(self, outputMolecule: SmallMolecule):
		"""
		### This function adds the given SmallMolecule to the output set.

		#### Parameters:
		outputMolecule (SmallMolecule): Molecule to be added to the output set.
		"""
		# If output set has not been yet initialized, raise exception
		if not self.__outputMolecules:
			raise AttributeError("Set of output molecules has not been initialized yet. You need to call setOutputMolecules first.")
		
		# Adding molecule to set
		self.__outputMolecules.append(outputMolecule)
	
	def getOutputMolecule(self, name: pwobj.String = '', fileName: pwobj.String = ''):
		"""
		### This function returns the molecule in the output set that matches all of the provided attributes (name and/or filename).

		#### Parameters:
		name (String): Optional. Name of the output molecule to look for.
		fileName (String): Optional. File name of the output molecule to look for.

		#### Returns:
		(SmallMolecule): The output molecule matching the given parameters. None if it was not found or if neither a name or a filename are provided.
		"""
		return self.__getMoleculeFromSet(self.getOutputMolecules(), name=name, fileName=fileName)
	
	################## ENERGY ATTRIBUTE FUNCTIONS ##################
	def getRequiredEnergy(self) -> Tuple[pwobj.Float, pwobj.String]:
		"""
		### This function returns the required energy of the chemical reaction.

		#### Returns:
		tuple(Float, String): Energy required and its units.
		"""
		return self.__energyRequired, self.__energyUnits

	def setRequiredEnergy(self, energy: pwobj.Float, units: pwobj.String):
		"""
		### This function sets the required energy of the chemical reaction and its units.

		#### Parameters:
		energy (Float): The amount of energy required for the chemical reaction to take place.
		units (String): The units for that amount of energy.
		"""
		self.__energyRequired = energy
		self.__energyUnits = units
	
	def printRequiredEnergy(self) -> str:
		"""
		### This function returns the required energy of the chemical reaction in a formatted string.

		#### Returns:
		str: Energy required and its units in a formatted string.
		"""
		# Getting stored data
		energy, units = self.getRequiredEnergy()

		# If units have not been set, leave blank
		units = f' {units}' if units else ''

		# Return formatted string
		return f'{energy}{units}'