# -*- coding: utf-8 -*-
#  **************************************************************************
# *
# * Authors:     Daniel Del Hoyo Gomez (ddelhoyo@cnb.csic.es)
# *              Martín Salinas Antón (martin.salinas@cnb.csic.es)
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
from .base import *
from .reactions import *


# ================================================================
# Custom ROI voting objects (Irene Sánchez, 2025)
# ================================================================

import pyworkflow.object as pwobj
from pwem.objects import EMObject, EMSet, Pointer


class ROIVote(EMObject):
    """Objeto individual: residuo + frecuencia + porcentaje."""
    _possibleAttributes = ['_residue', '_frequency', '_percentage']

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._residue = pwobj.String()
        self._frequency = pwobj.Integer()
        self._percentage = pwobj.Float()


class SetOfROIVotes(EMSet):
    """Set global para ROI Voting."""
    _targets = [ROIVote]
    _label = 'Set of ROI Votes'

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._itemType = ROIVote


# ================================================================
# Custom ROI color mapping objects (Irene Sánchez, 2025)
# ================================================================

class ROIColorItem(EMObject):
    """Single residue with color intensity information."""
    _possibleAttributes = ['_residue', '_frequency', '_percentage', '_color']

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._residue = pwobj.String()
        self._frequency = pwobj.Integer()
        self._percentage = pwobj.Float()
        self._color = pwobj.String()  # e.g., 'red', 'blue', etc.


class SetOfROIColorMap(EMSet):
    """Set for color-mapped residues."""
    _targets = [ROIColorItem]
    _label = 'Set of ROI Color Map'

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._itemType = ROIColorItem


# ================================================================
# Custom ROI frequency filter objects (Irene Sánchez, 2025)
# ================================================================

class ROIFilterItem(EMObject):
    """Single filtered residue with frequency and percentage."""
    _possibleAttributes = ['_residue', '_frequency', '_percentage']

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._residue = pwobj.String()
        self._frequency = pwobj.Integer()
        self._percentage = pwobj.Float()


class SetOfROIFiltered(EMSet):
    """Clean set with only 3 columns (residue, frequency, percentage)."""
    _targets = [ROIFilterItem]
    _label = 'Set of ROI Filtered Residues'

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._itemType = ROIFilterItem

