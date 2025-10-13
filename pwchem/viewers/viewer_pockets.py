# **************************************************************************
# *
# * Authors:     Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 3 of the License, or
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
from pathlib import Path

import pandas as pd
import pyworkflow.protocol.params as params
import pyworkflow.viewer as pwviewer
from matplotlib import pyplot as plt, cm
from matplotlib.patches import Circle, Patch

from pwem.viewers.mdviewer.viewer import MDViewer
from pyworkflow.protocol import Protocol

from pwchem.objects import SetOfStructROIs
from pwchem.protocols.VirtualDrugScreening.protocol_define_manual_structROIs import ProtDefineStructROIs
from pwchem.viewers.viewers_data import BioinformaticsDataViewer, PyMolViewer, VmdViewPopen
from pwchem.constants import *
from pwchem.protocols import ProtocolConsensusStructROIs


class StructROIPointsViewer(pwviewer.Viewer):
  _label = 'Viewer structural ROIs'
  _environments = [pwviewer.DESKTOP_TKINTER]

  # _targets = [SetOfPockets]

  def _visualize(self, obj, bBox=False, **kwargs):
    hetatmFile = obj.buildPDBhetatmFile()
    pmlFile = obj.createPML(bBox=bBox)

    pymolV = PyMolViewer(project=self.getProject())
    return pymolV._visualize(pmlFile, cwd=os.path.dirname(pmlFile))


class ContactSurfaceViewer(pwviewer.Viewer):
  _label = 'Viewer contact surface'
  _environments = [pwviewer.DESKTOP_TKINTER]

  # _targets = [SetOfPockets]

  def _visualize(self, obj, bBox=False, **kwargs):
    pmlFile = obj.createSurfacePml(bBox=bBox)

    pymolV = PyMolViewer(project=self.getProject())
    return pymolV._visualize(pmlFile, cwd=os.path.dirname(pmlFile))


VOLUME_PYMOL, VOLUME_PYMOL_SURF = 0, 1


class ViewerGeneralStructROIs(pwviewer.ProtocolViewer):
  _label = 'Viewer structural ROIs'
  _targets = [SetOfStructROIs, ProtDefineStructROIs]

  def __init__(self, **kwargs):
    pwviewer.ProtocolViewer.__init__(self, **kwargs)

  def _defineParams(self, form):
    form.addSection(label='Visualization of structural ROIs')
    group = form.addGroup('Pymol General Viewer')
    group.addParam('displayAtomStruct', params.EnumParam,
                  choices=['PyMol (ROI Points)', 'PyMol (Contact Surface)'],
                  default=VOLUME_PYMOL,
                  display=params.EnumParam.DISPLAY_HLIST,
                  label='Display output AtomStruct with',
                  help='*PyMol*: display AtomStruct and structural ROIs as points / surface.'
                  )
    group.addParam('displayBBoxes', params.BooleanParam,
                  default=False, label='Display ROIs bounding boxes',
                  help='Display the bounding boxes in pymol to check the size for the localized docking')
    group.addParam('pocketRadiusN', params.FloatParam, label='Grid radius vs pocket radius: ',
                  default=1.1, condition='displayBBoxes',
                  help='The radius * n of each ROI will be used as grid radius')

    form.addSection(label='Table view')
    form.addParam('displayTable', params.LabelParam,
                  label='Display ROIs set and attributes in table format: ',
                  help='Display the ROIs set in the set in table format with their respective attributes')

    form.addSection(label='Residue interaction view')
    form.addParam('distanceMinThreshold', params.FloatParam,
                  label='Min distance to display (Å)',
                  default=0.0,
                  help='Only show residue pairs with mean distance above this threshold.')
    form.addParam('distanceMaxThreshold', params.FloatParam,
                  label='Max distance to display (Å)',
                  default=5.0,
                  help='Only show residue pairs with mean distance below this threshold.')
    form.addParam('labelDistances', params.BooleanParam,
                  label='Label distances',
                  default=True,
                  help='Display mean distance labels on the connecting lines.')

  def _getVisualizeDict(self):
    return {
      'displayAtomStruct': self._showAtomStruct,
      'displayTable': self._viewSet,
      'labelDistances': self._viewResidueInteractions
    }

  def _viewSet(self, e=None):
    if type(self.protocol) == SetOfStructROIs:
      molSet = self.protocol
    elif hasattr(self.protocol, 'outputStructROIs'):
      molSet = getattr(self.protocol, 'outputStructROIs')
    else:
      print('Cannot find outputStructROIs')

    try:
      setV = MDViewer(project=self.getProject())
    except:
      setV = BioinformaticsDataViewer(project=self.getProject())
    return setV._visualize(molSet)

  def _validate(self):
    return []

  # =========================================================================
  # Display interacting residues
  # =========================================================================

  def _viewResidueInteractions(self, paramName=None):
      df = self.load_interaction_data(self.distanceMinThreshold.get(), self.distanceMaxThreshold.get())
      if df is None:
          return []

      allChains, xPositions, reisudesByChain, yPositionsByChain = self.build_positions(df)

      fig, ax = plt.subplots(figsize=(2 + 2 * len(allChains), 10))
      ax.set_xticks(list(xPositions.values()))
      ax.set_xticklabels([f'Chain {ch}' for ch in allChains])
      ax.set_ylabel('Residues')
      ax.set_title('Residue interactions across chains')

      nodeArtists, textArtists, chainColors = self.create_nodes(ax, reisudesByChain, xPositions,
                                                                   yPositionsByChain)
      connections = self.create_connections(ax, df, nodeArtists, self.labelDistances.get())

      legendHandles = [Patch(facecolor=color, edgecolor='black', label=f'Chain {chain}')
                        for chain, color in chainColors.items()]

      ax.legend(handles=legendHandles, title='Chains', loc='upper right', frameon=True)

      ax.set_xlim(-1, len(allChains))
      ax.set_ylim(-1, max(len(v) for v in reisudesByChain.values()) + 1)
      ax.set_aspect('auto')
      plt.tight_layout()

      self.enable_interactivity(fig, ax, nodeArtists, textArtists, connections)

      plt.show()
      return []

      # ----------------UTILS-----------------

  def load_interaction_data(self, minDistance, maxDistance):
    """Load and filter interacting residues."""
    csvPath = Path(self._targets[1]._getExtraPath("interacting_residues.csv")) #todo how do i access the csv
    if not os.path.exists(csvPath):
        raise FileNotFoundError(f"Interactions file not found: {csvPath}")

    df = pd.read_csv(csvPath)
    df['Mean distance'] = df['Mean distance'].astype(float)
    df = df[df['Mean distance'] >= minDistance]
    df = df[df['Mean distance'] <= maxDistance]

    if df.empty:
        print("No interactions within the specified distance threshold.")
        return None

    return df

  def build_positions(self, df):
    """Compute x/y positions for residues grouped by chain."""
    allChains = sorted(set(df['Chain1']).union(df['Chain2']))
    xPositions = {chain: i for i, chain in enumerate(allChains)}

    reisudesByChain = {}
    for chain in allChains:
        res_chain1 = df[df['Chain1'] == chain]['Residue1'].unique()
        res_chain2 = df[df['Chain2'] == chain]['Residue2'].unique()
        residues = sorted(set(res_chain1).union(res_chain2),
                          key=lambda r: int(''.join([c for c in r if c.isdigit()]) or 0))
        reisudesByChain[chain] = residues

    yPositionsByChain = {
        chain: {res: i for i, res in enumerate(residues)}
        for chain, residues in reisudesByChain.items()
    }

    return allChains, xPositions, reisudesByChain, yPositionsByChain

  def create_nodes(self, ax, reisudesByChain, xPositions, yPositionsByChain):
    """Draw residue nodes and labels. Colored by chain."""
    colors = cm.get_cmap('Set2', len(reisudesByChain))
    chainColors = {}

    nodeArtists, textArtists = {}, {}
    for i, (chain, residues) in enumerate(reisudesByChain.items()):
        color = colors(i)
        chainColors[chain] = color

        for res in residues:
            x, y = xPositions[chain], yPositionsByChain[chain][res]
            circle = Circle(
                (x, y),
                radius=0.25,
                facecolor=color,
                edgecolor='black',
                zorder=3,
                picker=True
            )
            ax.add_patch(circle)
            res_id = res.split(':')[-1]
            label = f"{chain}:{res_id}"
            text = ax.text(
                x, y, label,
                ha='center', va='center',
                fontsize=8, fontweight='bold',
                color='black', zorder=4
            )
            nodeArtists[(chain, res)] = circle
            textArtists[(chain, res)] = text

    return nodeArtists, textArtists, chainColors

  def create_connections(self, ax, df, nodeArtists, label_distances):
    """Establish connections between residues. Labels with distances if selected."""
    connections = []
    for _, row in df.iterrows():
        ch1, ch2 = row['Chain1'], row['Chain2']
        res1, res2 = row['Residue1'], row['Residue2']
        mean_dist = row['Mean distance']
        if (ch1, res1) not in nodeArtists or (ch2, res2) not in nodeArtists:
            continue
        x1, y1 = nodeArtists[(ch1, res1)].center
        x2, y2 = nodeArtists[(ch2, res2)].center
        line, = ax.plot([x1, x2], [y1, y2], linestyle=':', color='gray', alpha=0.6, zorder=1)

        label_obj = None
        if label_distances:
            label_obj = ax.text((x1 + x2) / 2, (y1 + y2) / 2,
                                f'{mean_dist:.1f}Å',
                                fontsize=8, color='black',
                                ha='center', va='center', zorder=2)

        connections.append((line, ch1, res1, ch2, res2, mean_dist, label_obj))
    return connections

  def enable_interactivity(self, fig, ax, nodeArtists, textArtists, connections):
    dragged_node = {'artist': None, 'chain': None, 'res': None, 'press': None}
    highlighted_node = {'chain': None, 'res': None}

    def reset_highlights():
        """Reset all connection lines and labels to default style."""
        for line, ch1, res1, ch2, res2, mean_dist, label_obj in connections:
            line.set_color('gray')
            line.set_linewidth(1)
            line.set_alpha(0.6)
            if label_obj:
                label_obj.set_fontweight('normal')

    def highlight_connections(chain, res):
        """Highlight lines (and labels) connected to a specific residue."""
        reset_highlights()
        for line, ch1, res1, ch2, res2, mean_dist, label_obj in connections:
            if (ch1, res1) == (chain, res) or (ch2, res2) == (chain, res):
                line.set_color('black')
                line.set_linewidth(2.5)
                line.set_alpha(1.0)
                if label_obj:
                    label_obj.set_fontweight('bold')
        fig.canvas.draw_idle()

    def on_press(event):
        if event.inaxes != ax:
            highlighted_node.update({'chain': None, 'res': None})
            reset_highlights()
            fig.canvas.draw_idle()
            return

        clicked_node = None
        for (chain, res), circ in nodeArtists.items():
            contains, _ = circ.contains(event)
            if contains:
                clicked_node = (chain, res)
                break

        if clicked_node:
            ch, rs = clicked_node
            if highlighted_node['chain'] == ch and highlighted_node['res'] == rs:
                # Clicking same node deselects it
                highlighted_node.update({'chain': None, 'res': None})
                reset_highlights()
            else:
                highlighted_node.update({'chain': ch, 'res': rs})
                highlight_connections(ch, rs)

            # Prepare dragging
            dragged_node.update({
                'artist': nodeArtists[(ch, rs)], 'chain': ch, 'res': rs, 'press': (event.xdata, event.ydata)
            })

        else:
            highlighted_node.update({'chain': None, 'res': None})
            reset_highlights()
            fig.canvas.draw_idle()

    def on_motion(event):
        if dragged_node['artist'] is None or event.inaxes != ax:
            return

        circ = dragged_node['artist']
        dx = event.xdata - dragged_node['press'][0]
        dy = event.ydata - dragged_node['press'][1]
        x0, y0 = circ.center
        circ.center = (x0 + dx, y0 + dy)
        dragged_node['press'] = (event.xdata, event.ydata)

        # Move text label with node
        textArtists[(dragged_node['chain'], dragged_node['res'])].set_position(circ.center)

        # Update connected lines and distance labels
        for line, ch1, res1, ch2, res2, mean_dist, label_obj in connections:
            if (ch1, res1) == (dragged_node['chain'], dragged_node['res']):
                x2, y2 = nodeArtists[(ch2, res2)].center
                line.set_data([circ.center[0], x2], [circ.center[1], y2])
                if label_obj:
                    label_obj.set_position(((circ.center[0] + x2) / 2, (circ.center[1] + y2) / 2))
            elif (ch2, res2) == (dragged_node['chain'], dragged_node['res']):
                x1, y1 = nodeArtists[(ch1, res1)].center
                line.set_data([x1, circ.center[0]], [y1, circ.center[1]])
                if label_obj:
                    label_obj.set_position(((x1 + circ.center[0]) / 2, (y1 + circ.center[1]) / 2))

        fig.canvas.draw_idle()

    def on_release(event):
        dragged_node.update({'artist': None, 'press': None})

    # Connect events
    fig.canvas.mpl_connect('button_press_event', on_press)
    fig.canvas.mpl_connect('motion_notify_event', on_motion)
    fig.canvas.mpl_connect('button_release_event', on_release)

  # =========================================================================
  # ShowAtomStructs
  # =========================================================================

  def getOutputAtomStructFile(self):
    return os.path.abspath(self.protocol.outputAtomStruct.getFileName())

  def _showAtomStruct(self, paramName=None):
    if self.displayAtomStruct == VOLUME_PYMOL:
      return self._showAtomStructPyMolPoints()

    elif self.displayAtomStruct == VOLUME_PYMOL_SURF:
      return self._showAtomStructPyMolSurf()

  def _showAtomStructPyMolPoints(self):
    bBox = self.displayBBoxes.get()
    if bBox:
      bBox = self.pocketRadiusN.get()

    pymolV = StructROIPointsViewer(project=self.getProject())
    if type(self.protocol) == SetOfStructROIs:
      return pymolV._visualize(self.protocol, bBox=bBox)
    elif hasattr(self.protocol, 'outputStructROIs'):
      return pymolV._visualize(getattr(self.protocol, 'outputStructROIs'), bBox=bBox)

  def _showAtomStructPyMolSurf(self):
    bBox = self.displayBBoxes.get()
    if bBox:
      bBox = self.pocketRadiusN.get()

    pymolV = ContactSurfaceViewer(project=self.getProject())
    if type(self.protocol) == SetOfStructROIs:
      return pymolV._visualize(self.protocol, bBox=bBox)
    elif hasattr(self.protocol, 'outputStructROIs'):
      return pymolV._visualize(getattr(self.protocol, 'outputStructROIs'), bBox=bBox)

  def _showAtomStructPyMol(self, pmlFile, outDir):
    pymolV = PyMolViewer(project=self.getProject())
    return pymolV._visualize(pmlFile, cwd=outDir)


MIXED, FPOCKET, P2RANK, AUTOLIGAND, SITEMAP = 'Mixed', 'FPocket', 'P2Rank', 'AutoLigand', 'Sitemap'
VOLUME_VMD = 2

class ViewerConsensusStructROIs(pwviewer.ProtocolViewer):
    _label = 'Viewer consensus structural ROIs'
    _targets = [ProtocolConsensusStructROIs]

    def __init__(self, **kwargs):
      pwviewer.ProtocolViewer.__init__(self, **kwargs)

    def _defineParams(self, form):
      form.addSection(label='Visualization of consensus Structural ROIs')
      pGroup = form.addGroup('Pymol General Viewer')
      pGroup.addParam('outputSet', params.EnumParam,
                    choices=self.getChoices(), default=0,
                    label='Set of pockets output: ',
                    help='Set of pockets output to visualize'
                    )
      pGroup.addParam('setClass', params.StringParam,
                    label='Set of Pockets Class: ', default='-',
                    help='Pocket class in chosen set')

      pGroup.addParam('displayFPocket', params.EnumParam,
                    choices=['PyMol (Pocket Points)', 'PyMol (Contact Surface)', 'VMD'],
                    default=0, condition='setClass=="{}"'.format(FPOCKET),
                    display=params.EnumParam.DISPLAY_HLIST,
                    label='Display output Set of Pockets with: ',
                    help='*PyMol*: display Set of Pockets and pockets as points / surface.\n '
                         '*VMD*: display Set of Pockets with VMD.'
                    )

      pGroup.addParam('displayPyMol', params.EnumParam,
                    choices=['PyMol (Pocket Points)', 'PyMol (Contact Surface)'],
                    default=0, condition='setClass!="{}"'.format(FPOCKET),
                    display=params.EnumParam.DISPLAY_HLIST,
                    label='Display output Set of Pockets with: ',
                    help='*PyMol*: display Set of Pockets and pockets as points / surface'
                    )
      pGroup.addParam('displayBBoxes', params.BooleanParam,
                    default=False, label='Display pocket bounding boxes',
                    condition='not (setClass=="{}" and displayFPocket==0)'.format(FPOCKET),
                    help='Display the bounding boxes in pymol to check the size for the localized docking')
      pGroup.addParam('pocketRadiusN', params.FloatParam, label='Grid radius vs pocket radius: ',
                    default=1.1, condition='displayBBoxes and not (setClass=="{}" and displayFPocket==0)'.format(FPOCKET),
                    help='The radius * n of each pocket will be used as grid radius')
      
      tGroup = form.addGroup('Table view')
      tGroup.addParam('displayTable', params.LabelParam,
                      label='Display table view of set of Pockets: ',
                      help='Table view')

    def getChoices(self):
        outputLabels = []
        for oAttr in self.protocol.iterOutputAttributes():
          outputLabels.append(oAttr[0])
        outputLabels.sort()
        return outputLabels

    def _getVisualizeDict(self):
      return {
        'displayTable': self._showTable,
        'displayFPocket': self._showFPocketPockets,
        'displayPyMol': self._showStandardPockets,
      }

    def _validate(self):
      return []

    # =========================================================================
    # ShowAtomStructs
    # =========================================================================

    def getOutputAtomStructFile(self):
      return os.path.abspath(self.protocol.outputAtomStruct.getFileName())

    def _showFPocketPockets(self, paramName=None):
      if self.displayFPocket == VOLUME_PYMOL:
        return self._showAtomStructPyMolPoints()

      elif self.displayFPocket == VOLUME_PYMOL_SURF:
        return self._showAtomStructPyMolSurf()

      elif self.displayFPocket == VOLUME_VMD:
        return self._showAtomStructVMD()

    def _showStandardPockets(self, paramName=None):
      if self.displayPyMol == VOLUME_PYMOL:
        return self._showAtomStructPyMolPoints()

      elif self.displayPyMol == VOLUME_PYMOL_SURF:
        return self._showAtomStructPyMolSurf()

    #Display functions
    def _showTable(self, paramName=None):
      outPockets = getattr(self.protocol, self.getEnumText('outputSet'))
      try:
        setV = MDViewer(project=self.getProject())
      except:
        setV = BioinformaticsDataViewer(project=self.getProject())
      return setV._visualize(outPockets)

    def _showAtomStructPyMolPoints(self):
      bBox = self.displayBBoxes.get()
      if bBox:
        bBox = self.pocketRadiusN.get()

      outPockets = getattr(self.protocol, self.getEnumText('outputSet'))
      pymolV = StructROIPointsViewer(project=self.getProject())
      return pymolV._visualize(outPockets, bBox=bBox)

    def _showAtomStructPyMolSurf(self):
      bBox = self.displayBBoxes.get()
      if bBox:
        bBox = self.pocketRadiusN.get()

      outPockets = getattr(self.protocol, self.getEnumText('outputSet'))
      pymolV = ContactSurfaceViewer(project=self.getProject())
      return pymolV._visualize(outPockets, bBox=bBox)

    def _showAtomStructVMD(self):
      outPockets = getattr(self.protocol, self.getEnumText('outputSet'))
      outFile = outPockets.getProteinHetatmFile().split('/')[-1]
      tclFile = outFile.replace('_out.pdb', '.tcl')
      outDir = os.path.abspath(outPockets.getSetDir())
      args = '{} -e {}'.format(outFile, tclFile)

      return [VmdViewPopen(args, cwd=outDir)]
