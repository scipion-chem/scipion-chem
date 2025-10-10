# **************************************************************************
# *
# * Authors:     Blanca Pueche (blanca.pueche@cnb.csic.es)
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
from matplotlib import pyplot as plt
from matplotlib.patches import Circle
from pyworkflow.protocol.params import *
import pyworkflow.viewer as pwviewer

from pwchem.protocols.VirtualDrugScreening.protocol_define_manual_structROIs import ProtDefineStructROIs



class InteractingResViewer(pwviewer.ProtocolViewer):
  _label = 'Viewer interacting residues'
  _targets = [ProtDefineStructROIs]
  _environments = [pwviewer.DESKTOP_TKINTER]

  def __init__(self, **kwargs):
    pwviewer.ProtocolViewer.__init__(self, **kwargs)

  def _defineParams(self, form):
      form.addSection(label='Residue interaction view')
      form.addParam('distanceThreshold', FloatParam,
                    label='Max distance to display (Å)',
                    default=5.0,
                    help='Only show residue pairs with mean distance below this threshold.')
      form.addParam('labelDistances', BooleanParam,
                    label='Label distances',
                    default=True,
                    help='Display mean distance labels on the connecting lines.')

  def _getVisualizeDict(self):
      return {
          'labelDistances': self._viewResidueInteractions
      }


  def _viewResidueInteractions(self, paramName=None):
      csv_path = Path(self.protocol._getExtraPath("interacting_residues.csv"))
      if not os.path.exists(csv_path):
          raise FileNotFoundError(f"Interactions file not found: {csv_path}")

      df = pd.read_csv(csv_path)
      df['Mean distance'] = df['Mean distance'].astype(float)
      df = df[df['Mean distance'] <= self.distanceThreshold.get()]

      if df.empty:
          print("No interactions within the specified distance threshold.")
          return []

      all_chains = sorted(set(df['Chain1']).union(df['Chain2']))
      n_chains = len(all_chains)

      x_positions = {chain: i for i, chain in enumerate(all_chains)}

      residues_by_chain = {}
      for chain in all_chains:
          res_chain1 = df[df['Chain1'] == chain]['Residue1'].unique()
          res_chain2 = df[df['Chain2'] == chain]['Residue2'].unique()
          residues = sorted(set(res_chain1).union(res_chain2),
                            key=lambda r: int(''.join([c for c in r if c.isdigit()]) or 0))
          residues_by_chain[chain] = residues

      y_positions_by_chain = {
          chain: {res: i for i, res in enumerate(residues)}
          for chain, residues in residues_by_chain.items()
      }

      fig, ax = plt.subplots(figsize=(2 + 2 * n_chains, 10))
      ax.set_xticks(list(x_positions.values()))
      ax.set_xticklabels([f'Chain {ch}' for ch in all_chains])
      ax.set_ylabel('Residues')
      ax.set_title('Interactive residue interactions across chains')

      node_artists = {}  # (chain, resname) -> Circle
      text_artists = {}
      connections = []  # (line, chain1, res1, chain2, res2)

      for chain, residues in residues_by_chain.items():
          for res in residues:
              x, y = x_positions[chain], y_positions_by_chain[chain][res]
              circle = Circle((x, y), radius=0.3, facecolor='skyblue', edgecolor='black', zorder=3, picker=True) #todo color by chain
              ax.add_patch(circle)
              text = ax.text(x, y, res.split(':')[-1], ha='center', va='center', fontsize=8, zorder=4)
              node_artists[(chain, res)] = circle
              text_artists[(chain, res)] = text

      for _, row in df.iterrows():
          ch1, ch2 = row['Chain1'], row['Chain2']
          res1, res2 = row['Residue1'], row['Residue2']
          mean_dist = row['Mean distance']

          if (ch1, res1) not in node_artists or (ch2, res2) not in node_artists:
              continue
          x1, y1 = node_artists[(ch1, res1)].center
          x2, y2 = node_artists[(ch2, res2)].center

          line, = ax.plot([x1, x2], [y1, y2], linestyle=':', color='gray', alpha=0.6, zorder=1)
          if self.labelDistances.get():
              mid_x, mid_y = (x1 + x2) / 2, (y1 + y2) / 2
              ax.text(mid_x, mid_y, f'{mean_dist:.1f}Å', fontsize=9, color='black', ha='center', va='center')
          connections.append((line, ch1, res1, ch2, res2))

      ax.set_xlim(-1, n_chains)
      ax.set_ylim(-1, max(len(v) for v in residues_by_chain.values()) + 1)
      ax.set_aspect('auto')
      plt.tight_layout()

      # ------------------ Interactivity ------------------

      dragged_node = {'artist': None, 'chain': None, 'res': None, 'press': None}

      def on_press(event):
          if event.inaxes != ax:
              return
          for (chain, res), circ in node_artists.items():
              contains, _ = circ.contains(event)
              if contains:
                  dragged_node['artist'] = circ
                  dragged_node['chain'] = chain
                  dragged_node['res'] = res
                  dragged_node['press'] = (event.xdata, event.ydata)
                  break

      def on_motion(event):
          if dragged_node['artist'] is None or event.inaxes != ax:
              return
          circ = dragged_node['artist']
          dx = event.xdata - dragged_node['press'][0]
          dy = event.ydata - dragged_node['press'][1]
          x0, y0 = circ.center
          circ.center = (x0 + dx, y0 + dy)
          dragged_node['press'] = (event.xdata, event.ydata)

          text_artists[(dragged_node['chain'], dragged_node['res'])].set_position(circ.center)

          for line, ch1, res1, ch2, res2 in connections:
              if (ch1, res1) == (dragged_node['chain'], dragged_node['res']):
                  x2, y2 = node_artists[(ch2, res2)].center
                  line.set_data([circ.center[0], x2], [circ.center[1], y2])
              elif (ch2, res2) == (dragged_node['chain'], dragged_node['res']):
                  x1, y1 = node_artists[(ch1, res1)].center
                  line.set_data([x1, circ.center[0]], [y1, circ.center[1]])
          fig.canvas.draw_idle()

      def on_release(event):
          dragged_node['artist'] = None
          dragged_node['press'] = None

      fig.canvas.mpl_connect('button_press_event', on_press)
      fig.canvas.mpl_connect('motion_notify_event', on_motion)
      fig.canvas.mpl_connect('button_release_event', on_release)

      plt.show()
      return []

