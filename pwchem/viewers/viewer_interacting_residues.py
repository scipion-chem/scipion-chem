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
from matplotlib import pyplot as plt, cm
from matplotlib.patches import Circle, Patch
from pyworkflow.protocol.params import *
import pyworkflow.viewer as pwviewer

from pwchem.protocols.VirtualDrugScreening.protocol_define_manual_structROIs import ProtDefineStructROIs
from pwchem.viewers import ViewerGeneralStructROIs


class InteractingResViewer(ViewerGeneralStructROIs):
    _label = 'Viewer interacting residues'
    _targets = [ProtDefineStructROIs]
    _environments = [pwviewer.DESKTOP_TKINTER]

    def _defineParams(self, form):
        form.addSection(label='Residue interaction view')
        form.addParam('distanceMinThreshold', FloatParam,
                      label='Min distance to display (Å)',
                      default=0.0,
                      help='Only show residue pairs with mean distance above this threshold.')
        form.addParam('distanceMaxThreshold', FloatParam,
                      label='Max distance to display (Å)',
                      default=5.0,
                      help='Only show residue pairs with mean distance below this threshold.')
        form.addParam('labelDistances', BooleanParam,
                      label='Label distances',
                      default=True,
                      help='Display mean distance labels on the connecting lines.')

    def _getVisualizeDict(self):
        return {'labelDistances': self._viewResidueInteractions}

    def _viewResidueInteractions(self, paramName=None):
        df = self.load_interaction_data(self.distanceMinThreshold.get(), self.distanceMaxThreshold.get())
        if df is None:
            return []

        all_chains, x_positions, residues_by_chain, y_positions_by_chain = self.build_positions(df)

        fig, ax = plt.subplots(figsize=(2 + 2 * len(all_chains), 10))
        ax.set_xticks(list(x_positions.values()))
        ax.set_xticklabels([f'Chain {ch}' for ch in all_chains])
        ax.set_ylabel('Residues')
        ax.set_title('Residue interactions across chains')

        node_artists, text_artists, chain_colors = self.create_nodes(ax, residues_by_chain, x_positions, y_positions_by_chain)
        connections = self.create_connections(ax, df, node_artists, self.labelDistances.get())

        legend_handles = [Patch(facecolor=color, edgecolor='black', label=f'Chain {chain}')
                          for chain, color in chain_colors.items()]

        ax.legend(handles=legend_handles, title='Chains', loc='upper right', frameon=True)

        ax.set_xlim(-1, len(all_chains))
        ax.set_ylim(-1, max(len(v) for v in residues_by_chain.values()) + 1)
        ax.set_aspect('auto')
        plt.tight_layout()

        self.enable_interactivity(fig, ax, node_artists, text_artists, connections)

        plt.show()
        return []


#----------------UTILS-----------------
    def load_interaction_data(self, min_distance, max_distance):
        """Load and filter interacting residues."""
        csv_path = Path(self.protocol._getExtraPath("interacting_residues.csv"))
        if not os.path.exists(csv_path):
            raise FileNotFoundError(f"Interactions file not found: {csv_path}")

        df = pd.read_csv(csv_path)
        df['Mean distance'] = df['Mean distance'].astype(float)
        df = df[df['Mean distance'] >= min_distance]
        df = df[df['Mean distance'] <= max_distance]

        if df.empty:
            print("No interactions within the specified distance threshold.")
            return None

        return df


    def build_positions(self, df):
        """Compute x/y positions for residues grouped by chain."""
        all_chains = sorted(set(df['Chain1']).union(df['Chain2']))
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

        return all_chains, x_positions, residues_by_chain, y_positions_by_chain


    def create_nodes(self, ax, residues_by_chain, x_positions, y_positions_by_chain):
        """Draw residue nodes and labels. Colored by chain."""
        colors = cm.get_cmap('Set2', len(residues_by_chain))
        chain_colors = {}

        node_artists, text_artists = {}, {}
        for i, (chain, residues) in enumerate(residues_by_chain.items()):
            color = colors(i)
            chain_colors[chain] = color

            for res in residues:
                x, y = x_positions[chain], y_positions_by_chain[chain][res]
                circle = Circle(
                    (x, y),
                    radius=0.2,
                    facecolor=color,
                    edgecolor='black',
                    zorder=3,
                    picker=True
                )
                ax.add_patch(circle)
                text = ax.text(x, y, res.split(':')[-1],
                               ha='center', va='center', fontsize=9, color='black', zorder=4)
                node_artists[(chain, res)] = circle
                text_artists[(chain, res)] = text

        return node_artists, text_artists, chain_colors


    def create_connections(self, ax, df, node_artists, label_distances):
        """Establish connections between residues. Labels with distances if selected."""
        connections = []
        for _, row in df.iterrows():
            ch1, ch2 = row['Chain1'], row['Chain2']
            res1, res2 = row['Residue1'], row['Residue2']
            mean_dist = row['Mean distance']
            if (ch1, res1) not in node_artists or (ch2, res2) not in node_artists:
                continue
            x1, y1 = node_artists[(ch1, res1)].center
            x2, y2 = node_artists[(ch2, res2)].center
            line, = ax.plot([x1, x2], [y1, y2], linestyle=':', color='gray', alpha=0.6, zorder=1)

            label_obj = None
            if label_distances:
                label_obj = ax.text((x1 + x2) / 2, (y1 + y2) / 2,
                                    f'{mean_dist:.1f}Å',
                                    fontsize=8, color='black',
                                    ha='center', va='center', zorder=2)

            connections.append((line, ch1, res1, ch2, res2, mean_dist, label_obj))
        return connections


    def enable_interactivity(self, fig, ax, node_artists, text_artists, connections):
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
            for (chain, res), circ in node_artists.items():
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
                    'artist': node_artists[(ch, rs)], 'chain': ch, 'res': rs, 'press': (event.xdata, event.ydata)
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
            text_artists[(dragged_node['chain'], dragged_node['res'])].set_position(circ.center)

            # Update connected lines and distance labels
            for line, ch1, res1, ch2, res2, mean_dist, label_obj in connections:
                if (ch1, res1) == (dragged_node['chain'], dragged_node['res']):
                    x2, y2 = node_artists[(ch2, res2)].center
                    line.set_data([circ.center[0], x2], [circ.center[1], y2])
                    if label_obj:
                        label_obj.set_position(((circ.center[0] + x2) / 2, (circ.center[1] + y2) / 2))
                elif (ch2, res2) == (dragged_node['chain'], dragged_node['res']):
                    x1, y1 = node_artists[(ch1, res1)].center
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
