import os
import sys
import pickle as pkl
import matplotlib.pyplot as plt
import webbrowser
import argparse
import pandas as pd
import seaborn as sns
from rdkit import DataStructs

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='View ProLIF analysis results.')
    parser.add_argument('pkl_path', type=str, help='Path to the _fingerprint.pkl file')
    parser.add_argument('--mode', choices=['barcode', 'network', 'matrix'], default='barcode',
                        help='Choose to show the barcode plot or the interactive network')

    args = parser.parse_args()
    fp_path = args.pkl_path

    # --- Mode: Network ---
    if args.mode == "network":
        # Derive the HTML path from the pkl path
        network_html = fp_path.replace('_fingerprint.pkl', '_network.html')

        if os.path.exists(network_html):
            print(f"Opening interactive network: {network_html}")
            # Use absolute path to ensure the browser finds it
            webbrowser.open(f'file://{os.path.abspath(network_html)}')
            sys.exit(0)
        else:
            print(f"Error: HTML file not found at {network_html}")
            sys.exit(1)

    # --- Mode: Barcode (or anything else requiring the pkl) ---
    elif args.mode == "barcode":
        if not os.path.exists(fp_path):
            print(f"Error: Pickle file not found at {fp_path}")
            sys.exit(1)

        barcode_png = fp_path.replace('_fingerprint.pkl', '_barcode.png')

        with open(fp_path, 'rb') as f:
            fp = pkl.load(f)

        df = fp.to_dataframe()
        if df.empty:
            print("No interactions detected in the loaded file.")
            sys.exit(0)

        fp.plot_barcode()
        plt.title("Ligand-Target Interaction Barcode")
        plt.tight_layout()
        plt.show()

    # --- Mode: Matrix ---
    elif args.mode == "matrix":
        with open(fp_path, 'rb') as f:
            fp = pkl.load(f)
        df = fp.to_dataframe()

        bitvectors = fp.to_bitvectors()
        similarity_data = [DataStructs.BulkTanimotoSimilarity(bv, bitvectors) for bv in bitvectors]
        similarity_matrix = pd.DataFrame(similarity_data, index=df.index, columns=df.index)

        fig, ax = plt.subplots(figsize=(7, 6))
        colormap = sns.diverging_palette(300, 145, s=90, l=80, sep=30, center="dark", as_cmap=True)

        sns.heatmap(
            similarity_matrix,
            ax=ax,
            square=True,
            cmap=colormap,
            vmin=0, vmax=1, center=0.5,
            xticklabels=max(1, len(df) // 10),
            yticklabels=max(1, len(df) // 10)
        )
        ax.invert_yaxis()
        plt.yticks(rotation="horizontal")
        plt.title("Tanimoto Similarity Heatmap")
        plt.tight_layout()
        plt.show()