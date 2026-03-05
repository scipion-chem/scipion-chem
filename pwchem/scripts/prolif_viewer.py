import os
import sys
import pickle as pkl
import matplotlib.pyplot as plt
import webbrowser
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='View ProLIF analysis results.')
    # The first argument is the path to the pkl
    parser.add_argument('pkl_path', type=str, help='Path to the _fingerprint.pkl file')
    # The second argument defines what to show
    parser.add_argument('--mode', choices=['barcode', 'network'], default='barcode',
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
    if not os.path.exists(fp_path):
        print(f"Error: Pickle file not found at {fp_path}")
        sys.exit(1)

    with open(fp_path, 'rb') as f:
        fp = pkl.load(f)

    df = fp.to_dataframe()
    if df.empty:
        print("No interactions detected in the loaded file.")
        sys.exit(0)

    print("Generating Barcode Plot...")
    fp.plot_barcode()
    plt.title("ProLIF Interaction Barcode")
    plt.tight_layout()
    plt.show()