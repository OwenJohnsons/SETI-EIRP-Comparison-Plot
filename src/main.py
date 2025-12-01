'''
Code Purpose: Main script to plot SETI EIRP limits from survey definitions in a YAML file with survey configs. 
Author: Owen A. Johnson 
Date: December 2025 
'''
import matplotlib.pyplot as plt
import argparse
from functions import plot_surveys

def get_args():
    parser = argparse.ArgumentParser(
        description="Plot SETI EIRP limits from survey definitions in a YAML file."
    )
    parser.add_argument(
        "-i", "--input", type=str, default="surveys.yaml",
        help="Input YAML file with survey definitions (default: surveys.yaml)"
    )
    parser.add_argument(
        "-o", "--output", type=str, default="SETI-EIRP-limits-Comp.pdf",
        help="Output filename for the plot (default: SETI-EIRP-limits-Comp.pdf)"
    )
    parser.add_argument("-pub", "--publish", action="store_true",
        help="Format the plot for publication using science plot."
    )
    
    parser.add_argument("--twocolumn", action="store_true",
        help="Format the plot for two-column layout.")
    
    parser.add_argument("--no-show", action="store_true",
        help="Do not display the plot interactively.")
    
    return parser.parse_args()

def main():
    args = get_args()
    plot_surveys("surveys.yaml", output=args.output, publish=args.publish, twocolumn=args.twocolumn)

if __name__ == "__main__":
    main()
