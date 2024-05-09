"""
This module annotates plots using matplotlib's Pyplot text function based on specified annotations.
__author__ = "Jacob, Noah"
"""

import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime

def annotate_plot(ax, annotations):
    """
    Annotates a plot using the provided annotations dictionary, positioning based on axes coordinates.

    Parameters:
    ax (AxesSubplot): The matplotlib axes on which to apply the annotations.
    annotations (dict): A dictionary where keys are labels to be annotated and values are dictionaries containing:
        'position' (ndarray): The x, y coordinates for the position of the text box in axes coordinates.
        'alignment' (list of strings): Contains horizontalalignment and verticalalignment values.
        'fontsize' (float): The font size of the text.
    """
    for label, properties in annotations.items():
        ax.text(
            properties['position'][0], properties['position'][1],
            label,
            horizontalalignment=properties['alignment'][0],
            verticalalignment=properties['alignment'][1],
            fontsize=properties['fontsize'],
            transform=ax.transAxes  # This makes the position relative to the axes
        )

# Example usage within a plotting function
if __name__ == "__main__":
    fig, ax = plt.subplots()
    x = np.linspace(0, 10, 100)
    y = np.sin(x)
    ax.plot(x, y, label='Sine wave')
    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    ax.set_title('Test Plot with Annotations')

    annotations = {
        "Created by Jacob, Noah ({})".format(datetime.now().date().isoformat()): {
            'position': np.array([0, 0]),  # Relative to axes coordinates
            'alignment': ['left', 'top'],
            'fontsize': 10
        }
    }
    annotate_plot(ax, annotations)
    ax.set_ylim(-1.5, 1.5)  # Adjusting y limits if necessary
    ax.legend()
    ax.grid(True)
    plt.show()
