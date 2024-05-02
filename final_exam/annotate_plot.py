"""
This module annotates plots using matplotlib's Pyplot text function based on specified annotations.
__author__ = "Jacob, Noah"
"""

import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime

def annotate_plot(annotations):
    """
    Annotates a plot using the provided annotations dictionary.

    Parameters:
    annotations (dict): A dictionary where keys are labels to be annotated and values are dictionaries containing:
        'position' (ndarray): The x, y coordinates for the position of the text box.
        'alignment' (list of strings): Contains horizontalalignment and verticalalignment values.
        'fontsize' (float): The font size of the text.
    """
    for label, properties in annotations.items():
        plt.text(
            properties['position'][0], properties['position'][1],
            label,
            horizontalalignment=properties['alignment'][0],
            verticalalignment=properties['alignment'][1],
            fontsize=properties['fontsize']
        )

if __name__ == "__main__":
    # Create a simple plot to demonstrate annotation
    plt.figure()
    x = np.linspace(0, 10, 100)
    y = np.sin(x)
    plt.plot(x, y, label='Sine wave')
    plt.xlabel('X axis')
    plt.ylabel('Y axis')
    plt.title('Test Plot with Annotations')

    # Annotations dictionary
    annotations = {
        "Created by Jacob, Noah ({})".format(datetime.now().date().isoformat()): {
            'position': np.array([0, -1.3]),  # Relative to axes coordinates
            'alignment': ['left', 'top'],
            'fontsize': 10
        }
    }

    # Call the annotate function
    annotate_plot(annotations)

    # Adjust plot to fit annotations and show plot
    plt.ylim(-1.5, 1.5)  # Adjusting y limits to show annotation clearly
    plt.legend()
    plt.grid(True)
    plt.show()
