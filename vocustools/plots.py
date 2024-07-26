import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from typing import List, Union

def color_oplot(x: List[float], y: List[float], z: List[float], z2: List[float], 
                markersize: float, xlegend: float = 0.85, ylegend: float = 0.75, 
                nlabels: int = 5, lsize: float = 1.5, lthick: float = 2):
    # Line 13058-13060: Set default values for optional parameters
    if ylegend < 0.26:
        ylegend = 0.26

    # Line 13061: Get length of input data
    length = len(x)

    # Lines 13062-13063: Normalize z and z2 values
    z_resize = [(val - min(z)) * 255.0 / (max(z) - min(z)) for val in z]
    z2_resize = [(val - min(z2)) * 2.0 / (max(z2) - min(z2)) + 1.0 for val in z2]

    # Lines 13064-13069: Plot colored markers
    fig, ax = plt.subplots()
    for i in range(length - 1):
        markersize2 = markersize
        color = plt.cm.viridis(z_resize[i] / 255.0)
        ax.plot(x[i:i+2], y[i:i+2], 'o-', color=color, markersize=markersize2, 
                linewidth=max(3, markersize2 * 4.0))

    # Lines 13072-13094: Create color legend
    sm = plt.cm.ScalarMappable(cmap=plt.cm.viridis, norm=plt.Normalize(vmin=min(z), vmax=max(z)))
    cbar = plt.colorbar(sm, ax=ax, location='right', pad=0.1)
    cbar.set_ticks(np.linspace(min(z), max(z), nlabels))
    cbar.set_ticklabels([f"{val:.2f}" for val in np.linspace(min(z), max(z), nlabels)])

    # Set plot limits and labels
    ax.set_xlim(min(x), max(x))
    ax.set_ylim(min(y), max(y))
    ax.set_xlabel('X-axis')
    ax.set_ylabel('Y-axis')
    ax.set_title('Color Plot')

    plt.tight_layout()
    plt.show()

def SaveImage(filename: str, dataset: List[List[float]], colortable: int) -> None:
    """
    Save a 2D array as a JPEG image.
    
    Lines 25638-25656: Function logic
    """
    import numpy as np
    from PIL import Image

    dataset = np.array(dataset)
    dataset = (dataset - dataset.min()) / (dataset.max() - dataset.min()) * 255
    dataset = dataset.astype(np.uint8)

    if colortable > -0.5:
        # Apply color table (this would need a separate implementation)
        # For now, we'll just use a grayscale image
        img = Image.fromarray(dataset, mode='L')
    else:
        img = Image.fromarray(dataset, mode='L')

    img.save(filename)
