"""
For generating human-readable output reports
"""
import sys
import os
import math

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure

from fsl.data.image import Image

class RstContent(object):
    """
    Simple helper class for creating documents containing ReStructuredText
    """

    def __init__(self):
        self._content = ""
        self._heading_chars = "=-~+"

    def image(self, fname):
        """
        Add a block-level image
        """
        self._content += ".. image:: %s\n\n" % fname
        
    def text(self, txt):
        """
        Add a line of text content
        """
        self._content += txt + "\n\n"

    def heading(self, txt, level=0):
        """
        Add a heading
        """
        if level >= len(self._heading_chars):
            raise ValueError("Unsupported heading level: %i" % level)
        self._content += txt + "\n"
        self._content += self._heading_chars[level] * len(txt) + "\n\n"
        
    def __str__(self):
        return self._content

class Report(object):
    """
    A report consisting of .rst documents and associated images
    which can be turned into a document (HTML, PDF etc)
    """

    def __init__(self, report_dir):
        self._dir = report_dir
        if not os.path.exists(report_dir):
            os.makedirs(report_dir)
        elif not os.path.isdir(report_dir):
            raise ValueError("%s exists but is not a directory" % report_dir)

    def add_rst(self, fname, src):
        """
        Write a .rst text file

        :param fname: File name to write to
        :param src: Object whose string representation (``str(srouce)`` provides the text content
                    Typically this will be an RstContent object but it does not have to be.
        """
        if not fname.endswith(".rst"):
            fname += ".rst"
        with open(os.path.join(self._dir, fname), "w") as rstfile:
            rstfile.write(str(src))

    def add_lightbox_img(self, fname, *imgs, **kwargs):
        """
        Write a .png image file showing a lightbox view of one or more Image instances

        :param fname: File name to write to
        :param *imgs: One or more ``fsl.data.Image`` instances. Later images will be 
                      overlaid onto earlier images
        """
        if not imgs:
            raise ValueError("At least one image required")

        shape = None
        for img in imgs:
            if not isinstance(img, Image):
                raise ValueError("Images must be instances of fsl.data.Image")
            if shape is None:
                shape = img.shape
            #print(img.axisMapping(img.worldToVoxMat))
            if img.ndim != 3:
                raise ValueError("Images must be 3D") 
            if img.shape != shape:
                raise ValueError("Images do not have consistent shapes")
        
        num_slices = min(9, shape[2])
        grid_size = int(math.ceil(math.sqrt(num_slices)))

        fig = Figure()
        FigureCanvas(fig)
        for nslice in range(num_slices):
            axes = fig.add_subplot(grid_size, grid_size, nslice+1)
            axes.set_yticklabels([])
            axes.set_xticklabels([])
            axes.set_xticks([])
            axes.set_yticks([])
            for img in imgs:
                slice_idx = int(float(shape[2]*nslice)/num_slices)
                axes.imshow(img.data[:, :, slice_idx].T)
            
        fig.subplots_adjust(wspace=0, hspace=0.05)
        fig.savefig(os.path.join(self._dir, fname))

def main():
    """
    Simple command line for testing
    """
    report = Report(".")
    report.lightbox_img("test.png", *[Image(fname) for fname in sys.argv[1:]])

if __name__ == "__main__":
    main()
