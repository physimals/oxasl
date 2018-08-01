"""
For generating human-readable output reports

Example usage:

    report = Report()

    page = ReportPage()
    page.heading("Registration to standard space")
    page.heading("Structural image", level=1)
    page.text("This image was used as the registration source")

    page.heading("Standard image", level=1)
    page.text("This image was used to obtain the registration target")

    page.heading("Transformation", level=1)
    page.heading("Affine transformation structural->standard", level=2)
    page.maths()
    page.heading("Affine transformation standard->structural", level=2)
    page.maths()
    
    report.add_rst("registration"< page)
    report.generate_html("my-report")

"""
import sys
import os
import math
import datetime
import shutil
import six

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

    def maths(self, content):
        """
        Write mathematical content
        """
        self.text(".. math::")
        if isinstance(content, six.string_types):
            self.text("    " + content)
        else:
            for line in content:
                self.text("    " + line)

    def heading(self, txt, level=0):
        """
        Add a heading
        """
        if level >= len(self._heading_chars):
            raise ValueError("Unsupported heading level: %i" % level)
        self._content += txt + "\n"
        self._content += self._heading_chars[level] * len(txt) + "\n\n"
        
    def table(self, name, tabdata):
        self._content += ".. csv-table:: " + name + "\n\n"
        import io
        import csv
        csvtxt = io.BytesIO()
        writer = csv.writer(csvtxt)
        for row in tabdata:
            writer.writerow(row)
        for line in csvtxt.getvalue().splitlines():
            self._content += "    " + line + "\n"
        self._content += "\n"

    def __str__(self):
        return self._content

REPORT_CONF = """
# This file is execfile()d with the current directory set to its
# containing dir.

extensions = ['sphinx.ext.mathjax',]
templates_path = ['_templates']
source_suffix = '.rst'
master_doc = 'index'
project = u'oxasl'
copyright = u'2018, oxasl'
author = u'oxasl'
version = u''
release = u''
language = None
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']
pygments_style = 'sphinx'
todo_include_todos = False

html_theme = 'alabaster'
# html_theme_options = {}
html_static_path = ['_static']
html_sidebars = {
    '**': [
        'relations.html',  # needs 'show_related': True theme option to display
        'searchbox.html',
    ]
}
htmlhelp_basename = 'oxasldoc'

# -- Options for LaTeX output ---------------------------------------------

latex_elements = {
    # 'papersize': 'letterpaper',
    # 'pointsize': '10pt',
    # 'preamble': '',
    # 'figure_align': 'htbp',
}

latex_documents = [
    (master_doc, 'oxasl.tex', u'oxasl Documentation',
     u'oxasl', 'manual'),
]
"""

INDEX_TEMPLATE = """
OXASL processing report
=======================

Start time: %s

End time: %s

.. toctree::
   :maxdepth: 1
   :caption: Contents:

%s
"""

class Report(object):
    """
    A report consisting of .rst documents and associated images
    which can be turned into a document (HTML, PDF etc)
    """

    def __init__(self, report_dir):
        self._dir = report_dir
        self._rst_files = []
        self._start_time = datetime.datetime.now()
        self._end_time = None

        if not os.path.exists(report_dir):
            os.makedirs(report_dir)
        elif not os.path.isdir(report_dir):
            raise ValueError("%s exists but is not a directory" % report_dir)

    def generate_html(self, dest_dir):
        """
        Generate an HTML report
        """
        self._end_time = datetime.datetime.now()
        with open(os.path.join(self._dir, "conf.py"), "w") as conffile:
            conffile.write(REPORT_CONF)

        with open(os.path.join(self._dir, "index.rst"), "w") as indexfile:
            rst_files = "\n".join(["  %s" % rst_file for rst_file in self._rst_files])
            indexfile.write(INDEX_TEMPLATE % (self._start_time.strftime("%Y-%m-%d %H:%M:%S"), self._end_time.strftime("%Y-%m-%d %H:%M:%S"), rst_files))

        os.system('sphinx-build -M html "%s" "%s"' % (self._dir, os.path.join(self._dir, "_build")))
        if not os.path.exists(dest_dir):
            os.makedirs(dest_dir)
        elif not os.path.isdir(dest_dir):
            raise ValueError("Report destination is not a directory: %s" % dest_dir)

        shutil.copytree(os.path.join(self._dir, "_build", "html"), os.path.join(dest_dir, "report"))

    def add_rst(self, fname, src):
        """
        Write a .rst text file

        :param fname: File name to write to
        :param src: Object whose string representation (``str(srouce)`` provides the text content
                    Typically this will be an RstContent object but it does not have to be.
        """
        if not fname.endswith(".rst"):
            fname += ".rst"
        self._rst_files.append(fname[:-4])
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
