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
    
    report.add("registration", page)
    report.generate_html("my-report")
"""
import sys
import os
import math
import datetime
import shutil
import tempfile
import warnings

import six
import numpy as np

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure

from fsl.data.image import Image

class LightboxImage(object):
    """
    A .png image file showing a lightbox view of one or more Image instances
    """

    def __init__(self, img, bgimage=None, mask=None, **kwargs):
        """
        :param *imgs: One or more ``fsl.data.Image`` instances. Later images will be 
                      overlaid onto earlier images
        :param zeromask: If True, treat zero values as transparent
        """
        self._img = img
        self._bgimage = bgimage
        self._mask = mask
        self._zeromask = kwargs.get("zeromask", True)
        self.extension = ".png"

    def _slicerange(self, img, shape):
        if img is not None:
            nonzero_slices = [idx for idx in range(shape[2]) if np.count_nonzero(img.data[:, :, idx]) > 0]
            if nonzero_slices:
               return min(nonzero_slices), max(nonzero_slices)
        else:
            return 0, shape[2]-1

    def tofile(self, fname):
        """
        Write image to a file
        """
        shape = None
        for img in [self._img, self._bgimage, self._mask]:
            if img is None:
                continue
            if not isinstance(img, Image):
                raise ValueError("Images must be instances of fsl.data.Image: %s" % img)
            if shape is None:
                shape = img.shape
            if img.ndim != 3:
                raise ValueError("Images must be 3D") 
            if img.shape != shape:
                raise ValueError("Images do not have consistent shapes")
        
        min_slice, max_slice = self._slicerange(self._img, shape)
        num_slices = min(16, max_slice - min_slice + 1)
        grid_size = int(math.ceil(math.sqrt(num_slices)))

        fig = Figure(figsize=(5, 5), dpi=200)
        FigureCanvas(fig)
        for nslice in range(num_slices):
            axes = fig.add_subplot(grid_size, grid_size, nslice+1)
            axes.set_yticklabels([])
            axes.set_xticklabels([])
            axes.set_xticks([])
            axes.set_yticks([])
            slice_idx = int(float((max_slice - min_slice + 1)*nslice)/num_slices) + min_slice

            if self._bgimage:
                data = self._bgimage.data[:, :, slice_idx].T
                axes.imshow(data, cmap='gray')
            
            if self._img:
                data = self._img.data[:, :, slice_idx].T
                if issubclass(data.dtype.type, np.integer):
                    cmap = "Set1"
                else:
                    cmap= "viridis"

                if self._mask:
                    data = np.ma.masked_array(data, self._mask.data[:, :, slice_idx].T == 0)
                elif self._zeromask:
                    data = np.ma.masked_array(data, data == 0)
                
                axes.imshow(data, cmap=cmap)
            
        fig.subplots_adjust(wspace=0, hspace=0.05)
        fig.savefig(fname, bbox_inches='tight')
        
class ReportPage(object):
    """
    Simple helper class for creating documents containing ReStructuredText
    """

    def __init__(self):
        self._content = ""
        self._heading_chars = "=-~+"
        self.extension = ".rst"
        
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
            content = content.splitlines()

        for line in content:
            self._content += "    " + line + "\n"
        self._content += "\n"

    def _latex_float(self, f, sf=3):
        """
        Format float in format suitable for Latex - nicked off StackOverflow!
        """
        pyformat = "{0:.%ig}" % sf
        float_str = pyformat.format(f)
        if "e" in float_str:
            base, exponent = float_str.split("e")
            return r"{0} \times 10^{{{1}}}".format(base, int(exponent))
        else:
            return float_str

    def matrix(self, mat, sf=3):
        """
        Add a matrix of numbers
        """
        matrix_latex = "\\begin{bmatrix}\n"
        for row in mat:
            matrix_latex += " & ".join([self._latex_float(v, sf) for v in row]) + " \\\\\n"
        matrix_latex += "\\end{bmatrix}\n"
        self.maths(matrix_latex)

    def heading(self, txt, level=0):
        """
        Add a heading
        """
        if level >= len(self._heading_chars):
            raise ValueError("Unsupported heading level: %i" % level)
        self._content += txt + "\n"
        self._content += self._heading_chars[level] * len(txt) + "\n\n"
        
    def table(self, tabdata, name="", headers=None):
        self._content += ".. csv-table:: " + name + "\n"
        if headers:
            self._content += "    :header: " + ",".join(['"%s"' % h for h in headers]) + "\n"
        self._content += "\n"

        import io
        import csv
        csvtxt = io.BytesIO()
        writer = csv.writer(csvtxt)
        for row in tabdata:
            writer.writerow(row)
        for line in csvtxt.getvalue().splitlines():
            self._content += "    " + line + "\n"
        self._content += "\n"

    def dicttable(self, dictionary):
        tabdata = dictionary.items()
        self.table(tabdata, headers=("Key", "Value"))
        
    def tofile(self, fname):
        """
        Write RST content to a file
        """
        with open(fname, "w") as rstfile:
            rstfile.write(self._content)

    def __str__(self):
        return self._content

class Report(object):
    """
    A report consisting of .rst documents and associated images
    which can be turned into a document (HTML, PDF etc)
    """

    def __init__(self, title="OXASL processing report"):
        self._contents = []
        self._files = {}
        self._start_time = datetime.datetime.now()
        self._end_time = None
        self.title = title
        self.extension = ""

    def generate_html(self, dest_dir, build_dir=None):
        """
        Generate an HTML report
        """
        self._end_time = datetime.datetime.now()
        duration = (self._end_time - self._start_time).total_seconds()

        if build_dir:
            if os.path.exists(build_dir):
                if not os.path.isdir(build_dir):
                    raise ValueError("Report build directory %s exists but is not a directory" % build_dir)
                else:
                    warnings.warn("Report build directory %s already exists" % build_dir)
            is_temp = False
        else:
            build_dir = tempfile.mkdtemp("_report")
            is_temp = True

        try:
            self.tofile(build_dir)

            with open(os.path.join(build_dir, "conf.py"), "w") as conffile:
                conffile.write(REPORT_CONF)

            os.system('sphinx-build -M html "%s" "%s"' % (build_dir, os.path.join(build_dir, "_build")))

            if os.path.exists(dest_dir):
                if not os.path.isdir(dest_dir):
                    raise ValueError("Report destination directory %s exists but is not a directory" % dest_dir)
                else:
                    warnings.warn("Report destination directory %s already exists - removing" % dest_dir)
                    shutil.rmtree(dest_dir)
                
            shutil.copytree(os.path.join(build_dir, "_build", "html"), dest_dir)
        finally:
            if is_temp:
                shutil.rmtree(build_dir)

    def tofile(self, build_dir):
        if not os.path.exists(build_dir):
            os.makedirs(build_dir)

        with open(os.path.join(build_dir, "index.rst"), "w") as indexfile:
            indexfile.write(self.title + "\n")
            indexfile.write("=" * len(self.title) + "\n\n")
            self._timings(indexfile)
            self._toc(indexfile)
            
        for fname, content in self._files.items():
            content.tofile(os.path.join(build_dir, fname))

    def add(self, name, content, overwrite=False):
        """
        Add content to a report

        :param name: Name of the content.
        :param content: Content object which has ``extension`` attribute and supports ``tofile()`` method
        :param overwrite: If True, and content already exists with the same ``name`` and extension, 
                          replace content. Otherwise an exception is thrown.
        """
        fname = name + content.extension
        if not overwrite and fname in self._files:
            raise ValueError("Content with file name '%s' already exists and overwrite=False" % fname)

        self._files[fname] = content
        if isinstance(content, ReportPage):
            self._contents.append(name)
        if isinstance(content, Report):
            self._contents.append(name + "/index")

    def _timings(self, indexfile):
        if self._start_time:
            indexfile.write("Start time: %s\n\n" % self._start_time.strftime("%Y-%m-%d %H:%M:%S"))
        if self._end_time: 
            indexfile.write("End time: %s\n\n" % self._end_time.strftime("%Y-%m-%d %H:%M:%S"))

    def _toc(self, indexfile):
        indexfile.write(".. toctree::\n")
        indexfile.write("  :maxdepth: 1\n")
        indexfile.write("  :caption: Contents:\n\n")
        for rst_file in self._contents:
            indexfile.write("  %s\n" % rst_file)

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

html_theme = 'sphinxdoc'
# html_theme_options = {}
html_static_path = []
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

def main():
    """
    Simple command line for testing
    """
    report = Report()
    report.add("lbox", LightboxImage(*[Image(fname) for fname in sys.argv[1:]]))

    page = ReportPage()
    page.image("lbox.png")
    report.add("test", page)

    report.generate_html("testreport")

if __name__ == "__main__":
    main()
