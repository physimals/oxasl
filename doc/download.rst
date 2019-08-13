Getting the OXASL software
==========================

To use OXASL you will need FSL - version 6.0 or later is strongly recommended.
See `FSL installation <https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FslInstallation>`_ 
for installation instructions.

I have FSL 6.0 or later
-----------------------

To install into the fslpython environment use::

    fslpython -m pip install oxasl oxasl_ve oxasl_mp --user

This installs the main package and the vessel-encoding plugin.
To check it is working, try running the main executable::

    oxasl --version

I have an older version of FSL
------------------------------

You will need to download a pre-built bundle package containing the OXASL code and
also the required updated FSL dependencies. This can be found on our GitHub
release page:

https://github.com/ibme-qubic/oxasl/releases

