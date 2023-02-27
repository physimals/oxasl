"""
OXASL - Quantification using Vaby (Python based VB)

Unlike the fabber method this does not yet support PVC

Copyright (c) 2008-2020 Univerisity of Oxford
"""
from . import multistep_fit

def run(wsp):
    # Basic non-PVC run
    multistep_fit.run(wsp.sub("basil"), impl="vaby")
    wsp.quantify_wsps.append("basil")
    