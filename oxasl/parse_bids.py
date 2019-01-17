"""
BIDS support

Initially we simply want to parse a BIDS data set and use it to set up our
ASL image and other workspace attributes.

In the longer run we could base the workspace around BIDS concepts entirely
"""
from __future__ import absolute_import

import sys

import bids

def parse(fname):
    dataset = bids.BIDSLayout(fname)
    print(dataset.as_data_frame())
    for subject in dataset.get_subjects():
        print("Subject: ", subject)
        for session in dataset.get_sessions():
            print("Session: ", session)
        else:
            fs = dataset.get(subject=subject, extensions=".nii.gz")
            for f in fs:
                print("File: ", f)
                print(dataset.get_metadata(f.filename))

if __name__ == "__main__":
    parse(sys.argv[1])
