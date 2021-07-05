#!/bin/sh

DATADIR=/Users/ctsu0221/data/asl/fsl_course/
#DEBUG=--debug

oxasl -i ${DATADIR}/mpld_asltc --casl --iaf=tc --ibf=tis --slicedt=0.0452 \
      --tis=1.65,1.9,2.15,2.4,2.65,2.9 --bolus=1.4 \
      --fslanat=${DATADIR}/T1.anat --senscorr \
      -c ${DATADIR}/aslcalib --tr=4.8 --cmethod=single --csf=${DATADIR}/csfmask \
      --cblip=${DATADIR}/aslcalib_PA --echospacing=0.00952 --pedir=y \
      --mc --pvcorr \
      -o mpld_pvc_out --overwrite ${DEBUG}


