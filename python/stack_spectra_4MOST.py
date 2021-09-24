#! /usr/bin/env python

"""
This script produces the stacks for samples defined by a list of spectra identifiers.

nohup python3 stack_spectra_ELG.py > stack_spectra_ELG.log &

"""
import sys
import os 
from os.path import join
import glob
import numpy as n
import SpectraStackingEBOSS as sse

stack_dir = join(os.environ['HOME'],"SDSS/stacks")

def stack_it( specList ):
	outfile = join(stack_dir, os.path.basename(specList)+".stack")
	print(stack_dir, outfile)
	stack=sse.SpectraStackingEBOSS(specList, outfile )
	stack.createStackMatrix()
	stack.stackSpectra()

list_2_stack = n.array(glob.glob(join(stack_dir, "*.ascii")))
for el in list_2_stack:
	stack_it(el)

