# qmost_templates

Imports and dependencies
```
pip install numpy, astropy, matplotlib, speclite, spectres, scipy, pyneb
```

Clone the repository
```
cd ~/software/linux
git clone git@github.com:JohanComparat/qmost_templates.git
```

Go to the directory 
```
cd ~/software/linux/qmost_templates/python
```

Step 1

Create the list of spectra to be stacked together
```
python create_stack_list_4MOST.py
```

Stack the spectra (long step). It uses the library SpectraStackingEBOSS.py and SDSS spectra from the DR16 https://dr16.sdss.org/
You might need to add this python directory to your python path to be able to use the library.
```
python stack_spectra_4MOST.py
```

As a function of redshift, stitch together the spectra. 
For each stitched stack, it is needed to input by hand the starting and ending wavelength used in each individual stack.
Then a routine makes a figure of the stack obtained (it uses lineListVac.py).
```
python stack_spectra_4MOST_stitch_full_wavelength.py
python plot_qmost_stack.py
```

Finally convert to a 4MOST template following these specifications https://4most.mpe.mpg.de/static/manual/spectral_template_format_specification.htm

```
python convert_stack_2_qmost_template.py
```

Then copy the results to this repository (intermediate files are large)
```
cd ~/software/linux/qmost_templates
rsync -avz /data43s/SDSS/stacks/X_AGN/*stitched* data/
```

and commit to the github.

