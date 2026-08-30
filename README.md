# FJX_addX_Rewrite_inR
I've reworked Michael Prather's original fortran addX code (from UCI_fastJ_addX_73cx) in R, with more useful variable names and comments. I tried to do this a few years ago and never quite got it right but I came back to it today with chatGPT's help and figured out how to fix my code. 

The original code was written in an extremely condensed form for what I presume must have been a character-limited version of fortran. The primary use case I envision for this tool is for GEOS-Chem users who want to add new cross-sections and quantum yields to GEOS-Chem. 

The original code can be found at https://github.com/geoschem/Cloud-J/tree/main/tools/AddXs

The FAST-JX photolysis mechanism (Wild et al., 2000) efficiently represents tropospheric photolysis, and is commonly implemented in CTMs like GEOS-Chem (e.g. https://wiki.seas.harvard.edu/geos-chem/index.php/FAST-JX_v7.0_photolysis_mechanism).

The code required to format new quantum yield and cross-section data into FJX is complex and somewhat opaque. I've tried to re-render it into more comprehensible R code that mimics the behavior of the original fortran while being comprehensible. This code requires an input file formatted in the same way as Prather's original.

Note that Prather's code contains a feature that is note immediately obvious. For an example of this, please look at the file XMeVK_JPL11X.dat. In that file, both line 6 and line 13 contain three temperature values in Kelvin. IF THESE TEMPERATURE VALUES ARE NOT IDENTICAL, the code will interpolate between the temperature values, while keeping values within the Temperature range set in line 13. I do not understand why you would want this behavior in the code, but I have included it in this to be consistent. However, I strongly recommend against using it. Users may wish to verify that past input files they've used include consistent temperatures in both lines.

Wild, O., X. Zhu, and M. J. Prather, Fast-J: Accurate simulation of in- and below-cloud photolysis in tropospheric chemical models, J. Atmos. Chem., 37, 245–282, 2000.
