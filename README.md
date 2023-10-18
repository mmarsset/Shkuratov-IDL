# Shkuratov-IDL
IDL script used to model S-type asteroid spectra and measure their ol/(ol+opx) silicate ratio

For scientific reproducibility of our analysis, we provide the IDL code SpexTot045-190Fin_NIR.pro used to model the mineralogy of asteroid spectra presented in the manuscript entitled "The Massalia asteroid family as the origin of ordinary L chondrites" submitted to Nature by Marsset et al. This is an adaptation of the Shkuratov model (Shkuratov et al. 1999, Icarus, vol. 137, no. 2, pp. 235–246) similar to the one used, e.g., by Vernazza et al. (2008) and Binzel et al. (2019).

We also provide:
- the optical constants of materials used to model the spectra (located in folder "ocs"),
- the averaged spectra of ordinary H, L and LL chondrites comptuted from spectral data retrieved from the RELAB database; https://sites.brown.edu/relab/ (located in folder "OC_spectra").

Asteroid spectra presented in our paper are available at http://smass.mit.edu.

This code was successfully run with IDL version 8.8.2 on macOS Monterey Version 12.1 and should be compatible with most OS.

To run the code, the following paths in the script must be updated:
- Lines 4, 194: variable ocdir should point to the folder of optical constants
- Line 205: variable dir should point to the folder of asteroid family spectra 

In addition, dir should contain a file named "list_nir.dat" listing the asteroid spectra that you wish the code to process. 

The output is a file named "out_nir_binzel.dat" listing the best-fit parameters of each spectral modelling. The second column of that file provides the ol/(ol+opx) silicate ratio discussed in the paper. In addition, a png figure is automatically generated for spot check. 
