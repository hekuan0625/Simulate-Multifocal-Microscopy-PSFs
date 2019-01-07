# Simulate-Multifocal-Microscopy-PSFs
Simulate 3D PSFs of a multifocal micrsocpe using a diffractive optical element, whose numerical data is stored in the subfolder.

`calcPsfMFM.m` is used to simulate MFM's 3D PSFs, and `calcPsfMicrsocope.m` is used to simulate the PSFs of a conventional microscope.

### Run Example

`$ psf = calcPsfMFM(psfSz, voxSz, lambdas, 'verbose', verbose, 'M', M, 'NA', NA, 'f_tl', f_tl)` <br/>
The meaning of those input parameters can be found inside the file.
