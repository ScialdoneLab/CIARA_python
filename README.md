# IDENTOM <br /> (Identification of rare cell types using ENTropy Of Mixing)

Python implementation of the entropy of mixing algorithm that integrates into scanpy's analysis with the AnnData format.

The package can be installed via pip:

`python -m pip install entropyofmixing`

## Tutorial

The package contains the two main functions `get_full_background()` and `entropy_mixing()` which should be imported via:

`from entropyofmixing import get_full_background, entropy_mixing`

The functions are designed to work with scanpy's AnnData objects. For an example tutorial check out the Human Gastrula IPython Notebook which is also part of this repository.

## R package

Based on original R package by Gabriele Lubatti: 

https://github.com/ScialdoneLab/EntropyMixingRnew
