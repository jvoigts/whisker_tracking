# whisker_tracking

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4315318.svg)](https://doi.org/10.5281/zenodo.4315318)

Unsupervised whisker tracking in freely behaving mice/rats

This repository includes a copy of ConvNet by Sergey Demyanov ( https://github.com/sdemyanov/ConvNet ) - if you'd rather use a recent version please get the current code from github.

This is _not_ the code described in our paper 'Unsupervised whisker tracking in unrestrained behaving animals', J Voigts, B Sakmann, T Celikel, Journal of neurophysiology 100 (1), 504-515 ( http://jn.physiology.org/content/100/1/504.short ), but a more recent, and simpler method. Specifically, this method only runs a CNN to label whiskers and misses all the code in the j neurophys method that turns such a labeled images into spline representations of individual whiskers.

See also the Clack et al. 2012 method ( https://openwiki.janelia.org/wiki/display/MyersLab/Whisker+Tracking ) if you're interested in precise tracking of individual whiskers, and Per M Knutsen's whisker tracker ( https://github.com/pmknutsen/whiskertracker ),
for tracking of vibrissae, head position and angle, and limbs.

There are two main analysis steps in this implementation: Identifying animal presence and nose tracking, and whisker tracking.


# Nose tracking
This code identifies periods where a mouse is present and tracks rough animal outlines.

# Whisker tracking
This code runs a CNN to label whiskers and then runs a very simple hough transform to extract whisker positions and angles.
