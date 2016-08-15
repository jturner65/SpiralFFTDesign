# SpiralFFTDesign
(Java/Processing) Spiral/Drawn Stroke Embossed on Spiral Designer, using FFT window as guide

This repo is an eclipse project, with all libraries it needs to run (on a windows-based system).  It should also work on other systems
although it hasn't been tested.

This program allows the user to draw trajectories or design spirals and have them embossed on a parent spiral, which can be controlled
via 3 control points.  Included is an FFT-based visualization tool that reveals the regularity of the resultant image and is intended to
allow for free-hand design uniformity - when the resultant design is slightly "off" in its proporitions, the FFT image will reveal it 
with many bright lines radiating from the center.  By adjusting the design until the lines in the FFT image themselves are aligned the
artist can be sure of regularity in the resultant image.

The motivation behind this process was the hypothesis that the aesthetic value of a drawn image is how well it represents the intention
behind its creation, and one measure of that intention might be some measurable regularity and consistency in the orientation of the 
components within the image.  A visualized FFT is a great tool to quickly measure these quantities.  
