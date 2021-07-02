# PolygonGrainTrack

Track regular polygons in granular materials using correlation template matching, discussed in Particle dynamics in two-dimensional point-loaded granular media composed of circular or pentagonal grains, Ryan  Kozlowski, Hu  Zheng, Karen E.  Daniels, Joshua E. S.  Socolar, EPJ Web Conf. 249 06010 (2021), DOI: 10.1051/epjconf/202124906010

Code was written originally in Matlab R2019a and has been tested in Matlab R2017b as well. The Computer Vision Toolbox is necessary for use of "insertShape" to visualize detected polygons over experimental image, but this could be bypassed by displaying the experimental image (imread) and then plotting the detected polygons or disks (plot, viscircles).

Sample images from experimental data of a slider pulled over a gravity-packed bed of bidisperse pentagons or disks is provided. The cameras are close to the imaging plane and a lens distortion correction has been applied; because of the proximity of the cameras, grain edges can be seen, which make detection more difficult. The detection is most successful when the particle edges are not visible, in the center of the image.

**To run a sample**: Download the entire directory with the code MaskProduction.m, GrainTrack.m, and the folder ProcessedImages/. Run GrainTrack from this directory, with either option 1 or 2 in Line 13 for setnumber = ;. Option 1 processes a pentagon image, Option 2 a disk image. Polygons take longer to track given that a stacked correlation has to be computed with polygons at varying angles. 

Modifications for your data:
(1) MaskProduction.m is written to produce polygon or disk masks, in this case for a bidisperse system. Lower angularity polygons than heptagons have to be written manually, and resolution of mask needs to be adjusted for resolution of experimental images. Angular resolution will also need to be adjusted.
(2) The GrainTrack.m code is written for bidisperse polygon packings, that is, to detect the "large" grains first and then "small" grains. For more size-disperse systems, one must be careful to write code that first finds the largest grains and iteratively detectes smaller ones. Otherwise the correlation will yield false peaks.
(3) GrainTrack.m can be modified in the case that the number of grains is known. The current example does not have a fixed number of grains since the cameras are moving along a channel of grains; a threshold in the correlation image is instead used to determine when to stop searching for more grains. If the number of grains is known, one can straightforwardly implement a conditional to stop the search once the quota has been detected.
(4) Masks other than polygons can be made, such as with ellipses/rods. Utilize the symmetry of the shape to avoid adding computation time to correlations with masks at different angles.
(5) Post-processing of the experimental image prior to or after binarization might help with edge-edge contacts, such as edge detection with a Sobel algorithm. 

