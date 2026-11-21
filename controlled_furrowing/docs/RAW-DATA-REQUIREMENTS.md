# Data needed for a specimen-resolved geometry test

No author has been contacted. This is a concrete data specification, not a claim that all conditional modeling is impossible without raw microscopy. Public numerical tables and an illustrative movie already support useful, explicitly limited tests.

For the stronger specimen-resolved test, retain the following for each small and larger activation record:

1. A raw-file/embryo identifier crosswalk into source Fig. 2 and Fig. 4; repeated observations must stay grouped.
2. Pre-activation 3D apical/basal surfaces, cell height and cell-size/topology information where resolvable, with length calibration and AP/DV orientation. Initial curvature is an input, not a parameter chosen from the final bulge.
3. Time-matched apical and basal profiles through the response peak and recovery, including nonactivated reference tissue used for cellularization correction. Separate maxima are insufficient for instantaneous shape reconstruction.
4. Full illumination geometry, scan times, depth range, and first-scan recruitment intensities/normalization metadata. Measured recruitment is distinct from commanded laser power and absolute active stress.
5. For bending labels: acquisition slice height/reference convention, optical axial resolution, time window, registration operations, and the image/record used to assign the plug-and-annulus category. The actual slice intersects the whole body, not just its apical displacement.
6. Mechanical perturbation/control records that constrain basal stiffness, prestress, outer-tissue response and any vitelline contact assumptions. Ablation velocity alone is not an absolute elastic modulus without a mechanical/drag model.

A conditional idealized finite-height model can be built now, carrying unknown geometry and material assumptions. These data are needed to assess whether its successful or failed predictions describe the actual specimens, rather than chosen idealized geometries.
