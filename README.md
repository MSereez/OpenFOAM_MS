# OpenFOAM_MS
These are some source codes and makefiles that will be useful in the OpenFOAM for incompressible flow in Aerodynamics

1) Spalart Allamaras Low Re version was taken from a publication in AIAA (year 2020). See the header file, for more details. This version
improves the prediction of skin friction coefficient and therefore, the dynamic stall of 2D,3D configurations (atleast for my test cases).

2) Conical Motion - A solid body motion function for conical motion i.e. rotation around wind velocity vector

3) forces - This is used for extracting the forces which are aligned with the body in a rotating grid. The function uses the angular velocity in rad/s and 
a rotation vector i.e. the axis around which rotation occurs.. to create a quaternion and then rotate the forces with respect to this quaternion. A sample control dict file has been provided. Note that, the moments and forces will be correctly extracted only if the rotation rates and axis in the dynamicMeshDict and the control dict match each other. Furthermore, the CofR, liftDir,dragDir, etc.. that you usually provide with openfoam force definition needs still to be provided. 

4) I also plan to upload the dual time stepping solver for incompressible flow which uses the transient SIMPLE algorithm along with local time stepping, but as it is an ongoing project there will be some delays.

please contact mohamed.sereez@dmu.ac.uk if you find any bugs or have some comments.

