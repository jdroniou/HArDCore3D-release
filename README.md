# HArDCore3D
HArD::Core3D (Hybrid Arbitrary Degree::Core 3D) - Library to implement schemes with face and cell polynomial unknowns on 3D generic meshes.

This is the 3D version of the HArD::Core library (https://github.com/jdroniou/HArDCore). The implementation principles are the same as those of the 2D version (https://github.com/jdroniou/HArDCore2D-release); see in particular the README.md file in that repository. The documentation can be found at https://jdroniou.github.io/HArDCore3D-release/

The meshing (src/Mesh/StemMesh) and quadrature (src/Quadratures) modules in HArDCore3D are partially based on Marco Manzini's code available at https://github.com/gmanzini-LANL/PDE-Mesh-Manager.

The implementations in this library follow general principles described in the appendix of the book "The Hybrid High-Order Method for Polytopal Meshes: Design, Analysis, and Applications" (D. A. Di Pietro and J. Droniou. Number 19 in Modeling, Simulation and Applications, Springer International Publishing, 2020. ISBN 978-3-030-37202-6 (Hardcover) 978-3-030-37203-3 (eBook). url: https://hal.archives-ouvertes.fr/hal-02151813). High-order methods with hybrid unknowns have certain specificities which sometimes require fine choices, e.g. of basis functions (hierarchical, orthonormalised or not), etc. We refer to the aformentioned manuscript for discussion on these specificities. If using the code provided here, or part thereof, for a scientific publication, please refer to this book for details on the implementation choices.

This library was developed with the direct help and indirect advice of several people. Many thanks to them: Daniel Anderson, Hanz Martin Cheng, Daniele Di Pietro, Lorenzo Botti.

The development of this library was partially supported by Australian Government through the Australian Research Council's Discovery Projects funding scheme (project number DP170100605).
