# Problem sets
Example problems (without solutions) in pdf format.  You can download problem sets by chapter, or as a full set.   Solutions to the problems will be published as a companion volume to the second edition of Applied Mechanics of Solids

# Example Finite Element Analysis codes
Example MATLAB codes from the text "Applied Mechanics of Solids" 2nd edition, Allan Bower

1. FEM_conststrain_triangles.m  Demonstration static linear elastic code using constant strain triangles.   The code can be run with input files FEM_conststrain_input.txt or FEM_conststrain_holeplate.txt;

2. FEM_1D_Static.m  Simple static analysis of a 1D bar subjected to axial body force;

3. FEM_1D_newmark.m Simple dynamic analysis of a 1D bar subjected to axial body force, using Newmark time integration.

4. FEM_1D_modal.m Simple dynamic analysis of a 1D bar subjected to axial body force, using modal time integration.

5. The following files all solve 2D or 3D static linear elastic problems, but illustrate various refinements of the finite element method:

		(a) FEM_2Dor3D_linelast_standard.m 2D(plane strain/stress) or 3D static linear elasticity code with fully integrated elements.  The code can be run with the following input files.
		Linear_elastic_triangles.txt: 2D plane strain problem with two triangular element; 
		Linear_elastic_quad4.txt: 2D plane strain problem with 4 noded quadrilateral element;
		Linear_elastic_quad8.txt: 2D plane strain problem with 8 noded quadrilateral element;
		Linear_elastic_brick4.txt: 3D problem with 4 noded brick element;
		Linear_elastic_brick20.txt: 3D problem with 20 noded brick element
		Linear_elastic_pressurized_cylinder.txt: 2D simulation of a pressurized cylinder;
		(b) FEM_shear_locking_demo.m - Solves the beam bending problem discussed in Section 8.6.2, and compares the FEM solution with the exact solution to illustrate shear locking.
			  This version of the code must be run with shear_locking_demo_linear.txt (solution with 4 noded quad elements)
		(c) FEM_incompatible_modes.m  - Solves the beam bending problem discussed in Section 8.6.2 using incompatible mode elements, and compares the FEM solution with the exact solution to demonstrate that the elements avoid shear locking.
			  This version of the code must be run with shear_locking_demo_linear.txt (solution with 4 noded quad elements)
		(d) FEM_volumetric_locking_demo.m  - Solves the pressurized cylindrical cavity problem discussed in Section 8.6.2, and compares the FEM solution with the exact solution.
			  This version of the code must be run with volumetric_locking_demo_linear.txt (solution with 4 noded quad elements) or volumetric_locking_demo_quadratic.txt (solution with 8 noded quadrilateral elements)
		(e) FEM_hourglassing_demo.m  - Solves the pressurized cylindrical cavity problem discussed in Section 8.6.2 with reduced integration elements, demonstrating hourglassing.
			  This version of the code must be run with volumetric_locking_demo_linear.txt (solution with 4 noded quad elements) 
		(f) FEM_selective_reduced_integration.m  - Solves the pressurized cylindrical cavity problem discussed in Section 8.6.2 using selectively reduced integration, and compares the FEM solution with the exact solution.
			  This version of the code must be run with volumetric_locking_demo_quadratic.txt (solution with 8 noded quadrilateral elements)
		(g) FEM_hourglasscontrol.m  - Solves the pressurized cylindrical cavity problem using reduced integration elements with hourglass control, and compares the FEM solution with the exact solution.
			  This version of the code must be run with volumetric_locking_demo_linear.txt (solution with 4 noded quad elements)
		(h) FEM_Bbar.m – Solves the pressurized cylinder problem discussed in Section 8.6.2 using the B-bar method, and compares the solution with the exact solution.
			  This version of the code must be run with volumetric_locking_demo_linear.txt or volumetric_locking_demo_quadratic.txt.
		(i) FEM_hybrid.m – Solves the pressurized cylinder problem discussed in Section 8.6.2 using hybrid elements, and compares the FEM solution with the exact solution.
			  This version of the code must be run with volumetric_locking_demo_linear.txt or volumetric_locking_demo_quadratic.txt

6. FEM_2Dor3D_linelast_dynamic.m:  Solves 2D or 3D dynamic linear elasticity problems, using Newmark time integration.   The code can be run with the input file Linear_elastic_dynamic_beam.txt.

7. FEM_2Dor3D_modeshapes.m:  Calculates mode shapes and natural frequencies for a linear elastic solid.   The code can be run with the input file Linear_elastic_dynamic_beam.txt.

8. FEM_hypoelastic.m: Solves 2D (plane strain only) or 3D static problems for a hypoelastic material, as discussed in Section 8.3.9.  The input files are hypoelastic_quad4.txt (2D plane strain) and hypoelastic_brick8.txt (3D).   

9. FEM_hyperelastic.m: Solves 2D (plane strain only) or 3D static problems for a hyperelastic (Neo-Hookean) material with fully integrated elements. An input file that calculates the deformed shape of a pressurized hyperelastic cylinder is provided in hyperelastic_cylinder.txt.  The code can also be run with the input file Bbar_hpere_demo.txt, which solves the same problem but for a near-incompressible material.   The solution is incorrect because of volumetric locking.

10. FEM_bbar_hyperelastic.m: Solves 2D (plane strain only) or 3D static problems for a hyperelastic (Neo-Hookean) material with locking-resistant ‘B-bar’ elements. The code can be run with the input file Bbar_hpere_demo.txt, which demonstrates that the B-bar method successfully prevents locking.

11. FEM_viscoplastic.m: Solves 2D (plane strain only) or 3D static problems for a small strain viscoplastic material. Input files for single-element tests are provided in viscoplastic_quad4.txt (2D plane strain) and viscoplastic_brick8.txt (3D)

12. FEM_truss.m; FEM_truss_hyperelastic.m; FEM_truss_plastic.m illustrate truss elements for linear elastic, hyperelastic, and viscoplastic materials, respectively.   The codes should be run with their corresponding input files FEM_truss_input.txt; FEM_truss_hyperelastic_input.txt and FEM_truss_plastic.txt.

13. FEM_beam_Euler.m and FEM_beam_Timoshenko.m calculate deflections in elastic beams with small deflections, using the Euler-Bernoulli and Timoshenko formulastions, respectively.   They can be run with input files FEM_beam_euler_input.txt and FEM_beam_timoshenko_input.txt.

14. FEM_plate_Kirchhoff.m and FEM_plate_Mindlin.m  calculate deflections in an elastic plate that experiences small deflections, using the Kirchhoff and Mindlin formulations, respectively..   The codes are set up to calculate the deflection of a circular plate loaded by a uniform pressure on its face.   No input files are needed as the mesh is generated directly by the code.

15. FEM_Penaltyconstraint.m and FEM_LagrangeMultiplierconstraint illustrates the use of the penalty method and Lagrange multipliers to enforce constraints in static computations.  They can be run with input files Penalty_example.txt and LagrangeMultiplier_example.txt

16. FEM_dynamic_constraints illustrates the predictor-corrector Newmark timestep for an explicit dynamic simulation with Lagrange multipliers.   It should be run with the input file dynamic_constrained_beam.txt

17. FEM_Cohesive_Zones.m illustrates the use of cohesive zone elements to model interface fracture.   It should be run with input files CohesiveZone_example.txt and CohesiveZone_example_2D.txt

18. FEM_Contact.m illustrates 2D contact elements.   It should be run with the input file Contact_example.txt 
