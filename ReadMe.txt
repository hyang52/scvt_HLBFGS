About:
	SCVT-HLBFGS is a Parallel SCVT Generator written in C++, using the Limited-memory BFGS for optimization.
	It makes use of boost's mpi and serialization libraries.
	It was written by Huanhuan Yang, while the Overlapping domain decomposition method for parallelization in 
	scvt computation is following the scvt-mpi (Lloyd method) program by Doug Jacobsen and Geoff Womeldorff

Compiling:
	To compile, you need the boost mpi and serialization libraries
	to be in your LD_LIBRARY_PATH, and the header files in your PATH variables.
	Afterwards, 
		make							build scvt_HLBFGS.exe
		make partition_gen				build partition_gen.exe			(check density function before you build)
		make num_estimate				build num_estimate.exe

Boost:
	Boost can be downloaded from:
	http://www.boost.org/

	The only libraries that need to be setup at mpi, and serialization. Please refer to boost's documentation
	to set up correctly.

Running:
	To run, simply type:
		./scvt_HLBFGS.exe
	
	To run using mpi, simply type:
		mpirun -n P ./scvt_HLBFGS.exe
	Where P is replaced with the number of processors used. The number of Processors has to be either
	1 or the number of points in RegionList as described below.

	run: ./partition_gen -p 16 
		 to create the partition map that will be used for 16 processors

Refinement ex: 12 42 162 642 2562 10242 40,962 163,842 655,362 2,621,442	---	 num_pts after refinement

Input:
	parameters: (Required, load automatically)
		This file contains parameters required to run the program 

	RegionList: (Required, load automatically from the subfolder region/)
		This file contains a list of region centers, for parallelization. Some example lists are provided in the
		original source package. Each region center needs to be given as x, y, z with one region per line.
		One region should be provided for each processor to be used. A simple example two processor setup (the 
		smallest number MPI-SCVT is capable of using) is the north and south pole, with a RegionList file containing
		the following data:
			0 0 1
			0 0 -1

	RegionTriangulation: (Required, load automatically from the subfolder region/)
		This file contains the connectivity of the regions in RegionList. Each line contains one "triangle" given as 
		three indices, where each index relates to a point number from the RegionList file. Contradictory to the name, 
		the connectivity does not need to be a Triangulation, though a Triangulation is easy to come by. As an example, 
		the connectivity for the example provided for RegionList would be:
			0 1 0

		A triangulation (in triangles.dat) as output from scvt_HLBFGS.exe can be used as well

	SaveVertices: (Optional, if use by config option then load automatically)
		This file contains an initial point set for iterating on. These points should be given as x, y, z
		with one point per line. This is required if the method of generating points is set to 0 in Params.

	SaveBoundaries: (not considered at this moment)
		This file lists all points which lie along some boundary segment or segments. Generator points can
		be projected onto these segments during generation, but do not have to be.

	SaveLoopCounts: (not considered at this moment)
		This file defines the continous segements of boundary lines. The first
		column is a beginning index, and the second column is the number of
		points to read from SaveBoundaries to create a loop.

Output:
	After a successful run of scvt_HLBFGS, files are created.

	point_initial.dat:
		This file contains initial generator pointset read from file or created online, given as x, y, z with one point per line

	end_points.dat:
		This file contains the ending pointset from MPI-SCVT, given as x, y, z with one point per line

	point_density.dat:
		This file contains the density function values of the ending pointset.

	boundary_points.dat: (not considered at this moment)
		This file contains all the boundary points: including bdry pts read from file and bdry points possiblely plugged in online to reach certain boundary resolution.

	triangles.dat:
		This file contains the final triangulation of the points in end_points.dat, given as one triangle
		per line, with each number being the index of a point in end_points.dat.

	bisected_points.dat:
		This file contains the bisected pointset from end_points.dat. For a full sphere grid, it should have
		4*n - 6 points, where n is the number of points in end_points.dat. 

	if save_before_bisect: end_points.dat.12	triangles.dat.12	point_density.dat.12 	(for instance)

	if printRegions:	RegionCenters.dat	  RegionRadii.dat	   RegionProperties.dat (local number of points)


