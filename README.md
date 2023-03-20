------------------------------------------------------------------------

 sFFG4MC: a Geant4 framework for simulating the JLab sFF experiment

------------------------------------------------------------------------

	http://github.com/dhamilton-glasgow/sFFG4MC
	(developed and tested with Geant4.10.4 and root6.22 on centos7)

	Scattering chamber, exit beamline and target 
	modified from original NPS simulation:
  	https://github.com/gboon18/HallC_NPS

------------------------------------------------------------------------
 Compilation

	mkdir build
  	cd build
  	cmake ../ (or cmake3 ../ on some systems)
  	make -jX (X specifies number of cores to use)

------------------------------------------------------------------------
 Running the simulation in visualisation mode

  	./sFFG4MC
  	then in the gui: 
	/control/execute macros/vis_beam.mac 
	or /control/execute macros/vis_root.mac 

------------------------------------------------------------------------
 Running the simulation in batch mode

  	./sFFG4MC macros/batch_beam.mac 
  	or 
  	./sFFG4MC macros/batch_root.mac 

------------------------------------------------------------------------
 DetectorConstruction options (always call update following any changes)

  	/sFFG4MC/detector/setShieldThickness  
  	Set Harm lead shield thickness  (value and unit)

  	/sFFG4MC/detector/setWindowThickness
  	Set scattering chamber window thickness (value and unit)

  	/sFFG4MC/detector/update	 
  	Update the detector geometry with changed values
  	Must be run before beamOn if detector has been changed  

------------------------------------------------------------------------
 PhysicsList options (always call in macro before /run/initialize)

  	/sFFG4MC/physics/addPhysics 
  	Add physics list (standard_opt3, QGSP_BIC_EMY, QGSP_BIS_HP, QGSP_BERT_HP, FTFP_BERT_HP)

------------------------------------------------------------------------
 Generator options 

  	/sFFG4MC/generator/Mode 
  	Set the mode of the generator (0 for BEAM or 1 for ROOT)

  	BEAM mode using geant4 gps to generate a 6.6 GeV electron beam with a 2x2mm raster.
  	ROOT mode requires a root file with a tree called TGen, which contains branches:

		Nparticles	-- number of primary particles in this event
  		weight		-- event weight (for consistency this should be per uA beam current)
		flag		-- an integer to represent generator (reaction) type, ie beam=0, elastic=1, ...
		vx[Nparticles]	-- vertex position vector
		vy[Nparticles]
		vz[Nparticles]
		px[Nparticles]	-- momentum direction vector
		py[Nparticles]
		pz[Nparticles]
		E[Nparticles]	-- energy in MeV
		pdg[Nparticles]	-- pdg code

	/sFFG4MC/generator/Nevents
	Set the number of primary events to be generated

  	/sFFG4MC/generator/InputFile
  	Set the full name and path of the file with the input ROOT ntuple (in ROOT mode)  

------------------------------------------------------------------------
 OutputManager options 

  	/sFFG4MC/output/setOutputFile
  	Set the full name and path of the output file

  	The output tree is called TOut and has the following branches:
		Event_weight	-- event weight (for consistency this should be per uA beam current)
		Event_flag	-- an integer to represent generator (reaction) type, ie beam=0, elastic=1, ...
		Primary_Nhits	-- number of primary particles in this event
		Primary_*[]	-- arrays of primary variables, including pdg, energy, position, direction
		Virtual_Nhits	-- number of virtual detector hits in this event
		Virtual_*[]	-- arrays of virtual hit variables, including pdg, energy, position, direction,
				   particle vertex, time, track id, mother track id and detector id variables:
					Virtual_det (0 for Earm, 1 for hodoscope, 2 for Harm)
					Virtual_mod (sector between 0 and 5)
					Virtual_row and Virtual_col 
		Real_Nhits	-- number of real detector hits in this event
		Real_*[]	-- arrays of real hit variables, including energy deposit, position, time 
				   and detector id variables: 
					Real_det (0 for Earm, 1 for hodoscope, 2 for Harm)
					Real_mod (sector between 0 and 5)
					Real_row and Real_col


