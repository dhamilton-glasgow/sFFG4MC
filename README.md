------------------------------------------------------------------------

 sFFG4MC: a Geant4 framework for simulating the JLab sFF experiment

------------------------------------------------------------------------

	http://github.com/dhamilton-glasgow/sFFG4MC
	(developed and tested with Geant4.10.4 and root6.22 on centos7)

	Scattering chamber, exit beamline and target modified from original NPS simulation:
  	https://github.com/gboon18/HallC_NPS

	Default geometry:
	Target is 10cm LH2
	Scattering chamber window thickness is 0.05 cm
	Earm is an array of 2 x 2 x 20 cm3 PbW04 blocks in six sectors (960 total)
	Earm is at 15.5 degrees and a distance of 239 cm
	Harm is an array of 15 x 15 x 101.2 cm3 scint-Fe blocks (288 total)
	Harm blocks are 44 pairs of plastic scintillator (10mm) and iron (13mm)
	Harm is at 42.5 degrees and a distance of 351 cm
	Hodoscope is an plastic scintillator array of 3 x 3 x 10 cm3 detectors (7200 total)
	Harm has a 10 cm thick lead shield

	Default units: MeV, cm, ns

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

	/sFFG4MC/detector/setBeamline
  	Set an integer flag for whether to include scattering chamber and beamline (0 or 1)

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
		flag		-- an integer to represent generator (reaction) type, ie elastic=1
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
		Event_flag	-- an integer to represent generator (reaction) type, ie elastic=1
		Primary_Nhits	-- number of primary particles in this event
		Primary_*[]	-- arrays of primary variables, including pdg, energy, position, direction
		Virtual_Nhits	-- number of virtual detector hits in this event
		Virtual_*[]	-- arrays of virtual hit variables, including pdg, energy, position, direction,
				   particle vertex, time, track id, mother track id and detector id variables:
					Virtual_det[] -- which detector (0 for Earm, 1 for hodoscope, 2 for Harm)
					Virtual_mod[] -- which sector (0 to 5)
					Virtual_row[] -- which row (0-31 Earm, 0-79 hodoscope, 0-16 Harm)
					Virtual_col[] -- which col (0-5 Earm, 0-15 hodoscope, 0-3 Harm)
		Real_Nhits	-- number of real detector hits in this event
		Real_*[]	-- arrays of real hit variables, including energy deposit, position, time,
				   and detector id variables:
					Real_det[] -- which detector (0 for Earm, 1 for hodoscope, 2 for Harm)
					Real_mod[] -- which sector (0 to 5)
					Real_row[] -- which row (0-31 Earm, 0-79 hodoscope, 0-16 Harm)
					Real_col[] -- which col (0-5 Earm, 0-15 hodoscope, 0-3 Harm)



