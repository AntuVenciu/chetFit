// Build the CHET detector
#include <vector>

#include "DetectorConstruction.hh" 
#include "TrackerSD.hh"

#include "G4VPhysicalVolume.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Colour.hh"
#include "G4Tubs.hh"
#include "G4VisAttributes.hh"
#include "G4UnionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SDManager.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4GDMLParser.hh"

namespace B1 {

  DetectorConstruction::DetectorConstruction() : 
    fNbOfPlanes(0), fLogicPlanes(NULL) {
    fNbOfPlanes = 7; // 30 petals + 4 cylinders or 7 cylinders, three instead of petals
    fLogicPlanes = new G4LogicalVolume*[fNbOfPlanes];
  }

  G4VPhysicalVolume* DetectorConstruction::Construct() {
    G4NistManager* nist = G4NistManager::Instance();

    // Define materials
    G4Material* scint_mat = nist->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
    G4Material* air_mat = nist->FindOrBuildMaterial("G4_AIR");
    
    // World
    G4double world_sizeZ = 105. * cm;
    G4double world_sizeXY = 30. * cm;
    G4Material* world_mat = nist->FindOrBuildMaterial("G4_Galactic");
    
    G4Box* solidWorld = new G4Box("World", world_sizeXY/2., world_sizeXY/2., world_sizeZ/2.);
    G4LogicalVolume* logicWorld = new G4LogicalVolume(solidWorld, world_mat, "World");
    G4VPhysicalVolume* physWorld = new G4PVPlacement(nullptr, G4ThreeVector(), logicWorld, "World", nullptr, false, 0, true);
    
    // Magnetic field
    ConstructSDandField(logicWorld);
    ConstructElectrodes(logicWorld);
    ConstructWFCoil(logicWorld);
    ConstructPulseFrame(logicWorld);
    
    // Create SD
    G4String detectorSDname = "chetSD";
    TrackerSD* aTrackerSD = new TrackerSD(detectorSDname, "TrackerHitsCollection");
    G4SDManager::GetSDMpointer()->AddNewDetector(aTrackerSD);

    // Detector planes
    
    const int numPlanes = 30;
    /*
    const double planeWidth = 0.1 * cm;
    const double planeLength = 4. * cm;
    const double planeHeightMin = 6.5 * cm;
    const double planeHeightMax = 9.5 * cm;
    const double deltaPhi = twopi / numPlanes;
    
    // Color of petals
    G4Colour blue(0., 0., 1., 0.3);
    G4VisAttributes* petalVisAttributes = new G4VisAttributes(blue);

    for (int i = 0; i < numPlanes; ++i) {
      double phi = i * deltaPhi;
      
      double centerX = std::cos(phi) * (planeHeightMin + (planeHeightMax - planeHeightMin) / 2);
      double centerY = std::sin(phi) * (planeHeightMin + (planeHeightMax - planeHeightMin) / 2);
      double centerZ = 0.;
      
      G4RotationMatrix* rot = new G4RotationMatrix();
      rot->rotateZ(phi);
      
      G4ThreeVector translation(centerX, centerY, centerZ);
      G4Transform3D transform(*rot, translation);
      
      G4Box* solidPlane = new G4Box("Plane", (planeHeightMax - planeHeightMin) / 2., planeWidth / 2., planeLength / 2.);
      fLogicPlanes[i] = new G4LogicalVolume(solidPlane, scint_mat, "Plane");
      fLogicPlanes[i]->SetSensitiveDetector(aTrackerSD);
      fLogicPlanes[i]->SetVisAttributes(petalVisAttributes);
      new G4PVPlacement(transform, fLogicPlanes[i], "Plane", logicWorld, false, i, true);
    }
    */

    bool use_disks = false;
    bool use_cylinders = true;

    if (use_disks) {
      // End Cap disks:
      // We place 8 disks US ans 8 disks DS
      // 4 disks are placed at radius 2 cm -> 2.5 cm;
      // 4 disks at radius 3.5 cm -> 4 cm
      // They are numbered like this: 0->7 US; 7->15 DS; even = inner, odd = outer 
      // If you change the diskRadialExtension and diskThickness you can make them into cylinders
      
      const int numDisks_perSide = 4; // Number of disks per side
      double diskInnerRadius = 2.0 * cm; // in cm
      double diskRadialExtension = 0.05 * cm; // in cm
      double diskOuterRadius = 6.0 * cm; // in cm
      double diskThickness = 6 * cm; // in cm (100 µm)
      double diskSpacing = 1 * cm; // Distance between disks along z-axis
      double zmin = 4.0 * cm; // cm
      double zmax = zmin + ((numDisks_perSide/2) - 1) * diskSpacing;
      
      for (int i = 0; i < numDisks_perSide; i++) {
	
	// calculate z position of disks
	double z = (i < numDisks_perSide/2.) ? - zmax + i * diskSpacing : zmin + (i - numDisks_perSide/2) * diskSpacing;
	
	// place inner and outer radii disks
	for (int j=0; j<2; j++) {
	  
	  int diskID = 2 * i + j + numPlanes; // ID of disks
	  G4Tubs* solidDisk;
	  
	  if (i%2 == 0) { // disk geometry at inner radius
	    solidDisk = new G4Tubs("Plane", diskInnerRadius, diskInnerRadius + diskRadialExtension, diskThickness / 2., 0., twopi);
	  }
	  else { // disk geometry at outer radius
	    solidDisk = new G4Tubs("Plane", diskOuterRadius - diskRadialExtension, diskOuterRadius, diskThickness / 2., 0., twopi);
	  }
	  
	  fLogicPlanes[diskID] = new G4LogicalVolume(solidDisk, scint_mat, "Plane");
	  fLogicPlanes[diskID]->SetSensitiveDetector(aTrackerSD);
	  G4ThreeVector position(0., 0., z);
	  new G4PVPlacement(0, position, fLogicPlanes[diskID], "Plane", logicWorld, false, diskID, true);
	  
	}
      }
    }

    else if (use_cylinders) {

      std::vector<double> radii = {2.05 * cm, 2.35 * cm, 3.60 * cm, 4.55 * cm, 6.55 * cm, 7.55 * cm, 8.55 * cm};
      double cylRadialExtension = 0.1 * cm; // in cm
      double cylThickness = 40 * cm; // in cm (100 µm)

      // Color of cylinders
      G4Colour orange(0.65, 0.35, 0., 0.3);
      G4VisAttributes* cylinderVisAttributes = new G4VisAttributes(orange);
      
      int cylID = 0;

      // place cylinder at various radii
      for (auto r : radii) {
	  
	double radiusIn = r - cylRadialExtension / 2.;
	double radiusOut = r + cylRadialExtension / 2.;
	G4Tubs* cylinder = new G4Tubs("Plane", radiusIn, radiusOut, cylThickness / 2., 0., twopi);

	// Add to logic plane
	fLogicPlanes[cylID] = new G4LogicalVolume(cylinder, scint_mat, "Plane");
	fLogicPlanes[cylID]->SetSensitiveDetector(aTrackerSD);
	fLogicPlanes[cylID]->SetVisAttributes(cylinderVisAttributes);
	G4ThreeVector position(0., 0., 0.);
	new G4PVPlacement(0, position, fLogicPlanes[cylID], "Plane", logicWorld, false, cylID, true);
	cylID++;
      }

    }

    // Write the Geometry to a GDML file
    G4GDMLParser parser;
    parser.Write("detectorGeometry_z20_fullGeo_7Cyl_Nopetals.gdml", physWorld);

    return physWorld;
  }
  
  
  // -------------------------------
  //  SD Field
  // -------------------------------
  void DetectorConstruction::ConstructSDandField(G4LogicalVolume *fLogicWorld) {
    // Create a uniform magnetic field of 2.9 Tesla oriented along the z-axis
    magneticField = new G4UniformMagField(G4ThreeVector(0, 0, 2.2 * tesla));
    
    // Create a field manager and set the magnetic field
    fieldManager = new G4FieldManager(magneticField);
    G4TransportationManager::GetTransportationManager()->GetFieldManager()->SetDetectorField(magneticField);
    G4TransportationManager::GetTransportationManager()->GetFieldManager()->CreateChordFinder(magneticField);

    // Aluminum cylinder (PSC magnet)
    G4NistManager* nist = G4NistManager::Instance();
    G4Material* aluminum_mat = nist->FindOrBuildMaterial("G4_Al");
    G4double aluminumInnerRadius = 10. * cm;
    G4double aluminumOuterRadius = (10. + 15.) * cm;
    G4double aluminumHalfLength = 50. * cm;
    
    // Color
    G4Colour grey(0.33, 0.33, 0.34, 0.1);
    G4VisAttributes* magnetVisAttributes = new G4VisAttributes(grey);

    G4Tubs* solidAluminumCylinder = new G4Tubs("AluminumCylinder", aluminumInnerRadius, aluminumOuterRadius, aluminumHalfLength, 0., twopi);
    G4LogicalVolume* logicAluminumCylinder = new G4LogicalVolume(solidAluminumCylinder, aluminum_mat, "AluminumCylinder");
    logicAluminumCylinder->SetVisAttributes(magnetVisAttributes);
    new G4PVPlacement(0, G4ThreeVector(), logicAluminumCylinder, "AluminumCylinder", fLogicWorld, false, 0);

  }

  // -------------------------------

    // ---------------------------------- //

  // -------- Electrodes ------------- //
  void DetectorConstruction::ConstructElectrodes(G4LogicalVolume *fLogicWorld) {
    G4NistManager* nist = G4NistManager::Instance();
    
    G4Element* fH  = new G4Element( "H",  "H", 1., 1.01*g/mole);
    G4Element* fC  = new G4Element( "C",  "C", 6., 12.01*g/mole);
    G4Material * fBC400 = new G4Material("BC400", 1.023*g/cm3, 2);
    fBC400->AddElement(fC, 1000);
    fBC400->AddElement(fH, 1103);

    G4double alpha = 0.3 + 0.4;
    G4Material* Al = nist->FindOrBuildMaterial("G4_Al");
    G4Material* Ka = nist->FindOrBuildMaterial("G4_KAPTON");
    G4Material* BC = fBC400;
    
    //==========================================================================================
    // Cathode
    G4double CathodeThickBC = 1.0 * mm;
    G4double CathodeThickAl = 0.00002 * mm;
    G4double CathodeThickKa = 0.02 * mm;
    G4double cathodeHeight  = 400 * mm;
    G4double cathodeIR      = 35 * mm;
    
    G4ThreeVector posCathode = G4ThreeVector(0.0*cm, 0.0*cm, 0.0*cm);
    G4RotationMatrix* rotCathode = new G4RotationMatrix();
    rotCathode->rotateX(0.0*deg);

    // Kapton Layer
    G4Tubs* cathodeKa = new G4Tubs("CathodeKa", cathodeIR, cathodeIR + CathodeThickKa, cathodeHeight / 2.0, 0*deg, 360*deg);
    G4LogicalVolume* logicCathodeKa = new G4LogicalVolume(cathodeKa, Ka, "CathodeKa_LV");

    G4VisAttributes* va = new G4VisAttributes(G4Color(0.5, 0.5, 0.5, 0.5));
    va->SetForceSolid(true);
    logicCathodeKa->SetVisAttributes(va);

    new G4PVPlacement(rotCathode, posCathode, logicCathodeKa, "CathodeKa_PV", fLogicWorld, false, 0, true);

    // Al Layer
    cathodeIR += CathodeThickKa;
    G4Tubs* cathodeAl = new G4Tubs("CathodeAl", cathodeIR, cathodeIR + CathodeThickAl, cathodeHeight / 2.0, 0*deg, 360*deg);
    G4LogicalVolume* logicCathodeAl = new G4LogicalVolume(cathodeAl, Al, "CathodeAl_LV");

    logicCathodeAl->SetVisAttributes(va);

    new G4PVPlacement(rotCathode, posCathode, logicCathodeAl, "CathodeAl_PV", fLogicWorld, false, 0, true);


    //==========================================================================================
    //Anode
    G4double AnodeThickAl = 0.00002 * mm;
    G4double AnodeThickKa = 0.02 * mm;
    G4double AnodeHeight  = 400 * mm;
    G4double AnodeIR      = 25 * mm;

    G4ThreeVector posAnode = G4ThreeVector(0.0*cm, 0.0*cm, 0.0*cm);
    G4RotationMatrix* rotAnode = new G4RotationMatrix();
    rotAnode->rotateX(0. * deg);
    
    // Al Layer
    G4Tubs* AnodeAl = new G4Tubs("AnodeAl", AnodeIR, AnodeIR + AnodeThickAl, AnodeHeight / 2.0, 0*deg, 360*deg);
    G4LogicalVolume* logicAnodeAl = new G4LogicalVolume(AnodeAl, Al, "AnodeAl_LV");

    logicAnodeAl->SetVisAttributes(va);

    new G4PVPlacement(rotAnode, posAnode, logicAnodeAl, "AnodeAl_PV", fLogicWorld, false, 0, true);    

    // Kapton Layer
    AnodeIR += AnodeThickAl;
    G4Tubs* AnodeKa = new G4Tubs("AnodeKa", AnodeIR, AnodeIR + AnodeThickKa, AnodeHeight / 2.0, 0*deg, 360*deg);
    G4LogicalVolume* logicAnodeKa = new G4LogicalVolume(AnodeKa, Ka, "AnodeKa_LV");

    logicAnodeKa->SetVisAttributes(va);

    new G4PVPlacement(rotAnode, posAnode, logicAnodeKa, "AnodeKa_PV", fLogicWorld, false, 0, true);
  }
  // ------------------------------------ //

  // ---------- WF Coils --------------- //
  // Weakly focusing coil
  void DetectorConstruction::ConstructWFCoil(G4LogicalVolume *fLogicWorld) {
    G4NistManager* nist = G4NistManager::Instance();
    G4Material* WfCoilMaterial = nist->FindOrBuildMaterial("G4_Cu");
    G4ThreeVector posWfCoil = G4ThreeVector(0*cm, 0*cm, 0*cm);
    G4RotationMatrix* rotWfCoil = new G4RotationMatrix();
    /* rotWfCoil->rotateX(0.0*deg); */
    /* rotWfCoil->rotateY(0.0*deg); */

    G4double WfCoilHeight = 10 * mm;
    G4double WfCoilRadius = 60 * mm ;
    G4double WfCoilInnerRadius = 50 * mm;
    /* G4Cons* WfCoil = new G4Cons("WfCoil", WfCoilInnerRadius, WfCoilRadius, WfCoilInnerRadius, WfCoilRadius, WfCoilHeight / 2.0, 0*deg, 360*deg); */
    G4Tubs* WfCoil = new G4Tubs("WfCoil", WfCoilInnerRadius, WfCoilRadius, WfCoilHeight / 2.0, 0*deg, 360*deg);
    G4LogicalVolume* logicWfCoil = new G4LogicalVolume(WfCoil, WfCoilMaterial, "WfCoil");
    new G4PVPlacement(rotWfCoil, posWfCoil, logicWfCoil, "WfCoil", fLogicWorld, false, 0, true);

    G4VisAttributes* va = new G4VisAttributes(G4Color(0.722, 0.451, 0.2, 1.0));
    va->SetForceSolid(true);
    logicWfCoil->SetVisAttributes(va);
  }
  //  -------------------------------------- //

  // -------- Kicker Coils ----------------- //
  // Kicker coil
  void DetectorConstruction::ConstructPulseFrame(G4LogicalVolume *fLogicWorld) {
    G4NistManager* nist = G4NistManager::Instance();
    G4Material* FrameMat = nist->FindOrBuildMaterial("G4_Cu");
    G4Material* InStriveMat = nist->FindOrBuildMaterial("G4_MYLAR");
    G4double Layer2thickness = 1.0 * mm;
    G4double Layer3thickness = 0.1 * mm;
    G4double alpha = 0.3;

    G4double FrameOuterRadius = 41.12002 * mm;
    G4double FrameInnerRadius = 41.02002 * mm;
    G4double FrameSizeZ       = 65.0 * mm;
    G4double FramePositionX   = 0. * mm;
    G4double FramePositionY   = 0. * mm;
    G4double FramePositionZ   = 0. * mm;
    G4double RadiusL2 = FrameOuterRadius + Layer2thickness;
    G4double RadiusL3 = RadiusL2 + Layer3thickness;

    
    G4double StriveSizeXY    = 20.0 * mm;
    G4double StriveSizeZ     = 25.0 * mm;    
    G4double StriveAngle     = (StriveSizeXY/FrameOuterRadius) * 180/pi;

    G4VisAttributes* vaCu = new G4VisAttributes(G4Color(0.6, 0.4, 0.0, alpha));
    vaCu->SetForceSolid(true);
    G4VisAttributes* vaMy = new G4VisAttributes(G4Color(0.8, 8, 0.8, alpha-0.1));
    vaMy->SetForceSolid(true);

    G4double rotangle = 0;
    G4RotationMatrix* rotation = new G4RotationMatrix();
    rotation->rotateZ(rotangle * deg);

    // Frame
    G4Tubs* FrameUp    = new G4Tubs("FrameUp"  , FrameInnerRadius, FrameOuterRadius, (FrameSizeZ - StriveSizeZ) / 4.0, 0*deg, 360*deg);
    G4Tubs* FrameDown  = new G4Tubs("FrameDown", FrameInnerRadius, FrameOuterRadius, (FrameSizeZ - StriveSizeZ) / 4.0, 0*deg, 360*deg);
    
    G4LogicalVolume* FrameUp_LV = new G4LogicalVolume(FrameUp, FrameMat, "FrameUp_LV");
    G4LogicalVolume* FrameDown_LV = new G4LogicalVolume(FrameDown, FrameMat, "FrameDown_LV");

    FrameUp_LV->SetVisAttributes(vaCu);
    FrameDown_LV->SetVisAttributes(vaCu);

    G4ThreeVector posFUp = G4ThreeVector(0,0, StriveSizeZ / 2.0 + (FrameSizeZ - StriveSizeZ) / 4.0);
    new G4PVPlacement(rotation,                 // rotation
                      posFUp,                   // at (x,y,z)
                      FrameUp_LV,               // its logical volume
                      "FrameUp_PV",             // its name
                      fLogicWorld,              // its mother  volume
                      false,                    // no boolean operations
                      0,                        // copy number
                      true);                    // checking overlaps
    G4ThreeVector posFDown = G4ThreeVector(0,0, -StriveSizeZ / 2.0 - (FrameSizeZ - StriveSizeZ) / 4.0);
    new G4PVPlacement(rotation,                 // rotation
                      posFDown,                 // at (x,y,z)
                      FrameDown_LV,             // its logical volume
                      "FrameDown_PV",           // its name
                      fLogicWorld,              // its mother  volume
                      false,                    // no boolean operations
                      0,                        // copy number
                      true);                    // checking overlaps

    // Strives
    G4Tubs* Strive1L1  = new G4Tubs("Strive1L1", FrameInnerRadius, FrameOuterRadius, StriveSizeZ / 2.0, (StriveAngle/2.0)*deg, StriveAngle*deg);
    G4Tubs* Strive2L1  = new G4Tubs("Strive2L1", FrameInnerRadius, FrameOuterRadius, StriveSizeZ / 2.0, (StriveAngle/2.0 + 90)*deg, StriveAngle*deg);
    G4Tubs* Strive3L1  = new G4Tubs("Strive3L1", FrameInnerRadius, FrameOuterRadius, StriveSizeZ / 2.0, (StriveAngle/2.0 + 180)*deg, StriveAngle*deg);
    G4Tubs* Strive4L1  = new G4Tubs("Strive4L1", FrameInnerRadius, FrameOuterRadius, StriveSizeZ / 2.0, (StriveAngle/2.0 + 270)*deg, StriveAngle*deg);
    G4Tubs* Strive1L2  = new G4Tubs("Strive1L2", FrameOuterRadius, RadiusL2, FrameSizeZ / 2.0, (StriveAngle/2.0)*deg, StriveAngle*deg);
    G4Tubs* Strive2L2  = new G4Tubs("Strive2L2", FrameOuterRadius, RadiusL2, FrameSizeZ / 2.0, (StriveAngle/2.0 + 90)*deg, StriveAngle*deg);
    G4Tubs* Strive3L2  = new G4Tubs("Strive3L2", FrameOuterRadius, RadiusL2, FrameSizeZ / 2.0, (StriveAngle/2.0 + 180)*deg, StriveAngle*deg);
    G4Tubs* Strive4L2  = new G4Tubs("Strive4L2", FrameOuterRadius, RadiusL2, FrameSizeZ / 2.0, (StriveAngle/2.0 + 270)*deg, StriveAngle*deg);
    G4Tubs* Strive1L3  = new G4Tubs("Strive1L3", RadiusL2, RadiusL3, FrameSizeZ / 2.0, (StriveAngle/2.0)*deg, StriveAngle*deg);
    G4Tubs* Strive2L3  = new G4Tubs("Strive2L3", RadiusL2, RadiusL3, FrameSizeZ / 2.0, (StriveAngle/2.0 + 90)*deg, StriveAngle*deg);
    G4Tubs* Strive3L3  = new G4Tubs("Strive3L3", RadiusL2, RadiusL3, FrameSizeZ / 2.0, (StriveAngle/2.0 + 180)*deg, StriveAngle*deg);
    G4Tubs* Strive4L3  = new G4Tubs("Strive4L3", RadiusL2, RadiusL3, FrameSizeZ / 2.0, (StriveAngle/2.0 + 270)*deg, StriveAngle*deg);

    G4UnionSolid* CuStrives = new G4UnionSolid("CuStrives",Strive1L1,Strive2L1);
    CuStrives = new G4UnionSolid("CuStrives",CuStrives,Strive3L1);
    CuStrives = new G4UnionSolid("CuStrives",CuStrives,Strive4L1);
    CuStrives = new G4UnionSolid("CuStrives",CuStrives,Strive1L3);
    CuStrives = new G4UnionSolid("CuStrives",CuStrives,Strive2L3);
    CuStrives = new G4UnionSolid("CuStrives",CuStrives,Strive3L3);
    CuStrives = new G4UnionSolid("CuStrives",CuStrives,Strive4L3);

    G4UnionSolid* MyStrives = new G4UnionSolid("MyStrives",Strive1L2,Strive2L2);
    MyStrives = new G4UnionSolid("MyStrives",MyStrives,Strive3L2);
    MyStrives = new G4UnionSolid("MyStrives",MyStrives,Strive4L2);

    G4LogicalVolume* CuStrives_LV = new G4LogicalVolume(CuStrives, FrameMat, "CuStrives_LV");
    G4LogicalVolume* MyStrives_LV = new G4LogicalVolume(MyStrives, InStriveMat, "MyStrives_LV");

    CuStrives_LV->SetVisAttributes(vaCu);
    MyStrives_LV->SetVisAttributes(vaMy);

    G4ThreeVector posStrives = G4ThreeVector(0,0,0);
    new G4PVPlacement(rotation,                 // rotation
                      posStrives,               // at (x,y,z)
                      CuStrives_LV,             // its logical volume
                      "CuStrives_PV",           // its name
                      fLogicWorld,              // its mother  volume
                      false,                    // no boolean operations
                      0,                        // copy number
                      true);                    // checking overlaps
    new G4PVPlacement(rotation,                 // rotation
                      posStrives,               // at (x,y,z)
                      MyStrives_LV,             // its logical volume
                      "MyStrives_PV",           // its name
                      fLogicWorld,              // its mother  volume
                      false,                    // no boolean operations
                      0,                        // copy number
                      true);                    // checking overlaps
  }
  //=================================================================

}
