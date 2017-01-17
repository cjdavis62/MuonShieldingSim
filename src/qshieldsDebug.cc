#include "qshieldsDebug.hh"
#include "qshieldsDetectorConstruction.hh"
#include "qshieldsConfigurator.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4VisExtent.hh"
#include "G4VSolid.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Transform3D.hh"
#include "G4String.hh"
#include "Randomize.hh"

#include "TGraph2D.h"
#include "TFile.h"

#include <vector>
#include <iomanip>
#include <iostream>
#include <fstream>

qshieldsDebug::qshieldsDebug()
{

    theC = qshieldsConfigurator::GetInstance();
    theD = theC->CDet;

    volumes   = theD->fPhysicalVolumeVector;
    transform = theD->fTransformVector;

    if(volumes.size() != transform.size())
    {
        G4cout << "fPhysicalVolumeVector size doesn't match fTransformVector size!" << G4endl;
        return;
    }
    nEvents = theC->DebugEvents;

    v1 = theC->v1;
    v2 = theC->v2;

	vNames[1] = "PE_Shield";
	vNames[2] = "Top_PE_hield";
	vNames[3] = "Boric_Acid_Shield";
	vNames[4] = "External_Pb_Shield";
	vNames[5] = "300K_Vessel";
	vNames[6] = "300K_Flange";
	vNames[7] = "40K_Vessel";
	vNames[8] = "40K_Flange";
	vNames[9] = "4K_Vessel";
	vNames[10] = "4K_Flange";
	vNames[11] = "600mK_Vessel";
	vNames[12] = "600mK_Flange";
	vNames[13] = "Roman_Pb_Shield";
	vNames[14] = "Roman_Pb_Cradle";
	vNames[15] = "Roman_Pb_Cradle_Steel_Rods";
	vNames[16] = "50mK_Vessel";
	vNames[17] = "50mK_Flange";
	vNames[18] = "10mK_Vessel";
	vNames[19] = "10mK_Flange";
	vNames[20] = "Pb_Plug_Top_Plate";
	vNames[21] = "Pb_Plug";
	vNames[22] = "10mK_Tiles";
	vNames[24] = "Pb_Plug_Bottom_plate";
	vNames[25] = "TSP";
	vNames[26] = "Detector_Top_Plate";
	vNames[28] = "Detector_Bottom_Plate";
	vNames[30] = "Internal_DCS_Tubes";
	vNames[31] = "Internal_Source_String_Teflon";
	vNames[32] = "Internal_Source_String_Copper";
	vNames[34] = "Internal_Source_String_Teflon";
	vNames[35] = "Internal_Source_String_Copper_Int";
	vNames[37] = "CUORE_Detector";
	vNames[38] = "Detector_Cradle";
	vNames[39] = "Wire_Trays";
	vNames[40] = "Wire_Pads";
	vNames[41] = "NTDs";
	vNames[42] = "PTFE_Supports";
	vNames[43] = "CuBe_Joints";
	vNames[44] = "PEN_Tapes";
	vNames[46] = "SuperInsulation";
	vNames[54] = "Internal_DCS_Weights";
	vNames[55] = "Internal_DCS_Capsules_21";
	vNames[56] = "Internal_DCS_Capsules_4";
	vNames[57] = "Internal_DCS_Capsules_1";
	vNames[58] = "External_DCS_Weight";
	vNames[59] = "External_DCS_Capsules_20";
	vNames[60] = "External_DCS_Capsules_5";

    if(v1 != 0 && v2 != 0)
        CheckSomeVolumes();
    else
        CheckAllVolumes();

}

void qshieldsDebug::CheckSomeVolumes()
{
   const int vSize = volumes.size();

    char dfName[256];
    sprintf(dfName, "Debug_%s_%s.root", vNames[v1], vNames[v2]);

    TFile* debugFile = new TFile(dfName, "RECREATE");

    TGraph2D* v1_h  = new TGraph2D();
    TGraph2D* v2_h  = new TGraph2D();
    TGraph2D* int_h = new TGraph2D();

    int counter1 = 0, counter2 = 0, counterInt = 0;

    for(int i = 0; i < vSize; i++)
    {
        int firstVolume = CheckThisVolume(volumes[i]);
        if(!firstVolume)
            continue;

        G4VSolid*   firstSolid = volumes[i]->GetLogicalVolume()->GetSolid();
        G4VisExtent firstVExt  = firstSolid->GetExtent();

        G4RotationMatrix firstRot = transform[i].getRotation   ();
        G4ThreeVector    firstPos = transform[i].getTranslation();
    
        for(int event = 0; event < nEvents; event++)
        {
            if(event % (nEvents/10) == 0)
                G4cout << "Checking volume " << std::setw(15) << volumes[i]->GetName() << " ... " << std::setw(8) << event << "/" << nEvents << '\xd' << std::flush;

            Position = GetPointInVolume(firstVExt);
            while( firstSolid->Inside(Position) != kInside)
                Position = GetPointInVolume(firstVExt);			

            // Position in hall frame of reference
            Position *= firstRot;
            Position += firstPos;

            if     (firstVolume == v1) v1_h->SetPoint(counter1++, Position[0], Position[1], Position[2]);
            else if(firstVolume == v2) v2_h->SetPoint(counter2++, Position[0], Position[1], Position[2]);

            for(int j = 0; j < vSize; j++)
            {
                int secondVolume = CheckThisVolume(volumes[i]);
                if(!secondVolume)
                    continue;

                if( j == i )   continue;

                G4VSolid*   secondSolid  = volumes[j]->GetLogicalVolume()->GetSolid();
		
                G4RotationMatrix secondRot = transform[j].getRotation   ();
                G4ThreeVector    secondPos = transform[j].getTranslation();

                secondRot = (G4RotationMatrix) secondRot.inverse();

                G4ThreeVector newPosition = Position;			

                newPosition -= secondPos;
                newPosition *= secondRot;

                if( secondSolid->Inside(newPosition) == kInside)
                    int_h->SetPoint(counterInt++, Position[0], Position[1], Position[2]);
            }
        }
        G4cout << "Checking volume " << std::setw(15) << volumes[i]->GetName() << " ..." << std::setw(9) << nEvents 
               << "/" << nEvents << " DONE" << G4endl;
    }
    G4cout << G4endl;
    v1_h->Write(vNames[v1]);
    v2_h->Write(vNames[v2]);

    char intName[20];
    sprintf(intName, "%s_%s", vNames[v1], vNames[v2]);
    int_h->Write(intName);

    debugFile->Close();
    
    G4cout << "Created file " << dfName << G4endl << G4endl;
}

void qshieldsDebug::CheckAllVolumes()
{
    const int vSize   = volumes.size()   ;

    std::vector<bool> checked(vSize,0);

    std::vector<G4String     > volumeA, volumeB;
    std::vector<G4ThreeVector> intersections;

    G4cout << G4endl;

    int nInts;
    for(int i = 0; i < vSize; i++)
    {
        G4VSolid*   firstSolid = volumes[i]->GetLogicalVolume()->GetSolid();
        G4VisExtent firstVExt  = firstSolid->GetExtent();

        G4RotationMatrix firstRot = transform[i].getRotation   ();
        G4ThreeVector    firstPos = transform[i].getTranslation();

        nInts = 0;
        for(int event = 0; event < nEvents; event++)
        {
            if(event % (nEvents/10) == 0)
                G4cout << "Checking volume " << std::setw(15) << volumes[i]->GetName() << " ... " << std::setw(8) << event << "/" << nEvents << '\xd' << std::flush;

            Position = GetPointInVolume(firstVExt);
            while( firstSolid->Inside(Position) != kInside)
                Position = GetPointInVolume(firstVExt);			

            // Position in hall frame of reference
            Position *= firstRot;
            Position += firstPos;

            for(int j = 0; j < vSize; j++)
            {

                if(checked[j]) continue;
                if( j == i )   continue;

                G4VSolid*   secondSolid  = volumes[j]->GetLogicalVolume()->GetSolid();
		
                G4RotationMatrix secondRot = transform[j].getRotation   ();
                G4ThreeVector    secondPos = transform[j].getTranslation();

                secondRot = (G4RotationMatrix) secondRot.inverse();

                G4ThreeVector newPosition = Position;			

                newPosition -= secondPos;
                newPosition *= secondRot;

                if( secondSolid->Inside(newPosition) == kInside)
                {
                    volumeA.push_back(volumes[i]->GetName());
                    volumeB.push_back(volumes[j]->GetName());
                    intersections.push_back(Position);
                    nInts++;
                    checked[j] = 1;
                }
            }
        }
        std::fill(checked.begin(), checked.end(), 0);
        G4cout << "Checking volume " << std::setw(15) << volumes[i]->GetName() << " ..." << std::setw(9) << nEvents 
               << "/" << nEvents << " DONE ; FOUND " << nInts << " INTERSECTIONS" << G4endl;
    }
    G4cout << G4endl;

    if(volumeA.size() == 0)
        G4cout << "No intersections detected!" << G4endl;

    for(int i = 0; i < (int) volumeA.size(); i++)
    {
        G4cout << "Intersection at ( " <<
                  std::setw(10) << intersections[i][0] << " "  <<
                  std::setw(10) << intersections[i][1] << " "  <<
                  std::setw(10) << intersections[i][2] << " )" <<
                  " between volumes " <<
                  std::setw(15) << volumeA[i] << " and " <<
                  std::setw(15) << volumeB[i] << G4endl;
    }
    G4cout << G4endl;
}

G4ThreeVector qshieldsDebug::GetPointInVolume( const G4VisExtent VisExt )
{
    return ( G4ThreeVector( VisExt.GetXmin() + (VisExt.GetXmax()-VisExt.GetXmin())*G4UniformRand(),
                            VisExt.GetYmin() + (VisExt.GetYmax()-VisExt.GetYmin())*G4UniformRand(),
                            VisExt.GetZmin() + (VisExt.GetZmax()-VisExt.GetZmin())*G4UniformRand() ) );
}

int qshieldsDebug::CheckThisVolume(G4VPhysicalVolume* vol)
{
    G4String name = vol->GetName();

    for(int i = 0; i < 2; i++)
    {
        int val;
        if(i == 0) val = v1;
        else       val = v2;

        switch(val)
        {
			case(1):if(name=="PETLatX_P")
                return val; 
                break;
			case(2):if(name=="PETupLatX_P")
                return val; 
                break;
			case(3):if(name=="BAcidLatX_P")
                return val; 
                break;
			case(4):if(name=="ExtPbLatX_P")
                return val; 
                break;
			case(5):if(name=="CuRS1Lat_P")
                return val; 
                break;
			case(6):if(name=="CuRS1Flan_P")
                return val; 
                break;
			case(7):if(name=="CuRS2Lat_P")
                return val; 
                break;
			case(8):if(name=="CuRS2Flan_P")
                return val; 
                break;
			case(9):if(name=="CuRS3Lat_P")
                return val; 
                break;
			case(10):if(name=="CuRS3Flan_P")
                return val; 
                break;
			case(11):if(name=="CuRS4Lat_P")
                return val; 
                break;
			case(12):if(name=="CuRS4Flan_P")
                return val; 
                break;
			case(13):if(name=="PbRS4Lat_P" || name=="PbRS4Bottom_P")
                return val; 
                break;
			case(14):if(name=="PbRS4CuCrandle_P" || name=="PbRS4CuRing_P")
                return val; 
                break;
			case(15):if(name=="AcRods_P")
                return val; 
                break;
			case(16):if(name=="CuRS5Lat_P")
                return val; 
                break;
			case(17):if(name=="CuRS5Flan_P")
                return val; 
                break;
			case(18):if(name=="CuRS6Lat_P") 
                return val; 
                break;
			case(19):if(name=="CuRS6Flan_P")
                return val; 
                break;
			case(22):if(name=="CuRS6Tiles_P")
                return val; 
                break;
			case(20):if(name=="CuTopPlate_P")
                return val; 
                break;
			case(21):if(name=="StdPb_P")
                return val; 
                break;
			case(24):if(name=="CuBottomPlate_P")
                return val; 
                break;
			case(25):if(name=="CuTSP_P")
                return val; 
                break;
			case(26):if(name=="CuDetUpPlate_P")
                return val; 
                break;
			case(28):if(name=="CuDetBottomPlate_P")
                return val; 
                break;
			case(30):if(name=="IntCalTubeCu_P" || name=="StdPbDCS_P")
                return val; 
                break;
			case(31):if(name=="IntWeightCapsuleTeflonLayer_P" || name=="IntSourceCapsuleTeflonLayer_P")
                return val; 
                break;
			case(32):if(name=="IntWeightCapsuleCuLayer_P" || name=="IntSourceCapsuleCuLayer_P")
                return val; 
                break;
			case(34):if(name=="ExtWeightCapsuleTeflonLayer_P" || name=="ExtSourceCapsuleTeflonLayer_P" || name=="ExtWeightCapsuleCuLayer_P" || name=="ExtSourceCapsuleCuLayer_P")
                return val; 
                break;
			case(37):if(name=="TeO2_P")
                return val; 
                break;
			case(38):if(name=="FrameB_P" || name=="FrameC_P" || name=="FrameT_P" || name=="Column_P")
                return val; 
                break;
			case(39):if(name=="WireTray_P")
                return val; 
                break;
			case(40):if(name=="WirePad_P")
                return val; 
                break;
			case(41):if(name=="NTD_P")
                return val; 
                break;
			case(42):if(name=="PTFEBottom_P" || name=="PTFETop_P")
                return val; 
                break;
			case(43):if(name=="GiuntiCuBe_P")
                return val; 
                break;
			case(44):if(name=="PENTape_P")
                return val; 
                break;
			case(46):if(name=="SIRS2Lat_P" || name=="SIRS2Flan_P" || name=="SIRS3Lat_P" || name=="SIRS3Flan_P")
                return val; 
                break;
			case(50):if(name=="DummyLat_P")
                return val; 
                break;
			case(51):if(name=="DummyTop_P")
                return val; 
                break;
			case(54):if(name=="IntWeightTungsten_8_P")
                return val; 
                break;
			case(55):if(name=="IntWeightTungsten_21_P")
                return val; 
                break;
			case(56):if(name=="IntWeightTungsten_4_P")
                return val; 
                break;
			case(57):if(name=="IntWeightTungsten_1_P")
                return val; 
                break;
			case(58):if(name=="ExtWeightTungsten_8_P")
                return val; 
                break;
			case(59):if(name=="ExtSourceTungsten_20_P")
                return val; 
                break;
			case(60):if(name=="ExtSourceTungsten_5_P")
                return val; 
                break;
			default:
                break;
        }
    }
    return 0;
}
