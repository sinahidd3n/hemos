#include "hemocell.h"

#include "hemoCellParticleDataTransfer.h"
#include "rbcHighOrderModel.h"
#include "hemoCellField.h"
#include "immersedBoundaryMethod.h"
#include "config.h"
#include "logfile.h"
#include "hemoCellParticleField.h"
#include "readPositionsBloodCells.h"
#include "meshMetrics.h"
#include "meshGeneratingFunctions.h"
#include "pltSimpleModel.h"
#include "WBCModel.h"
#include "CTC.h"
#include "cellInfo.h"
#include "fluidInfo.h"
#include "particleInfo.h"
#include "writeCellInfoCSV.h"
#include <fenv.h>
#include <stdio.h> 
#include <stdlib.h> 
#include <time.h> 


int main(int argc, char *argv[]) {
  if(argc < 2) {
    cout << "Usage: " << argv[0] << " <configuration.xml>" << endl;
    return -1;
  }

  HemoCell hemocell(argv[1], argc, argv);
  Config * cfg = hemocell.cfg;
  
  time_t start_time = time(NULL);
  printf("%s", ctime(&start_time));

  hlogfile << "(PipeFlow) (Geometry) reading and voxelizing STL file " << (*cfg)["domain"]["geometry"].read<string>() << endl; 

  std::auto_ptr<MultiScalarField3D<int>> flagMatrix;
  std::auto_ptr<VoxelizedDomain3D<T>> voxelizedDomain; 
  getFlagMatrixFromSTL((*cfg)["domain"]["geometry"].read<string>(),  
                       (*cfg)["domain"]["fluidEnvelope"].read<int>(),  
                       (*cfg)["domain"]["refDirN"].read<int>(),  
                       (*cfg)["domain"]["refDir"].read<int>(),  
                       voxelizedDomain, flagMatrix,  
                       (*cfg)["domain"]["blockSize"].read<int>(),
                       (*cfg)["domain"]["particleEnvelope"].read<int>()); 

  param::lbm_pipe_parameters((*cfg),flagMatrix.get());
  param::printParameters();
  
  hemocell.lattice = new MultiBlockLattice3D<T, DESCRIPTOR>(
            voxelizedDomain.get()->getMultiBlockManagement(),
            defaultMultiBlockPolicy3D().getBlockCommunicator(),
            defaultMultiBlockPolicy3D().getCombinedStatistics(),
            defaultMultiBlockPolicy3D().getMultiCellAccess<T, DESCRIPTOR>(),
            new GuoExternalForceBGKdynamics<T, DESCRIPTOR>(1.0/param::tau));

  defineDynamics(*hemocell.lattice, *flagMatrix.get(), (*hemocell.lattice).getBoundingBox(), new BounceBack<T, DESCRIPTOR>(1.), 0);

  hemocell.lattice->toggleInternalStatistics(false);
  hemocell.lattice->periodicity().toggleAll(false);
  hemocell.latticeEquilibrium(1.,plb::Array<T, 3>(0.,0.,0.));

  //Driving Force
  T poiseuilleForce =  8 * param::nu_lbm * (param::u_lbm_max * 0.5) / param::pipe_radius / param::pipe_radius;
  setExternalVector(*hemocell.lattice, (*hemocell.lattice).getBoundingBox(),
                    DESCRIPTOR<T>::ExternalField::forceBeginsAt,
                    plb::Array<T, DESCRIPTOR<T>::d>(poiseuilleForce, 0.0, 0.0));

  hemocell.lattice->initialize();   

  //Adding all the cells
  hemocell.initializeCellfield();

  hemocell.addCellType<RbcHighOrderModel>("RBC_HO", RBC_FROM_SPHERE);
  hemocell.setMaterialTimeScaleSeparation("RBC_HO", (*cfg)["ibm"]["stepMaterialEvery"].read<int>());
  hemocell.setMinimumDistanceFromSolid("RBC_HO", 0.5); //Micrometer! not LU

  hemocell.addCellType<WBCModel>("WBC", MESH_FROM_STL);
  hemocell.setMaterialTimeScaleSeparation("WBC", (*cfg)["ibm"]["stepMaterialEvery"].read<int>());

  hemocell.addCellType<CTC>("CTC", MESH_FROM_STL);
  hemocell.setMaterialTimeScaleSeparation("CTC", (*cfg)["ibm"]["stepMaterialEvery"].read<int>());

  hemocell.addCellType<PltSimpleModel>("PLT", MESH_FROM_STL); // modified mesh has been added to Hemocell which gives us more accuracy but comes with higher computational cost (needed for investigating the effect of adhesion)
  hemocell.setMaterialTimeScaleSeparation("PLT", (*cfg)["ibm"]["stepMaterialEvery"].read<int>());
  
  
  hemocell.setParticleVelocityUpdateTimeScaleSeparation((*cfg)["ibm"]["stepParticleEvery"].read<int>());

  hemocell.setRepulsion((*cfg)["domain"]["kRep"].read<T>(), (*cfg)["domain"]["RepCutoff"].read<T>());
  hemocell.setRepulsionTimeScaleSeperation((*cfg)["ibm"]["stepMaterialEvery"].read<int>());

  //hemocell.setRepulsion((*cfg)["domain"]["kRep"].read<T>(), (*cfg)["domain"]["RepCutoff"].read<T>());
//  hemocell.setRepulsionTimeScaleSeperation((*cfg)["ibm"]["stepMaterialEvery"].read<int>());

  vector<int> outputs = {OUTPUT_POSITION,OUTPUT_TRIANGLES,OUTPUT_FORCE,OUTPUT_FORCE_VOLUME,OUTPUT_FORCE_BENDING,OUTPUT_FORCE_LINK,OUTPUT_FORCE_AREA,OUTPUT_FORCE_VISC,OUTPUT_FORCE_REPULSION,OUTPUT_CELL_ID, OUTPUT_VERTEX_ID};
  hemocell.setOutputs("RBC_HO", outputs);
  hemocell.setOutputs("WBC", outputs);
  hemocell.setOutputs("CTC", outputs);
  hemocell.setOutputs("PLT", outputs);

  outputs = {OUTPUT_VELOCITY,OUTPUT_DENSITY,OUTPUT_FORCE};
  hemocell.setFluidOutputs(outputs);

  // Turn on periodicity in the X direction
  hemocell.setSystemPeriodicity(0, true);

  // Enable boundary particles
  hemocell.enableBoundaryParticles((*cfg)["domain"]["BkRep"].read<T>(), (*cfg)["domain"]["BRepCutoff"].read<T>(),(*cfg)["ibm"]["stepMaterialEvery"].read<int>());
  
  //loading the cellfield
  if (not cfg->checkpointed) {
    hemocell.loadParticles();
    hemocell.writeOutput();
  } else {
    hemocell.loadCheckPoint();
  }

  //Restructure atomic blocks on processors when possible
  //hemocell.doRestructure(false); // cause errors
  
  if (hemocell.iter == 0) {
    hlog << "(PipeFlow) fresh start: warming up cell-free fluid domain for "  << (*cfg)["parameters"]["warmup"].read<plint>() << " iterations..." << endl;
    for (plint itrt = 0; itrt < (*cfg)["parameters"]["warmup"].read<plint>(); ++itrt) { 
      hemocell.lattice->collideAndStream(); 
    }
  }

  unsigned int tmax = (*cfg)["sim"]["tmax"].read<unsigned int>();
  unsigned int tmeas = (*cfg)["sim"]["tmeas"].read<unsigned int>();
  unsigned int tadh = (*cfg)["sim"]["tadh"].read<unsigned int>();
  unsigned int tcheckpoint = (*cfg)["sim"]["tcheckpoint"].read<unsigned int>();
  unsigned int tcsv = (*cfg)["sim"]["tcsv"].read<unsigned int>();

  hlog << "(PipeFlow) Starting simulation..." << endl;


 


  while (hemocell.iter < tmax ) {
    hemocell.iterate();
    

    
    //Set driving force as required after each iteration
    setExternalVector(*hemocell.lattice, hemocell.lattice->getBoundingBox(),
                DESCRIPTOR<T>::ExternalField::forceBeginsAt,
                plb::Array<T, DESCRIPTOR<T>::d>(poiseuilleForce, 0.0, 0.0));

    

    //if (hemocell.iter % tmeas == 0) {
    if (hemocell.iter % tcsv == 0){

   
        hlog << "(main) Stats. @ " <<  hemocell.iter << " (" << hemocell.iter * param::dt << " s):" << endl;
        hlog << "\t # of cells: " << CellInformationFunctionals::getTotalNumberOfCells(&hemocell);
        hlog << " | # of RBC: " << CellInformationFunctionals::getNumberOfCellsFromType(&hemocell, "RBC_HO");
        hlog << " | # of WBC: " << CellInformationFunctionals::getNumberOfCellsFromType(&hemocell, "WBC");
        hlog << " | # of CTC: " << CellInformationFunctionals::getNumberOfCellsFromType(&hemocell, "CTC");
        hlog << ", PLT: " << CellInformationFunctionals::getNumberOfCellsFromType(&hemocell, "PLT") << endl;
        FluidStatistics finfo = FluidInfo::calculateVelocityStatistics(&hemocell); T toMpS = param::dx / param::dt;
        hlog << "\t Velocity  -  max.: " << finfo.max * toMpS << " m/s, mean: " << finfo.avg * toMpS<< " m/s, rel. app. viscosity: " << (param::u_lbm_max*0.5) / finfo.avg << endl;
        ParticleStatistics pinfo = ParticleInfo::calculateForceStatistics(&hemocell); T topN = param::df * 1.0e12;
        hlog << "\t Force  -  min.: " << pinfo.min * topN << " pN, max.: " << pinfo.max * topN << " pN (" << pinfo.max << " lf), mean: " << pinfo.avg * topN << " pN" << endl;


        //previous modifications (Modified by Sina [sanvarin@uwaterloo.ca])
       /* hlog << endl<< endl<< endl<< endl<< endl<< endl<< endl<< endl<< endl<< endl<< endl<< endl;
        
        for (const HemoCellParticle & particle : particles) {
        	hlog<<((*cellFields)[particle.sv.celltype]->doInteriorViscosity)


        hlog << endl<< endl<< endl<< endl<< endl<< endl<< endl<< endl<< endl<< endl<< endl<< endl;

        */ 




        // Additional useful stats, if needed
        //finfo = FluidInfo::calculateForceStatistics(&hemocell);
        //Set force as required after this function;
        // setExternalVector(*hemocell.lattice, hemocell.lattice->getBoundingBox(),
        //           DESCRIPTOR<T>::ExternalField::forceBeginsAt,
        //           hemo::Array<T, DESCRIPTOR<T>::d>(poiseuilleForce, 0.0, 0.0));
        // pcout << "Fluid force, Minimum: " << finfo.min << " Maximum: " << finfo.max << " Average: " << finfo.avg << endl;
        // ParticleStatistics pinfo = ParticleInfo::calculateVelocityStatistics(&hemocell);
        // pcout << "Particle velocity, Minimum: " << pinfo.min << " Maximum: " << pinfo.max << " Average: " << pinfo.avg << endl;
        hemocell.writeOutput();

        
      // writeCellInfo_CSV(hemocell);

        hlog << "Saving simple mean cell values to CSV at timestep " << hemocell.iter << endl;

        ///////////////////Sina//////////////////
      
        //previous modifications by sina [sanvarin@uwaterloo.ca] 

        //hlog << endl<< "----------------------------------------------------------------" << endl<< endl ;

       //previous modifications (Modified by Sina [sanvarin@uwaterloo.ca])

    /*   map<int,CellInformation> info_per_cell;
        CellInformationFunctionals::calculateCellInformation(&hemocell,info_per_cell);
        HemoCellGatheringFunctional<CellInformation>::gather(info_per_cell);
        for (auto & pair : info_per_cell) {
           const plint cid = pair.first;
           CellInformation & cinfo = pair.second;
        hlog<<endl<<cinfo.position[1]<<endl<<endl;
      }
      */
        




        /*map<int,CellInformation> info_per_cell;
        CellInformationFunctionals::calculateCellInformation(&hemocell,info_per_cell);
                for (auto it = info_per_cell.cbegin(); it != info_per_cell.cend() ;) 
  {
    if (!it->second.centerLocal) {
      it = info_per_cell.erase(it);
    } else {
      ++it;
    }

  }


        CellInformationFunctionals::calculateCellPosition(&hemocell);
        CellInformationFunctionals::calculateCellVolume(&hemocell);


        for (auto & pair : info_per_cell) {
      const plint cid = pair.first;
      CellInformation & cinfo = pair.second;
      hlog <<endl <<cinfo.position[0] << ", " << cinfo.position[1] << ", " << cinfo.position[2] <<", "<<cid<<endl<<endl;
    }
        T volume3 = (CellInformationFunctionals::info_per_cell[2].volume);
         hlog<<volume3<<"  "<<endl; 
         */


      
    






         /*T volume = (CellInformationFunctionals::info_per_cell[0].volume);
         T volume2 = (CellInformationFunctionals::info_per_cell[1].volume);
         T volume3 = (CellInformationFunctionals::info_per_cell[2].volume);
         hlog<<volume<<"   "<<volume2<<"   "<<volume3<<"  "<<endl; */





  
         /*hemo::Array<T,3> position1 = CellInformationFunctionals::info_per_cell[0].position;
         hlog << "\t Cell center at: {" <<position1[0]<<","<<position1[1]<<","<<position1[2] << "} lu" <<endl; 

         hemo::Array<T,3> position2 = CellInformationFunctionals::info_per_cell[1].position;
         hlog << "\t Cell center at: {" <<position2[0]<<","<<position2[1]<<","<<position2[2] << "} lu" <<endl;

         hemo::Array<T,3> position3 = CellInformationFunctionals::info_per_cell[3].position;
         hlog << "\t Cell center at: {" <<position3[0]<<","<<position3[1]<<","<<position3[2] << "} lu" <<endl;


      //float distt01;
      //float distt02;
      
      distt01= sqrt(pow((poss[1]-poss[5]),2)+pow((poss[2]-poss[6]),2)+pow((poss[3]-poss[7]),2));
      distt02= sqrt(pow((poss[1]-poss[9]),2)+pow((poss[2]-poss[10]),2)+pow((poss[3]-poss[11]),2));

      hlog<< endl<< endl<< endl<< endl<<distt02<< endl  <<distt01<< endl<< endl<< endl;
*/
         
       
        

       

      //hlog << endl<< "----------------------------------------------------------------" << endl<< endl ;

      ///////////////////////////Sina////////////////


        
    }
    
    if (hemocell.iter % tcsv == 0) {
      


      

    }
    
    if (hemocell.iter % tcheckpoint == 0) {
      hemocell.saveCheckPoint();
    }
  }

  hlog << "(main) Simulation finished :) " << endl;
  time_t end_time = time(NULL);
  printf("%s", ctime(&end_time));

  return 0;
}
