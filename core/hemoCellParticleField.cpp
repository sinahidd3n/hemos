/*
This file is part of the HemoCell library

HemoCell is developed and maintained by the Computational Science Lab 
in the University of Amsterdam. Any questions or remarks regarding this library 
can be sent to: info@hemocell.eu

When using the HemoCell library in scientific work please cite the
corresponding paper: https://doi.org/10.3389/fphys.2017.00563

The HemoCell library is free software: you can redistribute it and/or
modify it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

The library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "hemoCellParticleField.h"
#include "hemocell.h"
#include "octree.h"
#include "mollerTrumbore.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <cstdlib>
#include <ctime>
#include "cellInfo.h"


namespace hemo { 
/* *************** class HemoParticleField3D ********************** */

HemoCellParticleField::HemoCellParticleField(plint nx, plint ny, plint nz)
    : AtomicBlock3D(nx,ny,nz,0), particleDataTransfer(*(new HemoCellParticleDataTransfer()))
{ 
    boundingBox = Box3D(0,this->getNx()-1, 0, this->getNy()-1, 0, this->getNz()-1);
    dataTransfer = &particleDataTransfer;
    particleDataTransfer.setBlock(*this);
    AddOutputMap(); 
}

HemoCellParticleField::HemoCellParticleField(HemoCellParticleField const& rhs)
    : AtomicBlock3D(rhs), particleDataTransfer(*(new HemoCellParticleDataTransfer()))
{
    boundingBox = Box3D(0,this->getNx()-1, 0, this->getNy()-1, 0, this->getNz()-1);
    dataTransfer = &particleDataTransfer;
    particleDataTransfer.setBlock(*this);
    for (const HemoCellParticle & particle : rhs.particles) {
      addParticle(particle.sv);
    }
    ppc_up_to_date = false;
    lpc_up_to_date = false;
    ppt_up_to_date = false;
    pg_up_to_date = false;
    AddOutputMap();
}

HemoCellParticleField::~HemoCellParticleField()
{
  //AtomicBlock3D::dataTransfer = new HemoCellParticleDataTransfer();
  if (particle_grid) {
    delete[] particle_grid;
    particle_grid = 0;
  }
  if(particle_grid_size) {
    delete[] particle_grid_size;
    particle_grid_size = 0;
  }
}

HemoCellParticleField& HemoCellParticleField::operator=(HemoCellParticleField const& rhs){
 HemoCellParticleField *copy = new HemoCellParticleField(rhs);
 delete this;
  return *copy;
}

HemoCellParticleField* HemoCellParticleField::clone() const
{
    return new HemoCellParticleField(*this);
}

const vector<vector<unsigned int>> & HemoCellParticleField::get_particles_per_type() { 
    if (!ppt_up_to_date) { update_ppt(); }
    return _particles_per_type;
  }
const map<int,vector<int>> & HemoCellParticleField::get_particles_per_cell() { 
    if (!ppc_up_to_date) { update_ppc(); }
    return _particles_per_cell;
  }

const map<int,bool> & HemoCellParticleField::get_lpc() { 
    if (!lpc_up_to_date) { update_lpc(); }
    return _lpc;
  }
void HemoCellParticleField::update_lpc() {
  _lpc.clear();
  for (const HemoCellParticle & particle : particles) {
     if (isContainedABS(particle.sv.position, localDomain)) {
       _lpc[particle.sv.cellId] = true;
     }
  }
  lpc_up_to_date = true;
}
void HemoCellParticleField::update_ppt() {
  _particles_per_type.clear();
  _particles_per_type.resize(cellFields->size());
  
  for (unsigned int i = 0 ; i <  particles.size() ; i++) { 
    _particles_per_type[particles[i].sv.celltype].push_back(i);
  }
  ppt_up_to_date = true;
}
void HemoCellParticleField::update_ppc() {
  _particles_per_cell.clear();
  
  for (unsigned int i = 0 ; i <  particles.size() ; i++) { 
     insert_ppc(&particles[i],i);
  }
  ppc_up_to_date = true;
}

void HemoCellParticleField::update_pg() {
  //Check if map exists, otherwise create
  if (!particle_grid) {
    particle_grid = new hemo::Array<unsigned int, 10>[this->atomicLattice->getNx()*this->atomicLattice->getNy()*this->atomicLattice->getNz()];
  }
  if (!particle_grid_size) {
    particle_grid_size = new unsigned int[this->atomicLattice->getNx()*this->atomicLattice->getNy()*this->atomicLattice->getNz()];
  }
  if (!this->atomicLattice) {
    return;
  }
  
  memset(particle_grid_size,0,sizeof(unsigned int)*this->atomicLattice->getNx()*this->atomicLattice->getNy()*this->atomicLattice->getNz());
  Dot3D const& location = this->atomicLattice->getLocation();
  hemo::Array<T,3> * pos;
  
  for (unsigned int i = 0 ; i <  particles.size() ; i++) {
    pos = &particles[i].sv.position;
    int x = pos->operator[](0)-location.x+0.5;
    int y = pos->operator[](1)-location.y+0.5;
    int z = pos->operator[](2)-location.z+0.5;
    if ((x >= 0) && (x < this->atomicLattice->getNx()) &&
	(y >= 0) && (y < this->atomicLattice->getNy()) &&
	(z >= 0) && (z < this->atomicLattice->getNz()) ) 
    {
      unsigned int index = grid_index(x,y,z);
      particle_grid[index][particle_grid_size[index]] = i;
      particle_grid_size[index]++;
    }
  }
  pg_up_to_date = true;
}

void HemoCellParticleField::addParticle(HemoCellParticle* particle) {
  addParticle(particle->sv);
}  
void HemoCellParticleField::addParticle(const HemoCellParticle::serializeValues_t & sv) {
  HemoCellParticle * local_sparticle, * particle;
  const hemo::Array<T,3> & pos = sv.position;
  const map<int,vector<int>> & particles_per_cell = get_particles_per_cell();

  if( this->isContainedABS(pos, this->getBoundingBox()) )
  {
    //check if we have particle already, if so, we must overwrite but not
    //forget to delete the old entry
    if ((!(particles_per_cell.find(sv.cellId) == particles_per_cell.end()))) { 
      if (particles_per_cell.at(sv.cellId)[sv.vertexId] != -1) {
        local_sparticle =  &particles[particles_per_cell.at(sv.cellId)[sv.vertexId]];

        //If our particle is local, do not replace it, envelopes are less important
        if (isContainedABS(local_sparticle->sv.position, localDomain)) {
          return;
        } else {
          //We have the particle already, replace it
          local_sparticle->sv = sv;
          particle = local_sparticle;
          particle->setTag(-1);

          //Invalidate lpc hemo::Array
          lpc_up_to_date = false;
          pg_up_to_date = false;

        }
      } else {
        goto outer_else;
      }
    } else {
outer_else:
      //new entry
      particles.emplace_back(sv);
      particle = &particles.back();
      
      //invalidate ppt
      ppt_up_to_date=false;
        if(this->isContainedABS(pos, localDomain)) {
          _lpc[particle->sv.cellId] = true;
        }
        if (ppc_up_to_date) { //Otherwise its rebuild anyway
         insert_ppc(particle, particles.size()-1);
        }
      
      if (pg_up_to_date) {
        Dot3D const& location = this->atomicLattice->getLocation();
        hemo::Array<T,3>  & pos = particle->sv.position;
        int x = pos[0]-location.x+0.5;
        int y = pos[1]-location.y+0.5;
        int z = pos[2]-location.z+0.5;
        if ((x >= 0) && (x <= this->atomicLattice->getNx()) &&
            (y >= 0) && (y <= this->atomicLattice->getNy()) &&
            (z >= 0) && (z <= this->atomicLattice->getNz()) ) 
        {
          unsigned int index = grid_index(x,y,z);
          particle_grid[index][particle_grid_size[index]] = particles.size()-1;
          particle_grid_size[index]++;
        }
      }
    }
  }
}

void inline HemoCellParticleField::insert_ppc(HemoCellParticle* sparticle, unsigned int index) {
  if (_particles_per_cell.find(sparticle->sv.cellId) == _particles_per_cell.end()) {
    _particles_per_cell[sparticle->sv.cellId].resize((*cellFields)[sparticle->sv.celltype]->numVertex,-1);
  }
  _particles_per_cell.at(sparticle->sv.cellId)[sparticle->sv.vertexId] = index;

}
void inline HemoCellParticleField::insert_preinlet_ppc(HemoCellParticle* sparticle, unsigned int index) {
  if (_preinlet_particles_per_cell.find(sparticle->sv.cellId) == _preinlet_particles_per_cell.end()) {
    _preinlet_particles_per_cell[sparticle->sv.cellId].resize((*cellFields)[sparticle->sv.celltype]->numVertex);
    for (unsigned int i = 0; i < _preinlet_particles_per_cell[sparticle->sv.cellId].size(); i++) {
      _preinlet_particles_per_cell[sparticle->sv.cellId][i] = -1;
    }
  }
  _preinlet_particles_per_cell.at(sparticle->sv.cellId)[sparticle->sv.vertexId] = index;

}

void HemoCellParticleField::removeParticles(plint tag) {
//Almost the same, but we save a lot of branching by making a seperate function

  const unsigned int old_size = particles.size();
  for (unsigned int i = 0 ; i < particles.size() ; i++) {
    if (particles[i].getTag() == tag) {
      particles[i] = particles.back();
      particles.pop_back();
  //   if (particles[i].sv.cellId==-1 || particles[i].sv.cellId==1) {cout<<endl<<"Warniiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiing removeParticles1 "<<endl;}
      i--;
    }
  }
  if (particles.size() != old_size) {
    lpc_up_to_date = false;
    ppt_up_to_date = false;
    ppc_up_to_date = false;
    pg_up_to_date = false;
  } 
}

void HemoCellParticleField::removeParticles(Box3D domain, plint tag) {
//Almost the same, but we save a lot of branching by making a seperate function
  Box3D finalDomain;
  
  intersect(domain, this->getBoundingBox(), finalDomain);

  const unsigned int old_size = particles.size();
  for (unsigned int i = 0 ; i < particles.size() ; i++) {
    if (particles[i].getTag() == tag && this->isContainedABS(particles[i].sv.position,finalDomain)) {
      particles[i] = particles.back();
      particles.pop_back();
   //   if (particles[i].sv.cellId==-1 || particles[i].sv.cellId==1) {cout<<endl<<"Warniiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiing removeParticles2 "<<endl;}
      i--;
    }
  }
  if (particles.size() != old_size) {
    lpc_up_to_date = false;
    ppt_up_to_date = false;
    ppc_up_to_date = false;
    pg_up_to_date = false;
  } 
}

void HemoCellParticleField::removeParticles(Box3D domain) {
//Almost the same, but we save a lot of branching by making a seperate function

  Box3D finalDomain;
  
  intersect(domain, this->getBoundingBox(), finalDomain);

  const unsigned int old_size = particles.size();
  for (unsigned int i = 0 ; i < particles.size() ; i++) {
    if (this->isContainedABS(particles[i].sv.position,finalDomain)) {
      particles[i] = particles.back();
      particles.pop_back();
  //   if (particles[i].sv.cellId==-1 || particles[i].sv.cellId==1) {cout<<endl<<"Warniiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiing removeParticles3 "<<endl;}
      i--;
    }
  }
  if (particles.size() != old_size) {
    lpc_up_to_date = false;
    ppt_up_to_date = false;
    ppc_up_to_date = false;
    pg_up_to_date = false;
  } 
}

//remove everything outside this domain
void HemoCellParticleField::removeParticles_inverse(Box3D domain) {
//Almost the same, but we save a lot of branching by making a seperate function

  Box3D finalDomain;
  
  intersect(domain, this->getBoundingBox(), finalDomain);

  const unsigned int old_size = particles.size();
  for (unsigned int i = 0 ; i < particles.size() ; i++) {
    if (!this->isContainedABS(particles[i].sv.position,finalDomain)) {
      particles[i] = particles.back();
      particles.pop_back();
    // if (particles[i].sv.cellId==-1 || particles[i].sv.cellId==1) {cout<<endl<<"Warniiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiing removeParticles4 "<<endl;}
      i--;
    }
  }
  if (particles.size() != old_size) {
    lpc_up_to_date = false;
    ppt_up_to_date = false;
    ppc_up_to_date = false;
    pg_up_to_date = false;
  } 
}

void HemoCellParticleField::syncEnvelopes() {
  removeParticles_inverse(localDomain);
}

void HemoCellParticleField::findParticles (
        Box3D domain, std::vector<HemoCellParticle*>& found )
{
    found.clear();
    PLB_ASSERT( contained(domain, this->getBoundingBox()) );
    for (HemoCellParticle & particle : particles) {
        if (this->isContainedABS(particle.sv.position,domain)) {
            found.push_back(&particle);
        }
    }
}

void HemoCellParticleField::findParticles (
        Box3D domain, std::vector<const HemoCellParticle*>& found ) const
{
    found.clear();
    PLB_ASSERT( contained(domain, this->getBoundingBox()) );
    for (const HemoCellParticle & particle : particles) {
        if (this->isContainedABS(particle.sv.position,domain)) {
            found.push_back(&particle);
        }
    }
}
void HemoCellParticleField::findParticles (
        Box3D domain, std::vector<HemoCellParticle*>& found, pluint type)
{
    
    found.clear();
    PLB_ASSERT( contained(domain, this->getBoundingBox()) );
    //hemo::Array<T,3> pos; 
    const vector<vector<unsigned int>> & particles_per_type = get_particles_per_type();
    if (!(particles_per_type.size() > type)) 
      {return;} 
    else {
      for (const unsigned int i : particles_per_type[type]) {
          if (this->isContainedABS(particles[i].sv.position,domain)) {
              found.push_back(&(particles[i]));
          }
      }
    }
    
}

inline plint HemoCellParticleField::nearestCell(T const pos) const {
  return int(pos + 0.5);
}

inline void HemoCellParticleField::computeGridPosition (
            hemo::Array<T,3> const& position,
                    plint* iX, plint* iY, plint* iZ ) const
{
      Dot3D const& location = this->getLocation();
      *iX = nearestCell(position[0]) - location.x;
      *iY = nearestCell(position[1]) - location.y;
      *iZ = nearestCell(position[2]) - location.z;
}

void HemoCellParticleField::computeGridPosition (
            hemo::Array<T,3> const& position,
                    plint& iX, plint& iY, plint& iZ ) const
{
      Dot3D const& location = this->getLocation();
      iX = nearestCell(position[0]) - location.x;
      iY = nearestCell(position[1]) - location.y;
      iZ = nearestCell(position[2]) - location.z;
}

void HemoCellParticleField::issueWarning(HemoCellParticle & p){
	cout << "(HemoCell) (Delete Cells) WARNING! Particle deleted from local domain. This means the whole cell will be deleted!" << endl;
  cout << "\t Particle ID:" << p.sv.cellId << endl;
  cout << "\t Position: " << p.sv.position[0] << ", " << p.sv.position[1] << ", " << p.sv.position[2] << "; vel.: " << p.sv.v[0] << ", " <<  p.sv.v[1] << ", " << p.sv.v[2] << "; force: " << p.sv.force[0] << ", " << p.sv.force[1] << ", " << p.sv.force[2] << endl;
}

int HemoCellParticleField::deleteIncompleteCells(pluint ctype, bool verbose) {
  int deleted = 0;

  const map<int,vector<int>> & particles_per_cell = get_particles_per_cell();
  //Warning, TODO, high complexity, should be rewritten 
  //For now abuse tagging and the remove function
  for ( const auto &lpc_it : particles_per_cell ) {
    int cellid = lpc_it.first;
    bool broken = false;
    for (pluint i = 0; i < particles_per_cell.at(cellid).size() ; i++) {
      if (particles_per_cell.at(cellid)[i] == -1) {
        broken = true;
        //if (particles[i].sv.cellId==-1 || particles[i].sv.cellId==1) { cout<<endl<<"Warniiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiing delete incomplete cell 1 "<<endl;}
        break;
      }
    }
    if (!broken) {continue;}

    bool warningIssued = false;
    for (pluint i = 0; i < particles_per_cell.at(cellid).size() ; i++) {
      if (particles_per_cell.at(cellid)[i] == -1) {continue;}

      //issue warning
      if (verbose) {
        if (!warningIssued) {
          if (isContainedABS(particles[particles_per_cell.at(cellid)[i]].sv.position,localDomain)) {
                  issueWarning(particles[particles_per_cell.at(cellid)[i]]);
            warningIssued = true;
          }
        }
      }
      
      //actually add to tobedeleted list
      particles[particles_per_cell.at(cellid)[i]].setTag(1);
      deleted++;
    }
  } 

  //We have our list, now abuse the removeall function
  removeParticles(1);

  return deleted; 
}

int HemoCellParticleField::deleteIncompleteCells(const bool verbose) {
  int deleted = 0;
  const map<int,vector<int>> & particles_per_cell = get_particles_per_cell();
  //Warning, TODO, high complexity, should be rewritten 
  //For now abuse tagging and the remove function
  for ( const auto &lpc_it : particles_per_cell ) {
    int cellid = lpc_it.first;
    bool broken = false;
    for (pluint i = 0; i < particles_per_cell.at(cellid).size() ; i++) {
      if (particles_per_cell.at(cellid)[i] == -1) {
        broken = true;
       //if (particles[i].sv.cellId==-1 || particles[i].sv.cellId==1) {cout<<endl<<"Warniiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiing delete incomplete cell 2 "<<endl;}
        break;
      }
    }
    if (!broken) {continue;}

    bool warningIssued = false;
    for (pluint i = 0; i < particles_per_cell.at(cellid).size() ; i++) {
      if (particles_per_cell.at(cellid)[i] == -1) {continue;}

      //issue warning
      if (verbose) {
        if (!warningIssued) {
          if (isContainedABS(particles[particles_per_cell.at(cellid)[i]].sv.position,localDomain)) {
                  issueWarning(particles[particles_per_cell.at(cellid)[i]]);
            warningIssued = true;
          }
        }
      }
      
      //actually add to tobedeleted list
      particles[particles_per_cell.at(cellid)[i]].setTag(1);
      deleted++;
    }
  } 

  //We have our list, now abuse the removeall function

  removeParticles(1);

  return deleted; 
}

void HemoCellParticleField::setlocalDomain(Box3D & localDomain_) {
  localDomain = localDomain_;
  localDomain.x0 -= this->getLocation().x;
  localDomain.x1 -= this->getLocation().x;
  localDomain.y0 -= this->getLocation().y;
  localDomain.y1 -= this->getLocation().y;
  localDomain.z0 -= this->getLocation().z;
  localDomain.z1 -= this->getLocation().z;
}


void HemoCellParticleField::advanceParticles() {
  for(HemoCellParticle & particle:particles){
    particle.advance();
    //By lack of better place, check if it is on a boundary, if so, delete it
    plb::Box3D const box = atomicLattice->getBoundingBox();
    plb::Dot3D const& location = atomicLattice->getLocation();
    plint x = (particle.sv.position[0]-location.x)+0.5;
    plint y = (particle.sv.position[1]-location.y)+0.5;
    plint z = (particle.sv.position[2]-location.z)+0.5;

    if ((x >= box.x0) && (x <= box.x1) &&
	(y >= box.y0) && (y <= box.y1) &&
	(z >= box.z0) && (z <= box.z1)) {
      if (atomicLattice->get(x,y,z).getDynamics().isBoundary()) {
        particle.tag = 1;
      }
    }
  }
  removeParticles(1);
  
  lpc_up_to_date = false;
  pg_up_to_date = false;
}

void HemoCellParticleField::separateForceVectors() {
  //Also save the total force, therfore recalculate in advance
  applyConstitutiveModel();

  for (HemoCellParticle & sparticle : particles) {
    //Save Total Force
    sparticle.force_total = sparticle.sv.force + sparticle.sv.force_repulsion;

    //Just repoint all possible outputs for now //TODO only repoint the ones we
    //want

    sparticle.force_volume = new hemo::Array<T,3>({0.0,0.0,0.0});
    allocated_for_output.push_back(sparticle.force_volume);
    sparticle.force_link = new hemo::Array<T,3>({0.0,0.0,0.0});
    allocated_for_output.push_back(sparticle.force_link);
    sparticle.force_area = new hemo::Array<T,3>({0.0,0.0,0.0});
    allocated_for_output.push_back(sparticle.force_area);
    sparticle.force_bending = new hemo::Array<T,3>({0.0,0.0,0.0});
    allocated_for_output.push_back(sparticle.force_bending);
    sparticle.force_visc = new hemo::Array<T,3>({0.0,0.0,0.0});
    allocated_for_output.push_back(sparticle.force_visc);
    sparticle.force_inner_link = new hemo::Array<T,3>({0.0,0.0,0.0});
    allocated_for_output.push_back(sparticle.force_inner_link);
  }
}

void HemoCellParticleField::unifyForceVectors() {
  for (const hemo::Array<T,3>* mem : allocated_for_output) {
    delete mem;
  }
  allocated_for_output.clear();
  for (HemoCellParticle& sparticle : particles) {
    sparticle.repoint_force_vectors();
  }
}

void HemoCellParticleField::applyConstitutiveModel(bool forced) {
  map<int,vector<HemoCellParticle*>> * ppc_new = new map<int,vector<HemoCellParticle*>>();
  const map<int,vector<int>> & particles_per_cell = get_particles_per_cell();
  map<int,bool> lpc;
  //Fill it here, probably needs optimization, ah well ...
  for (const auto & pair : particles_per_cell) {
    const int & cid = pair.first;
    const vector<int> & cell = pair.second; 
    (*ppc_new)[cid].resize(cell.size());
    for (unsigned int i = 0 ; i < cell.size() ; i++) {
      if (cell[i] == -1) {
        (*ppc_new).erase(cid); //not complete, remove entry
        goto no_add_lpc;
      } else {
        (*ppc_new)[cid][i] = &particles[cell[i]];
      }
    }
    lpc[cid]=true;
    no_add_lpc:;
  }
  
  for (pluint ctype = 0; ctype < (*cellFields).size(); ctype++) {
    if ((*cellFields).hemocell.iter % (*cellFields)[ctype]->timescale == 0 || forced) {
      vector<HemoCellParticle*> found;
      findParticles(getBoundingBox(),found,ctype);
      if (found.size() > 0) {
        //only reset forces when the forces actually point at it.
        if (found[0]->force_area == &found[0]->sv.force) {
          for (HemoCellParticle* particle : found) {
            particle->sv.force = {0.,0.,0.};
#ifdef INTERIOR_VISCOSITY
            particle->normalDirection = {0., 0., 0.};
#endif
          }
        }
      }
      (*cellFields)[ctype]->mechanics->ParticleMechanics(*ppc_new,lpc,ctype);
    }
  }
  
  delete ppc_new;
  
}

// modified by Sina Anvari Naeini (sanvarin@uwaterloo.ca)
/*
            bool stochastic_model=0;
            if (stochastic_model){
                srand(time(0));
                int p=rand()%10000;
      	        float k_on0=10000;
          	    float k_off0=20;
          	    float sigma_on=0.0000005;
          	    float sigma_off=0.00000005;
          	    float kbT=4.28*(10^(-21));
                float deltaT=0.0000001;
                float P_on;
                float P_off;
                float k_on = k_on0*exp((-1*sigma_on*distance*distance*25*(10^(-14)))/2*kbT);
                float k_off = k_off0*exp(((sigma_off-sigma_on)*distance*distance*25*(10^(-14)))/2*kbT);
            	P_on = (1-exp(-1*k_on*deltaT))*10000;
            	P_off = (1-exp(-1*k_off*deltaT))*10000;
                const hemo::Array<T, 3> rfm = ((P_on*(1-(P_off/10000)))/10000) * br_const * (distance)*5*(10^(-7))  * (dv/distance);
                lParticle.sv.force_repulsion = lParticle.sv.force_repulsion + rfm;
               } else {
               const hemo::Array<T, 3> rfm = br_const * (1/(distance/br_cutoff))  * (dv/distance);
               if (lParticle.sv.cellId==62 || lParticle.sv.cellId==63 || lParticle.sv.cellId==64 || lParticle.sv.cellId==65){  //this should be the cellId of the CTC (or PLT) which can be found in the result folder (.csv file of CTC)
               	lParticle.sv.force_repulsion = lParticle.sv.force_repulsion + rfm;} 
               	cout<<endl<<distance<<endl; \
            }
*/

//pluint cells1=lParticle.sv.celltype;
//cout<<endl<<lParticle.sv.v[0]<<endl<<lParticle.sv.v[1]<<endl<<lParticle.sv.v[2]<<endl<<cells1;
//(lParticle.sv.cellId==0 && nParticle.sv.cellId==1)||( lParticle.sv.cellId==0 && nParticle.sv.cellId==2 )||( lParticle.sv.cellId==0 && nParticle.sv.cellId==3 )||( lParticle.sv.cellId==1 && nParticle.sv.cellId==0)|| (lParticle.sv.cellId==2 && nParticle.sv.cellId==0)|| (lParticle.sv.cellId==3 && nParticle.sv.cellId==0)

//                cout<<"P:"<<P_on<<endl<<"dist2"<<dist2<<"distance:"<<distance<<endl<<"rfm1"<<rfm[1]<<endl; 
// ((P_on*(1-(P_off/100)))/100)
//(lParticle.sv.cellId==0 && nParticle.sv.cellId==1)||( lParticle.sv.cellId==0 && nParticle.sv.cellId==2 )||( lParticle.sv.cellId==0 && nParticle.sv.cellId==3 )||( lParticle.sv.cellId==1 && nParticle.sv.cellId==0)|| (lParticle.sv.cellId==2 && nParticle.sv.cellId==0)|| (lParticle.sv.cellId==3 && nParticle.sv.cellId==0)
// if (lParticle.sv.cellId==62 || lParticle.sv.cellId==63 || lParticle.sv.cellId==64 || lParticle.sv.cellId==65){ \  //this should be the cellId of the CTC (or PLT) which can be found in the result folder (.csv file of CTC)
#define inner_loop \
  bool stochastic_model=1; \
  const int & l_index = grid_index(x,y,z); \
  const int & n_index = grid_index(xx,yy,zz); \
  for (unsigned int i = 0; i < particle_grid_size[l_index];i++){ \
    for (unsigned int j = 0; j < particle_grid_size[n_index];j++){ \
      HemoCellParticle & lParticle = particles[particle_grid[l_index][i]]; \
      HemoCellParticle & nParticle = particles[particle_grid[n_index][j]]; \
      if (&nParticle == &lParticle) { continue; } \
      if (lParticle.sv.cellId == nParticle.sv.cellId) { continue; } \
      const hemo::Array<T,3> dv = lParticle.sv.position - nParticle.sv.position; \
      const T distance = sqrt(dv[0]*dv[0]+dv[1]*dv[1]+dv[2]*dv[2]); \
      float dist2=distance*5*pow(10,-7); \
      float df2=6.40625*pow(10,-9); \
      pluint cells1=lParticle.sv.celltype; \
      pluint cells2=nParticle.sv.celltype; \
      if (distance < r_cutoff) { \
      	if (stochastic_model){ \
      	        float k_on0=1000; \
          	    float k_off0=1; \
          	    float sigma_on=0.0000000005; \
          	    float sigma_off=0.00000000005; \
          	    float kbT=4.28*(pow(10,-21)); \
                float deltaT=0.0000001; \
                float P_on; \
                float P_off; \
                float k_on = k_on0*exp((-1*sigma_on*dist2*dist2)/(2*kbT)); \
                float k_off = k_off0*exp(((sigma_off-sigma_on)*dist2*dist2)/(2*kbT)); \
            	P_on = (1-exp(-1*k_on*deltaT)); \
            	P_off = (1-exp(-1*k_off*deltaT)); \
                const hemo::Array<T, 3> rfm = -1* pow(10,-3)*(P_on*(1-P_off))*(dist2-0.0000001)  * (dv/distance)/df2; \
                if ((cells1==2 && cells2==3 )){ \
                lParticle.sv.force_repulsion = lParticle.sv.force_repulsion + rfm; \
                nParticle.sv.force_repulsion = nParticle.sv.force_repulsion - (rfm/5); \
            } \
            else if ((cells1==3 &&cells2==2 )){ \
                lParticle.sv.force_repulsion = lParticle.sv.force_repulsion + (rfm/5); \
                nParticle.sv.force_repulsion = nParticle.sv.force_repulsion - rfm; \
            } \
               } \
       } \
    } \
  }

void HemoCellParticleField::applyRepulsionForce(bool forced) {
  const T r_const = cellFields->repulsionConstant;
  const T r_cutoff = cellFields->repulsionCutoff;
  if(!pg_up_to_date) {
    update_pg();
  }
  for (HemoCellParticle & particle : particles) {
    particle.sv.force_repulsion = {0.,0.,0.};
  }
  
  for (int x = 0; x < atomicLattice->getNx()-1; x++) {
    for (int y = 0; y < atomicLattice->getNy()-1; y++) {
      for (int z = 0; z < atomicLattice->getNz()-1; z++) {
        //Manual finding, we could make a map, but for now this should be fast enough
        //0, 0, 0 //0, 1, 0 //0, 0, 1 //0, 1, 1
        int xx = x;
        for (int yy = y; yy <= y+1; yy++) {
          for (int zz = z; zz <= z+1; zz++) {
            inner_loop
          }
        }
        //1, 0, 0
        //1, 1, 0
        //1, 0, 1
        //1, 1, 1
        //1, 0,-1
        //1,-1, 0
        //1,-1,-1
        //1, 1,-1
        //1,-1, 1
        xx = x+1;
        for (int yy = y-1; yy <= y+1; yy++) {
          for(int zz = z-1; zz <= z+1; zz++) {
            if (yy < 0) {continue;}
            if (zz < 0) {continue;}
            inner_loop
          }
        }
        xx = x;
        int yy = y+1, zz = z-1;
        if (zz<0) { continue; }
        //0, 1,-1
        inner_loop
        
      }
    }
  }
}

#ifdef INTERIOR_VISCOSITY
void HemoCellParticleField::internalGridPointsMembrane(Box3D domain) {
  // This could be done less complex I guess?
  for (const HemoCellParticle & particle : particles) { // Go over each particle
     if (!(*cellFields)[particle.sv.celltype]->doInteriorViscosity) { continue; }

    for (unsigned int i = 0; i < particle.kernelCoordinates.size(); i++) {
      const hemo::Array<T, 3> latPos = particle.kernelCoordinates[i]-(particle.sv.position-atomicLattice->getLocation());
      const hemo::Array<T, 3> & normalP = particle.normalDirection;

      if (computeLength(latPos) > (*cellFields)[particle.sv.celltype]->mechanics->cellConstants.edge_mean_eq) {continue;}
      
      T dot1 = hemo::dot(latPos, normalP);

      if (dot1 < 0.) {  // Node is inside
        particle.kernelLocations[i]->attributeDynamics((*cellFields)[particle.sv.celltype]->innerViscosityDynamics);
      } else {  // Node is outside
        particle.kernelLocations[i]->attributeDynamics(&atomicLattice->getBackgroundDynamics());
      }
    }
  }
}

// For performance reason, this is only executed once every n iterations to make
// sure that there are no higher viscosity grid points left after substantial movement
void HemoCellParticleField::findInternalParticleGridPoints(Box3D domain) {
  // Reset all the lattice points to the orignal relaxation parameter
  for (plint iX = atomicLattice->getBoundingBox().x0; iX <= atomicLattice->getBoundingBox().x1; iX++) {
    for (plint iY = atomicLattice->getBoundingBox().y0; iY <= atomicLattice->getBoundingBox().y1; iY++) {
      for (plint iZ = atomicLattice->getBoundingBox().z0; iZ <= atomicLattice->getBoundingBox().z1; iZ++) {
        for (HemoCellField * cellfield : cellFields->cellFields) {
          if ( &atomicLattice->get(iX, iY, iZ).getDynamics() == cellfield->innerViscosityDynamics) {
            atomicLattice->get(iX, iY, iZ).attributeDynamics(&atomicLattice->getBackgroundDynamics());
          }
        }
      }
    }
  }

  for (const auto & pair : get_lpc()) { // Go over each cell?
    const int & cid = pair.first;
    const vector<int> & cell = get_particles_per_cell().at(cid);
    const pluint ctype = particles[cell[0]].sv.celltype;

    // Plt and Wbc now have normal tau internal, so we don't have
    // to raycast these particles
    if (!(*cellFields)[ctype]->doInteriorViscosity) {
      continue;
    }
    
    hemo::OctreeStructCell octCell(3, 1, 30,
                                  (*cellFields)[ctype]->mechanics->cellConstants.triangle_list,
                                  particles, cell);

    vector<Cell<T,DESCRIPTOR>*> innerNodes;
    octCell.findInnerNodes(atomicLattice,particles,cell,innerNodes);
    for (Cell<T,DESCRIPTOR>* node : innerNodes) {
      node->attributeDynamics((*cellFields)[ctype]->innerViscosityDynamics);
    }
  }
}
#else
void HemoCellParticleField::findInternalParticleGridPoints(Box3D domain) {
  pcout << "(HemoCellParticleField) (Error) findInternalParticleGridPoints called, but INTERIOR_VISCOSITY not defined, exiting..." << endl;
  exit(1);
}
void HemoCellParticleField::internalGridPointsMembrane(Box3D domain) {
  pcout << "(HemoCellParticleField) (Error) internalGridPointsMembrane called, but INTERIOR_VISCOSITY not defined, exiting..." << endl;
  exit(1);
}
#endif

void HemoCellParticleField::interpolateFluidVelocity(Box3D domain) {
  //Prealloc is nice
  hemo::Array<T,3> velocity;
  plb::Array<T,3> velocity_comp;

  for (HemoCellParticle &particle:particles) {

    //Clever trick to allow for different kernels for different particle types.
    //(*cellFields)[particle.sv.celltype]->kernelMethod(*atomicLattice,particle);

    //We have the kernels, now calculate the velocity of the particles.
    //Palabos developers, sorry for not using a functional...
    velocity = {0.0,0.0,0.0};
    for (pluint j = 0; j < particle.kernelLocations.size(); j++) {
      //Yay for direct access
      particle.kernelLocations[j]->computeVelocity(velocity_comp);
      velocity += (velocity_comp * particle.kernelWeights[j]);
    }
    particle.sv.v = velocity;
  }

}

void HemoCellParticleField::spreadParticleForce(Box3D domain) {
  for( HemoCellParticle &particle:particles) {

    //Clever trick to allow for different kernels for different particle types.
    (*cellFields)[particle.sv.celltype]->kernelMethod(*atomicLattice,particle);

    // Capping force to ensure stability -> NOTE: this introduces error!
#ifdef FORCE_LIMIT
    const T force_mag = norm(particle.sv.force);
    if(force_mag > param::f_limit)
      particle.sv.force *= param::f_limit/force_mag;
#endif

    //Directly change the force on a node , Palabos developers hate this one
    //quick non-functional trick.
    for (pluint j = 0; j < particle.kernelLocations.size(); j++) {
      //Yay for direct access
      particle.kernelLocations[j]->external.data[0] += ((particle.sv.force_repulsion[0] + particle.sv.force[0]) * particle.kernelWeights[j]);
      particle.kernelLocations[j]->external.data[1] += ((particle.sv.force_repulsion[1] + particle.sv.force[1]) * particle.kernelWeights[j]);
      particle.kernelLocations[j]->external.data[2] += ((particle.sv.force_repulsion[2] + particle.sv.force[2]) * particle.kernelWeights[j]);
    }

  }
}

void HemoCellParticleField::populateBoundaryParticles() {
  //qqw2qswws32 <- Greatly appreciated input of Gábor

  for (int x = 0; x < this->atomicLattice->getNx()-1; x++) {
    for (int y = 0; y < this->atomicLattice->getNy()-1; y++) {
      for (int z = 0; z < this->atomicLattice->getNz()-1; z++) {
        if (this->atomicLattice->get(x,y,z).getDynamics().isBoundary()) {
          for (int xx = x-1; xx <= x+1; xx++) {
            if (xx < 0 || xx > this->atomicLattice->getNx()-1) {continue;}
            for (int yy = y-1; yy <= y+1; yy++) {
              if (yy < 0 || yy > this->atomicLattice->getNy()-1) {continue;}
              for (int zz = z-1; zz <= z+1; zz++) {
                if (zz < 0 || zz > this->atomicLattice->getNz()-1) {continue;}
                if (!this->atomicLattice->get(xx,yy,zz).getDynamics().isBoundary()) {
                  boundaryParticles.push_back({x,y,z}); 
                  goto end_inner_loop;
                }
              }
            }
          }
        }
end_inner_loop:;
      }
    }
  }
}

void HemoCellParticleField::applyBoundaryRepulsionForce() {
    if(!pg_up_to_date) {
    update_pg();
  }
 
  const T & br_cutoff = cellFields->boundaryRepulsionCutoff;
  const T & br_const = cellFields->boundaryRepulsionConstant;
  for (Dot3D & b_particle : boundaryParticles) {
    for (int x = b_particle.x-5; x <= b_particle.x+5; x++) {
      if (x < 0 || x > this->atomicLattice->getNx()-1) {continue;}
      for (int y = b_particle.y-5; y <= b_particle.y+5; y++) {
        if (y < 0 || y > this->atomicLattice->getNy()-1) {continue;}
        for (int z = b_particle.z-5; z <= b_particle.z+5; z++) {
          if (z < 0 || z > this->atomicLattice->getNz()-1) {continue;}
          const int & index = grid_index(x,y,z);
          for (unsigned int i = 0 ; i < particle_grid_size[index] ; i++ ) {
            HemoCellParticle & lParticle = particles[particle_grid[index][i]];
            const hemo::Array<T,3> dv = lParticle.sv.position - (b_particle + this->atomicLattice->getLocation()); 
            const T distance = sqrt(dv[0]*dv[0]+dv[1]*dv[1]+dv[2]*dv[2]);
            float dist2=distance*5*pow(10,-7); \
            float df2=6.40625*pow(10,-9); \
            
           // cout<<endl<<lParticle.sv.v[0]<<endl<<lParticle.sv.v[1]<<endl<<lParticle.sv.v[2]<<endl<<cells1;
            
            /*
            float dist2=distance*5*pow(10,-7); \
      float df2=6.40625*pow(10,-9); \
      if (distance < r_cutoff) { \
      	if (stochastic_model){ \
      	        float k_on0=1000; \
          	    float k_off0=0.1; \
          	    float sigma_on=0.0000000005; \
          	    float sigma_off=0.00000000005; \
          	    float kbT=4.28*(pow(10,-21)); \
                float deltaT=0.0000001; \
                float P_on; \
                float P_off; \
                float k_on = k_on0*exp((-1*sigma_on*dist2*dist2)/(2*kbT)); \
                float k_off = k_off0*exp(((sigma_off-sigma_on)*dist2*dist2)/(2*kbT)); \
            	P_on = (1-exp(-1*k_on*deltaT)); \
            	P_off = (1-exp(-1*k_off*deltaT)); \
                const hemo::Array<T, 3> rfm = -1* pow(10,-3) * (P_on*(1-P_off))*dist2  * (dv/distance)/df2; \
                if ((lParticle.sv.cellId==0 && nParticle.sv.cellId==1)||( lParticle.sv.cellId==0 && nParticle.sv.cellId==2 )||( lParticle.sv.cellId==0 && nParticle.sv.cellId==3 )||( lParticle.sv.cellId==1 && nParticle.sv.cellId==0)|| (lParticle.sv.cellId==2 && nParticle.sv.cellId==0)|| (lParticle.sv.cellId==3 && nParticle.sv.cellId==0)){ \
                lParticle.sv.force_repulsion = lParticle.sv.force_repulsion + rfm; \
                nParticle.sv.force_repulsion = nParticle.sv.force_repulsion - rfm; \
            */
            //probabilistic adhesion function (modified by Sina sanvarin@uwaterloo.ca)
            if ((0.2 <distance )&&(distance< br_cutoff)){
            

            bool stochastic_model=1; //turning on the stochastic model of boundary adhesion
            
            if (stochastic_model){
             
           	   // srand(time(0)); //helping to create random numbers
               // int p=rand()%10000; //random numbers between 0 and 99	
      	        float k_on0=800; //unstressed on rate                             	       (will be validated in future) xiao 2017
          	    float k_off0=0.1; //unstressed off rate                            			   (will be validated in future) xiao 2017
          	    float sigma_on=0.0000000005; // on strength                                 			   (will be validated in future) xiao 2017
          	    float sigma_off=0.00000000005; // off strength                              			   (will be validated in future) xiao 2017
          	    float kbT=4.28*(pow(10,-21)); //boltmann constant * temperature                  			   (will be validated in future) xiao 2017
                float deltaT=0.0000001; //time step of the simulation (time interval)        		   (will be validated in future) xiao 2017
                float P_on;
                float P_off;

                float k_on = k_on0*exp((-1*sigma_on*dist2*dist2)/(2*kbT)); //simplification: we used distance instead of (distance - equilibrium spring length )
                float k_off = k_off0*exp(((sigma_off-sigma_on)*dist2*dist2)/(2*kbT)); //simplification: we used distance instead of (distance - equilibrium spring length )

            	P_on = (1-exp(-1*k_on*deltaT));;
            	P_off = (1-exp(-1*k_off*deltaT));;
                 pluint ct=lParticle.sv.celltype;
               //cout << endl << "k_on: " << k_on << "       k_off: "<<k_off<< "      distance: "<< distance << "     "<<  endl; //having a sense about the k_off and k_on values
               //cout << endl << "P_on: " << P_on << "       P_off: "<<P_off<< endl; //having a sense about the P_off and P_on values
                hemo::Array<T, 3> vel1=lParticle.sv.v;
                float vel2 = sqrt(vel1[0]*vel1[0]+vel1[1]*vel1[1]+vel1[2]*vel1[2]);
                
                float momentum= 1.98831*pow(10,-10)*vel2*5*1.2;//mass=1.98831 = 1025(rho)*((4*1e-6)^3)*pi*(4/3)/1382[number of nodes]
                                                                    //1.98831*pow(10,-16)* vel2 *5 [is the coefficient of dx/dt ]* pow(10,-7)
                                                                    //for plt it should be multiply by 0.1 because of the ratio of nodes/volume
                if (ct==3){
               	momentum=momentum*0.2;}

                hemo::Array<T, 3> rep=  (dv/distance) * momentum/df2 ;
                rep[0]=0;
                 
                const hemo::Array<T, 3> rfm = -1* pow(10,-3) * (P_on*(1-P_off))*dist2  * (dv/distance)/df2 ;
                
                cout<<"rep:"<<endl<<rep[1]<<endl<<rep[2]<<endl<<"rfm"<<rfm[1]<<rfm[2]<<endl;

               // if (lParticle.sv.cellId==0 || lParticle.sv.cellId==-1 || lParticle.sv.cellId==1 ||lParticle.sv.cellId==-2 ){  //this should be the cellId of the CTC (or PLT) which can be found in the result folder (.csv file of CTC) [|| lParticle.sv.cellId==1 || lParticle.sv.cellId==2 || lParticle.sv.cellId==3 || lParticle.sv.cellId==-1 || lParticle.sv.cellId==-2 || lParticle.sv.cellId==-3 || lParticle.sv.cellId==-4 || lParticle.sv.cellId==-5 || lParticle.sv.cellId==-6 || lParticle.sv.cellId==-7 || lParticle.sv.cellId==-8]
               	if (distance<3){

                 if (ct==3){
                lParticle.sv.force_repulsion = lParticle.sv.force_repulsion + (rfm/5) + rep; 
                 }else{lParticle.sv.force_repulsion = lParticle.sv.force_repulsion + rfm+ rep;}
            } else{
                 if (ct==3){
                lParticle.sv.force_repulsion = lParticle.sv.force_repulsion + (rfm/5); 
                }else {lParticle.sv.force_repulsion = lParticle.sv.force_repulsion + rfm;}
               	 } 
               	  

               // cout<< endl<< " adhesion forces:  x  "<< rfm[0]<<"  y  "<<rfm[1]<<"   z  "<<rfm[2]<<endl;
              //  }
//


            } 
               
            
            
        
 

             
          
 

            } 
          }
        }
      }
    }
  }
}

void HemoCellParticleField::solidifyCells() {
#ifdef SOLIDIFY_MECHANICS
  for (HemoCellField * type : cellFields->cellFields) {
    ppc_up_to_date = false;
    if(type->doSolidifyMechanics) {
      type->mechanics->solidifyMechanics(get_particles_per_cell(),particles, this->atomicLattice, this->CEPAClattice, type->ctype);
    }
  }
  removeParticles(1);
  for (Dot3D & b_particle : boundaryParticles) {
    for (int x = b_particle.x-1; x <= b_particle.x+1; x++) {
      if (x < 0 || x > this->atomicLattice->getNx()-1) {continue;}
      for (int y = b_particle.y-1; y <= b_particle.y+1; y++) {
        if (y < 0 || y > this->atomicLattice->getNy()-1) {continue;}
        for (int z = b_particle.z-1; z <= b_particle.z+1; z++) {
          if (z < 0 || z > this->atomicLattice->getNz()-1) {continue;}
          const int & index = grid_index(x,y,z);
          for (unsigned int i = 0 ; i < particle_grid_size[index] ; i++ ) {
            HemoCellParticle & lParticle = particles[particle_grid[index][i]];
            const hemo::Array<T,3> dv = lParticle.sv.position - (b_particle + this->atomicLattice->getLocation()); 
            const T distance = sqrt(dv[0]*dv[0]+dv[1]*dv[1]+dv[2]*dv[2]); 
            const hemo::Array<T,9> shearrate;
            if (distance < (*cellFields)[lParticle.sv.celltype]->mechanics->cfg["MaterialModel"]["distanceThreshold"].read<T>()) { 
              lParticle.sv.solidify = true; 
            } 
          }
        }
      }
    }
  }
#else
  hlog << "(HemoCellParticleField) SolidifyCells called but SOLIDIFY_MECHANICS not enabled" << endl;
  exit(1);
#endif
}


HemoCellParticleDataTransfer& HemoCellParticleField::getDataTransfer() {
    return particleDataTransfer;
}
HemoCellParticleDataTransfer const& HemoCellParticleField::getDataTransfer() const {
    return particleDataTransfer;
}

std::string HemoCellParticleField::getBlockName() {
    return std::string("HemoParticleField3D");
}

HemoCellFields* HemoCellParticleField::cellFields=0;
}
