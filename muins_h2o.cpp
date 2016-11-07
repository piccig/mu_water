/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2014 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "Colvar.h"
#include "ActionRegister.h"
#include "tools/Communicator.h"
#include "tools/Pbc.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <string>
#include <cmath>
#include <cstdlib>

using namespace std;

namespace PLMD{
  namespace colvar{

   
    //+PLUMEDOC COLVAR DELTAU 
    /*
      Calculate the solute and solvent concentration in a planar shell of the box
      !!! The derivative contains only the outer boundary terms, the resulting bias force is non-conservative !!! Only for CmuMD

      \par Examples


      \verbatim
      d: DIPOLE GROUP=1-10
      PRINT FILE=output STRIDE=5 ARG=5
      \endverbatim
      (see also \ref PRINT)

      \attention 
      If the total charge Q of the group in non zero, then a charge Q/N will be subtracted to every atom,
      where N is the number of atoms. This implies that the dipole (which for a charged system depends
      on the position) is computed on the geometric center of the group.


    */
    //+ENDPLUMEDOC
   
    class mumeta : public Colvar {
      //      SwitchingFunction switchingFunction; //instance sw.f
      int gsz,func, Natoms, id, gsz2;
      double l_func, sigma_LJ, eps_LJ, R_min, R_max, u,  r_cut, r_skin, coff, kT, r0, zdist, eps_Wall, sigma_Wall, r_cut_Wall, wallv;
      bool vshift, issig, quad, nopbcz,gridfile;
      vector<double> LJpar;
      vector<double> Wallpar;
      vector<double> LJr;
      vector<Vector> grid_in;
      vector<double> grid_w;
      vector<int> griddim;
      vector<AtomNumber> at_list;
      double Z0;
      
      // verlet list structure
      
      struct  dulist_s{
	
	double rcut;
	double rskin;
	int  step;
	vector<int> nn;
	vector<vector<int> > ni;
	vector<Vector> posit;
      } dulist;
      
    public:

      
      mumeta(const ActionOptions&);
      virtual void calculate();
      string to_string(int s);
      double swfon(double z, double Coff);
      double swfoff(double z, double Coff);
      double dswf(double z, double Coff);
      double phi(double rq, double rmin, double rmax, double sigma, double eps, double Vmax, double Vshift);
      double phis(double rq, double r, double rmin, double rmax, double sigma, double eps, double Vmax, double Vshift, double r0, bool issig, double coff);
      double phiq(double rq, double r, double rmin, double rmax, double sigma, double eps, double Vmax, double Vshift, double dmin);
      double dphi(double rq, double r, double rmin, double rmax, double sigma, double eps);
      double dphis(double rq, double r, double rmin, double rmax, double sigma, double eps, double Vmax, double Vshift, double r0, bool issig, double coff);
      double dphiq(double rq, double r, double rmin, double rmax, double sigma, double eps, double dmin);
      void du_newlist(vector<AtomNumber> &list, vector<Vector> &grid, dulist_s &dulist);
      void du_checklist(vector<AtomNumber> &list, vector<Vector> &grid, dulist_s &dulist);
      
      static void registerKeywords(Keywords& keys);
      //debugfile DBG
      ofstream fdbg;
      ofstream rdbg;
      ifstream ingrid;
    };

    PLUMED_REGISTER_ACTION(mumeta,"MUMETA")

    void mumeta::registerKeywords(Keywords& keys){
      Colvar::registerKeywords(keys);
      keys.add("atoms","GROUP","the group of atoms we are calculating the dipole moment for");
      keys.add("optional","LJPAR","Lennard-Jones Parameters (default:  sigmaLJ = 0.3405 nm\t eps_LJ = 0.996 kJ/mol)");
      keys.add("optional","LJR","Lennard-Jones radii (default: R_min = 0.9*sigmaLJ \t R_max = 1.25 nm");
      keys.add("optional","R_CUT","Verlet lists, default r_cut = R_max");
      keys.add("optional","R_0","For continuous shift to soft-core");
      keys.addFlag("VSHIFT",false,"Set to TRUE if you want to shift the LJ potential to be 0 at the cut-off");
      keys.addFlag("SIG",false,"Set to TRUE if you want sigmoid switch for continuous derivatives (useless if R_0 is not set)");
      keys.addFlag("QUAD",false,"Set to TRUE if for the quadratic soft-core (R_0 will be ignored)");
      keys.add("optional","R_SKIN","Verlet lists, default r_skin = 1.25 R_max");
      keys.add("compulsory","GRID","space grid");
      keys.addFlag("GRIDFILE",false,"the space grid is read from an input file grid.dat");
      //keys.addFlag("W2",false,"use numeric weight instead of theorical one");
      keys.add("optional","XY","xy plane grid");
      keys.addFlag("NOPBCZ",false,"remove pbc along z");
      keys.add("optional","WALLJ","Wall potential parameters (for LJ walls)");
      keys.add("optional","WALLVAL","Wall potential value");
      keys.add("optional","KT","kT, default 0.73167 kJ/mol =>  88 K");
      keys.add("optional","COFF","force sigma cutoff");
      keys.remove("NOPBC");
    }
    
    mumeta::mumeta(const ActionOptions&ao):
      PLUMED_COLVAR_INIT(ao)
    {
      //Read atom groups
      
      parseAtomList("GROUP",at_list);
      
      Natoms=at_list.size();
      
      //log check
      log.printf("Number of atoms:\t%d\n",Natoms);
      griddim.resize(3);
      parseVector("GRID",griddim);
      
      gsz=1;
      
      //Z Walls and inhomogeneous features

      zdist=-1; //default, 3d grid
      parse("XY",zdist); 
      if(zdist!=-1){
	if(zdist<0){
	  log.printf("warning xy plane cannot be negative!!! Back to normal grid\n");
	  zdist = -1;
	}else if(zdist>=0.5){
	  log.printf("warning xy plane cannot be beyond 0.5*Z!!! Set to 0.5 Z\n");
	  zdist=0.5;
	  griddim[2]=1;
	}else{
	  log.printf("warning xy plane at %f Z\n",zdist);
	  griddim[2]=2;
	}
      }

      for(int i=0; i<3; i++) gsz=gsz*griddim[i];
      parseFlag("GRIDFILE",gridfile);
      if(gridfile){
	//parseFlag("W2",w2instead);
	ingrid.open("grid.dat"); 
	log.printf("Space grid provided from a file, defined grid and xy options not used\n");
	if(!ingrid.is_open()){
	  log.printf("grid.dat not found, using defined grid.\n");
	}else{//read file
	  zdist = -1;
	  gsz=0;
	  gsz2=0;	
	  string buff;
	  while(!ingrid.eof()){ //count lines
	    Vector v;
	    int ind,count;
	    double w1;
	    ingrid >> ind >> v[0] >> v[1] >> v[2] >> w1 >> count;
	    getline(ingrid,buff);
	    if(!ingrid.eof()){//fill arrays
	      grid_in.push_back(v);
	      grid_w.push_back(1.*count/w1); //already divided by the box volume
	      gsz++;
	      gsz2+=count; //count multiple weights
	    }
	  }
	}
	ingrid.close();
      }else{ //equally weighted
	grid_w.resize(gsz);
	fill(grid_w.begin(), grid_w.end(),1.); //already divided by the box volume
      }
      
      /*if(gridfile){
	ofstream dbgrid;
	dbgrid.open("dbgrid.dat"); 
	for(int i=0; i<gsz; i++){
	  dbgrid << setprecision(6);
	  dbgrid << i << "\t" << grid_in[i][0]  << "\t" << grid_in[i][1]  << "\t" << grid_in[i][2]  << "\t" << grid_w[i] << endl;
	}
	dbgrid << gsz2 <<endl;
	dbgrid.close();
	}*/	
      
      //log check
      
      log.printf("gridsize:\t%d\n",gsz);
                      
      //LJ parameters
      
      LJpar.resize(2);
      LJpar[0] = 0.996; //kJ/mol at 88K
      LJpar[1] = 0.3405; //LJ parameters 
      parseVector("LJPAR",LJpar);
      eps_LJ = LJpar[0];
      sigma_LJ = LJpar[1];
      
      LJr.resize(2);
      LJr[0] = 0.9*sigma_LJ; //nm
      LJr[1] = 1.25; //
      parseVector("LJR",LJr);
      R_min = LJr[0];
      R_max = LJr[1];
      
      //log check
      log.printf("LJ parameters:\teps_LJ = %.4f\tsigma_LJ = %.4f\n",eps_LJ,sigma_LJ);
      log.printf("LJ radii:\tR_min = %.4f\tR_max = %.4f\n",R_min,R_max);
      
      //LJ WALL parameters
      
      //check for z wall parameters
      Wallpar.resize(3);
      Wallpar[0] = 0.0; //If not found set to 0
      Wallpar[1] = 0.0; //LJ parameters 
      Wallpar[2] = 0.0; //LJ parameters 
      parseVector("WALLJ",Wallpar);
      eps_Wall = Wallpar[0];
      sigma_Wall = Wallpar[1];
      r_cut_Wall = Wallpar[2];
      if(eps_Wall!=0){
	log.printf("Wall parameters:\teps_w = %.4f\tsigma_w = %.4f\tcut_w = %.4f\n",eps_Wall,sigma_Wall,r_cut_Wall);
      }
      wallv=0;
      parse("WALLVAL",wallv);
      if(wallv!=0){
	eps_Wall = 0; //skip LJ wall calculation
	log.printf("Wall value at the plane:\tphi_W = %.4f\n",wallv);
      }
      
      
      
      //fdbg.open("dbg.dat");
      
      r0 = 0.0;
      parse("R_0",r0);
      if(r0>0.0 && r0<R_min){
	log.printf("derivable shift to soft-core: r0=%.4f\n",r0);
      }else if(r0>=R_min){
	log.printf("Error, r0>Rmin, derivable shift to soft-core inactive");
	r0=0.0;
      }else{
	log.printf("r0 = 0.0, non-derivable shift to soft-core");
      }
      
      parseFlag("SIG",issig);
      if(issig){
	log.printf("Sigmoidal switch\n");
      }else{
	log.printf("Sinusoidal switch\n");
      }
      
      parseFlag("QUAD",quad);
      if(quad){
	log.printf("Quadratic soft-core\n");
	r0= 0.0;
      }

      parseFlag("NOPBCZ",nopbcz);
      if(nopbcz){
	log.printf("No pbcs along z\n");
      }

      //init Verlet lists
      
      r_cut=R_max; 
      parse("R_CUT",r_cut);
      r_skin=1.25*R_max; 
      parse("R_SKIN",r_skin);
      
      parseFlag("VSHIFT",vshift);
      
      //log check
      log.printf("Verlet list:\tR_cut = %.4f\tR_skin = %.4f\n",r_cut,r_skin);
      log.flush();

      dulist.rcut=r_cut;  
      dulist.rskin=r_skin;  
      dulist.step=0;
      id=0;
      //create dulist.pos, array containing atom positions X
      Vector ze;
      ze.zero();

      //init neigh lists
      
      dulist.posit.resize(Natoms);
      fill(dulist.posit.begin(), dulist.posit.end(), ze);
      //create dulist.nn array containing the number of neighbors per gridpoint X
      dulist.nn.resize(gsz);
      fill(dulist.nn.begin(), dulist.nn.end(), 0);
      //create dulist.ni array containing the atom index of the neighbors X
      dulist.ni.resize(gsz);
      for (int i = 0; i < gsz; ++i){
	dulist.ni[i].resize(Natoms);
	fill(dulist.ni[i].begin(), dulist.ni[i].end(),0);
      }
      
      //kT [KJ/mol] (default = 0.73167 => 88K)
      
      kT=0.73167; 
      parse("KT",kT);
      
      //log check
      log.printf("kT = %.4f\tkJ/mol",kT);
      
      //exponential cutoff (default = 20.0)
      coff=20.0; 
      parse("COFF",coff);
      
      //log check
      log.printf("Cut-off = %.4f kT\n",coff);
      
      checkRead();
      addValueWithDerivatives(); 
      setNotPeriodic();
  
      //log atom list
      log.printf("Atom indexes\n");
      for(unsigned int i=0; i<at_list.size(); ++i){
	log.printf("  %d", at_list[i].serial());
      }
      log.printf("  \n");
      requestAtoms(at_list);
      log.printf("Atoms requested  \n");
      
      unsigned int rank; 
      
      rank=comm.Get_rank(); //Rank of pr
      
      //fdbg.flush(); //DBG
    }

    //fermi switchon (auto-width (rmin-r0)/20)
    double mumeta::swfon(double z, double Coff){
      double sig;
      if( z < -Coff){
	sig=0.0;
      }else if(z > Coff){
	sig=1.0;
      }else{
	sig=1.0/(exp(-z)+1.0);
      }
      return(sig);
    }
     
    //fermi switchoff (auto-width (rmin-r0)/20)
    double mumeta::swfoff(double z, double Coff){
      double sig;
      if( z < -Coff){
	sig=0.0;
      }else if(z > Coff){
	sig=1.0;
      }else{
	sig=1.0/(exp(z)+1.0);
      }
      return(sig);
    }
    
    //fermi switchon derivative  
    double mumeta::dswf(double z, double Coff){
      double dsig;
      if(fabs(z) > Coff){
	dsig=0.0;
      }else{
	dsig=0.5/(1.0+cosh(z));
      }
      return(dsig);
    }
    

    //LJ Potential
    double mumeta::phi(double rq, double rmin, double rmax, double sigma, double eps, double Vmax, double Vshift){
      double phi, sorq, sor6, sor12;
      if( rq <= rmin*rmin){
	phi=Vmax;
      }else if(rq >= rmax*rmax){
	phi=0.0;
      }else{
	sorq = sigma*sigma/rq;
	sor6 = sorq * sorq * sorq;
	sor12 = sor6 * sor6;
	phi=4*eps*(sor12-sor6)-Vshift;
      }
      return(phi);
    }

    
    //LJ derivative
    double mumeta::dphi(double rq, double r, double rmin, double rmax, double sigma, double eps){
      double dphi,sorq,sor6,sor12;
      if( rq <= rmin*rmin){
	dphi=0.0;
      }else if(rq >= rmax*rmax){
	dphi=0.0;
      }else{
	sorq = sigma*sigma/rq;
	sor6 = sorq * sorq * sorq;
	sor12 = sor6 * sor6;
	dphi=-24*eps*(2*sor12-sor6)/r;
      }
      return(dphi);
    }
    
    //LJ Potential
    double mumeta::phis(double rq, double r, double rmin, double rmax, double sigma, double eps, double Vmax, double Vshift, double r0, bool issig, double coff){
      double phi, ph, sorq, sor6, sor12;
      if( rq <= r0*r0){
	phi=Vmax;
      }else if(rq >= rmax*rmax){
	phi=0.0;
      }else if(rq > r0*r0 && rq <rmin*rmin){
	sorq = sigma*sigma/rq;
	sor6 = sorq * sorq * sorq;
	sor12 = sor6 * sor6;
	ph=4*eps*(sor12-sor6)-Vshift;
	if(issig){  //fermi switch 
	  double Z0=0.5*(r0+rmin);
	  double w=0.05*(rmin-r0);
	  double z=(r-Z0)/w;
	  double son=swfon(z,coff);
	  double soff=swfoff(z,coff);
	  phi=Vmax*soff+ph*son;
	}else{
	  //sinusoidal switch
	  double arg=0.5*M_PI*(r-r0)/(rmin-r0);
	  double cosine=cos(arg);
	  double sine=sin(arg);
	  phi=Vmax*cosine*cosine+ph*sine*sine;
	}
      }else{
	sorq = sigma*sigma/rq;
	sor6 = sorq * sorq * sorq;
	sor12 = sor6 * sor6;
	phi=4*eps*(sor12-sor6)-Vshift;
      }
      return(phi);
    }

    //LJ derivative
    double mumeta::dphis(double rq, double r, double rmin, double rmax, double sigma, double eps, double Vmax, double Vshift, double r0, bool issig, double coff){
      double dph, ph, dphi, sorq, sor6, sor12;
      if( rq <= r0*r0){
	dphi=0.0;
      }else if(rq >= rmax*rmax){
	dphi=0.0;
      }else if(rq > r0*r0 && rq <rmin*rmin){
	sorq = sigma*sigma/rq;
	sor6 = sorq * sorq * sorq;
	sor12 = sor6 * sor6;
	ph=4*eps*(sor12-sor6)-Vshift;
	dph=-24*eps*(2*sor12-sor6)/r;
	if(issig){ //fermi switch 
	  double Z0=0.5*(r0+rmin);
	  double inw=20/(rmin-r0);
	  double z=(r-Z0)*inw;
	  double son=swfon(z,coff);
	  double ds=dswf(z,coff)*inw;
	  dphi=-Vmax*ds+dph*son+ph*ds;
	}else{
	  //sinusoidal switch
	  double a=M_PI/(rmin-r0);
	  double arg=0.5*a*(r-r0);
	  double sine=sin(arg);
	  double cossin=cos(arg)*sin(arg);
	  dphi=-Vmax*a*cossin+dph*sine*sine+ph*a*cossin;
	}
      }else{
	sorq = sigma*sigma/rq;
	sor6 = sorq * sorq * sorq;
	sor12 = sor6 * sor6;
	dphi=-24*eps*(2*sor12-sor6)/r;
      }
      return(dphi);
    }
    
    double mumeta::phiq(double rq, double r, double rmin, double rmax, double sigma, double eps, double Vmax, double Vshift, double dmin){
      double phi, sorq, sor6, sor12;
      if(rq <= rmin*rmin){
	phi=0.5*dmin*(rq/rmin-rmin)+Vmax;
      }else if(rq >= rmax*rmax){
	phi=0.0;
      }else{
	sorq = sigma*sigma/rq;
	sor6 = sorq * sorq * sorq;
	sor12 = sor6 * sor6;
	phi=4*eps*(sor12-sor6)-Vshift;
      }
      return(phi);
    }
    
    
    //LJ derivative
    double mumeta::dphiq(double rq, double r, double rmin, double rmax, double sigma, double eps, double dmin){
      double dphi, sorq, sor6, sor12;
      if( rq <= rmin*rmin){
	dphi=r*dmin/rmin;
      }else if(rq >= rmax*rmax){
	dphi=0.0;
      }else{
	sorq = sigma*sigma/rq;
	sor6 = sorq * sorq * sorq;
	sor12 = sor6 * sor6;
	dphi=-24*eps*(2*sor12-sor6)/r;
      }
      return(dphi);
    }
    
    

    // verlet list formation
    

    void mumeta::du_newlist(vector<AtomNumber> &list, vector<Vector> &grid, dulist_s &dulist)
    {
      int i,j,gridsize,listsize;
      Vector rij, test;
      double mod_rij;
      
      gridsize=grid.size(); //get size of the grid
      listsize=list.size(); //get the size of the list

      //fill dulist.pos with new list positions

      for(j=0; j<listsize; ++j){
	dulist.posit[j] = getPosition(j);
      }
      
      unsigned int stride;
      unsigned int rank; 
      
      stride=comm.Get_size();  //Number of processes
      rank=comm.Get_rank(); //Rank of pr

      for(j=rank; j<gridsize; j+=stride) {     // sum over grid
	dulist.nn[j]=0; //reset nlistsize
	for(i=0; i<listsize; ++i) {                                           // sum over atoms
	  rij = pbcDistance(dulist.posit[i],grid[j]);
	  if(nopbcz) rij[2] = grid[j][2]-dulist.posit[i][2]; //nopbc z
	  mod_rij=rij.modulo();
	  if (mod_rij < dulist.rskin){ //if distance < rskin
	    dulist.ni[j][dulist.nn[j]]=i; //index in neighlist
	    dulist.nn[j]++; //increment nn 
	  }	  
	}
      }
    }

    void mumeta::du_checklist(vector<AtomNumber> &list, vector<Vector> &grid, dulist_s &dulist)
    {
      unsigned int j; 
      double dr=(dulist.rskin-dulist.rcut)*0.5;
      Vector rij;
        
      for (j=0; j<list.size(); ++j) { 
	rij = pbcDistance(getPosition(j),dulist.posit[j]); //check position variations
	if(nopbcz) rij[2] = dulist.posit[j][2]-getPosition(j)[2]; //nopbc z
       	if( fabs(rij[0])>dr || fabs(rij[1])>dr || fabs(rij[2])>dr ) { 
	  id++; //rebuild list counter
	  du_newlist(list,grid,dulist); 
	  break; 
	}
      }
    } 
    
    // calculator
    void mumeta::calculate()    
    {
      int i,j,k;
      //Parallel parameters
      
      unsigned int stride;
      unsigned int rank; 
      
      stride=comm.Get_size();  //Number of processes
      rank=comm.Get_rank(); //Rank of pr

      //Box size
      
      vector<double> LBC(3);
      for(i=0;i<3;++i) LBC[i]=getBox()[i][i];
      if(dulist.step==0) Z0=LBC[2];
            
      double Vbox=getBox().determinant(); //box volume
      double dVi=1./gsz;  //volume element, in [Vbox] units
      if(gridfile) dVi=1./gsz2; //divided by box volume for imported grid
      //Space grid
      
      ofstream fgrid; //DBG
      
      Vector di;

      for(i=0;i<3;++i){
	di[i]=1./(1.*griddim[i]); //grid2 pace
      }
      
      vector<Vector> xgrid(gsz);
      Vector gs;
      int index;

      
      //derivative vector
      vector<Vector> deriv(getNumberOfAtoms());
      
      //virial contribution
      Tensor virial;
      virial.zero();
      Vector ze;
      ze.zero();
      //init derivatives
      fill(deriv.begin(), deriv.end(), ze);
      
      
      if(gridfile){
	for(i=0; i<gsz; ++i){
	  xgrid[i]=matmul(transpose(getBox()),grid_in[i]);
	}
      }else{
	for(i=0; i<griddim[0]; ++i){
	  for(j=0; j<griddim[1]; ++j){
	    if(zdist == -1){ //normal grid
	      for(k=0; k<griddim[2]; ++k){
		index=griddim[1]*griddim[2]*i+griddim[2]*j+k; //index
		gs=Vector(di[0]*i,di[1]*j,di[2]*k);
		xgrid[index]=matmul(transpose(getBox()),gs);
	      }
	    }else{  //plane
	      index=griddim[1]*griddim[2]*i+griddim[2]*j; //index
	      gs=Vector(di[0]*i,di[1]*j,zdist);
	      xgrid[index]=matmul(transpose(getBox()),gs);
	      //fdbg<< setprecision(8);
	      //fdbg<< xgrid[index][0]<<"\t"<<xgrid[index][1]<<"\t"<<xgrid[index][2]<<endl;
	      if(zdist<0.5){
		index=griddim[1]*griddim[2]*i+griddim[2]*j+1; //index
		gs=Vector(di[0]*i,di[1]*j,1.-zdist);
		xgrid[index]=matmul(transpose(getBox()),gs);
		//fdbg<< setprecision(8);
		//fdbg<< xgrid[index][0]<<"\t"<<xgrid[index][1]<<"\t"<<xgrid[index][2]<<endl;
		//fdbg.flush();
	      }
	    }
	  }
	}
      }
      
      
      
      //initialize phi
      double Vshift;
      if(vshift){
	Vshift=phi(R_max*R_max, 0.0, R_max+1., sigma_LJ, eps_LJ, 0.0, 0.0);
      }else{
	Vshift=0;
      }
      
      double Vmax;
      if(r0>0.0){
	Vmax=phi(r0*r0, 0.0, R_max, sigma_LJ, eps_LJ, 0.0, Vshift);
      }else{
	Vmax=phi(R_min*R_min, 0.0, R_max, sigma_LJ, eps_LJ, 0.0, Vshift);
      }

      double phiwall=0.0;
      if(r_cut_Wall==0) r_cut_Wall=LBC[2]; //choose wall cut-off (even if useless)
      double wallshift=0.0;
      double wallmax=0.0;
      double rwallmax=0.1*sigma_Wall; // standard value
      if(eps_Wall !=0.0){  //if z-wall involved
	wallshift=phi(r_cut_Wall*r_cut_Wall, 0.0, r_cut_Wall+1., sigma_Wall, eps_Wall, 0.0, 0.0); //compute shift
	wallmax=phi(rwallmax*rwallmax, 0.0, r_cut_Wall+1., sigma_Wall, eps_Wall, 0.0, wallshift); //compute wallmax
	if(zdist!=-1){ //if fixed z compute fixed wall contribution
	  phiwall=phi(xgrid[0][2]*xgrid[0][2], rwallmax, r_cut_Wall, sigma_Wall*LBC[2]/Z0, eps_Wall, wallmax, wallshift)+phi((LBC[2]-xgrid[0][2])*(LBC[2]-xgrid[0][2]), rwallmax, r_cut_Wall, sigma_Wall*LBC[2]/Z0, eps_Wall, wallmax, wallshift);
	}
      }

      //if fixed wallval is given in input
      if(wallv!=0)  phiwall=wallv; //LJ computation inactive

      //fdbg<< setprecision(8);
      //fdbg<< dulist.step << "\t" << zdist*zdist*LBC[2]*LBC[2] << "\t" << phiwall <<endl;
      
      double dmin;
      if(quad) dmin=dphi(R_min*R_min, R_min, 0.0, R_max, sigma_LJ, eps_LJ);
      
      
      //create neigh lists  
      
      if(dulist.step==0) du_newlist(at_list,xgrid,dulist);
      ++dulist.step; 
      du_checklist(at_list,xgrid,dulist);
      
      double Su=0;
      double Du,modrij,modq,zu,fDu,dfDu, MDu, ph;
      Vector rij;
      int atomid;
      double beta=1./kT;
      
      //cycle on space grid
      
      for(i=rank; i<gsz; i+=stride){
	
	//fdbg.flush();
	
	//auxiliary array of vectors for derivatives
	vector<Vector> phiprime(dulist.nn[i]);
	//auxiliary array of tesors for derivatives
	vector<Tensor> phit(dulist.nn[i]);
	
	//z dependent wall 
	if(zdist == -1 && eps_Wall !=0.0) phiwall=phi(xgrid[i][2]*xgrid[i][2], rwallmax, r_cut_Wall, sigma_Wall, eps_Wall, wallmax, wallshift)+phi((LBC[2]-xgrid[i][2])*(LBC[2]-xgrid[i][2]), rwallmax, r_cut_Wall, sigma_Wall, eps_Wall, wallmax, wallshift);
	//if(dulist.step == 1){
	//  fdbg<< setprecision(8);
	//  fdbg << xgrid[i][2] << "\t" << phiwall <<endl;
	//}
	//auxiliary 
	Du=phiwall; //zero if no wall is present
	
	MDu=0;
	
	//cycle over neighbor particles
	
	for(j=0; j<dulist.nn[i]; ++j){
	  if(beta*Du>coff)   break;
	  atomid = dulist.ni[i][j]; //atomindex
	  //if(i==1){
	  //  fdbg<< atomid <<"  ";
	  //  fdbg.flush();
	  //}
	  //Delta U
	  rij = pbcDistance(xgrid[i],getPosition(atomid));
	  modrij=rij.modulo();
	  modq=modrij*modrij;
	  
	  double ph,dph;
	  
	  if(r0>0.0){
	    ph=phis(modq, modrij, R_min, R_max, sigma_LJ, eps_LJ, Vmax, Vshift, r0, issig, coff); //potential difference
	    dph= dphis(modq, modrij, R_min, R_max, sigma_LJ, eps_LJ, Vmax, Vshift, r0, issig, coff); //potential derivative
	  }else if(quad){
	    ph=phiq(modq, modrij, R_min, R_max, sigma_LJ, eps_LJ, Vmax, Vshift, dmin);
	    dph=dphiq(modq, modrij, R_min, R_max, sigma_LJ, eps_LJ, dmin);
	    //fdbg setprecision(8) << modrij << "\t" << scientific << ph << "\t" << dph << endl; //DBG
	  }else{
	    ph=phi(modq, R_min, R_max, sigma_LJ, eps_LJ, Vmax, Vshift);
	    dph=dphi(modq, modrij, R_min, R_max, sigma_LJ, eps_LJ);
	  }
	  Du+=ph; //update potential difference
	  phiprime[j] = rij*dph/modrij; //potential derivative
	  phit[j] = Tensor(phiprime[j],rij); //Virial component
	}
	
	//exponential (with cut off at coff)
	if(beta*Du>coff){
	  fDu=0.0;
	  dfDu=0.0;
	}else{
	  fDu=exp(-beta*(Du));
	  dfDu=-beta*fDu;
	}
	
	Su+=fDu*dVi*grid_w[i]; //evaluate integral
	double dSoV=dfDu*dVi*grid_w[i];//Vbox;
	
	//derivatives -- cycle again over neighbors

	for(j=0; j < dulist.nn[i]; ++j){
	  atomid = dulist.ni[i][j]; //atomindex
	  deriv[atomid]+=dSoV*phiprime[j];
	  
	  virial-=dSoV*phit[j];
	}
	vector<Vector>().swap(phiprime);
	vector<Tensor>().swap(phit);
      }
      
      comm.Sum(Su);
      comm.Sum(deriv);
      comm.Sum(virial);
      
      
      if(Su==0){
	setValue(kT*coff); //When the integral is 0 mu=coff* kT
      }else{
	setValue(-kT*(std::log(Su)));
      }
      for(j=0; j < Natoms; ++j)	{
	setAtomsDerivatives(j, -kT/Su*deriv[j]);
      }
      
      //Add total volume derivative
      
      setBoxDerivatives(-kT/Su*virial);//+SuT);
    }
  }
}
