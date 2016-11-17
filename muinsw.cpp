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
   
    class muinsw : public Colvar {
      //      SwitchingFunction switchingFunction; //instance sw.f
      int gsz,func, Natoms, id, gsz2, Nrot, Nconf;
      double l_func, sigma_LJ, eps_LJ, R_min, R_max, u,  r_cut, r_skin, coff, kT, zdist, eps_Wall, sigma_Wall, r_cut_Wall, wallv;
      bool vshift, quad, nopbcz,gridfile,rigid;
      vector<double> LJpar;
      vector<double> Wallpar;
      vector<double> LJr;
      vector<double> grid_w;
      vector<int> griddim;
      vector<AtomNumber> at_list;
      double Z0;
      
      //CHANGE (inserted molecule attributes)
      
      int Nat;
      Vector rcom;
      double Mtot;
      double Maxri;
      vector<double> ins_m;
      vector<Vector> ins_r;
      vector<Vector> ins_s;
      vector<Vector> ins_rot;

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

      
      muinsw(const ActionOptions&);
      virtual void calculate();
      string to_string(int s);
      double swfon(double z, double Coff);
      double swfoff(double z, double Coff);
      double dswf(double z, double Coff);
      double phi(double rq, double rmin, double rmax, double sigma, double eps, double Vmax, double Vshift);
      double phiq(double rq, double r, double rmin, double rmax, double sigma, double eps, double Vmax, double Vshift, double dmin);
      double dphi(double rq, double r, double rmin, double rmax, double sigma, double eps);
      double dphiq(double rq, double r, double rmin, double rmax, double sigma, double eps, double dmin);
      void newlist(vector<Vector> &grid, dulist_s &dulist);
      void checklist(vector<Vector> &grid, dulist_s &dulist);
      void rotate(Vector angle, Vector& vin, Vector& vout);
      
      static void registerKeywords(Keywords& keys);
      //debugfile DBG
      ofstream fdbg;
      ofstream fdbg2;
      ofstream rdbg;
      ifstream inmol;
      ifstream inangles;
    };

    PLUMED_REGISTER_ACTION(muinsw,"MUINSW")

    void muinsw::registerKeywords(Keywords& keys){
      Colvar::registerKeywords(keys);
      keys.add("atoms","GROUP","the group of atoms we are calculating the dipole moment for");
      keys.add("optional","LJPAR","Lennard-Jones Parameters (default:  sigmaLJ = 0.3405 nm\t eps_LJ = 0.996 kJ/mol)");
      keys.add("optional","LJR","Lennard-Jones radii (default: R_min = 0.9*sigmaLJ \t R_max = 1.25 nm");
      keys.add("optional","R_CUT","Verlet lists, default r_cut = R_max");
      keys.addFlag("VSHIFT",false,"Set to TRUE if you want to shift the LJ potential to be 0 at the cut-off");
      keys.addFlag("QUAD",false,"Set to TRUE if for the quadratic soft-core (R_0 will be ignored)");
      keys.add("optional","R_SKIN","Verlet lists, default r_skin = 1.25 R_max");
      keys.add("compulsory","GRID","space grid");
      keys.addFlag("GRIDFILE",false,"the space grid is read from an input file grid.dat");
      keys.addFlag("RIGID",false,"The distances in inserted vector are rigid");
      keys.add("optional","KT","kT, default 0.73167 kJ/mol =>  88 K");
      keys.add("optional","COFF","force sigma cutoff");
      
      //CHANGES
      keys.add("compulsory","INSERT","Molecule file containing the inserted molecule geometry and masses.");
      keys.add("optional","INANG","Angle file containing the rotations for the test insertion of the molecule.");


      keys.remove("NOPBC");
    }
    
    muinsw::muinsw(const ActionOptions&ao):
      PLUMED_COLVAR_INIT(ao)
    {
      //Read atom groups
      
      parseAtomList("GROUP",at_list);
      
      Natoms=at_list.size();
      //CHANGE add var number of molecules 
      //
      //log check
      log.printf("Number of atoms:\t%d\n",Natoms);

      //GRID DEFINITION (Uniform grid)

      griddim.resize(3);
      parseVector("GRID",griddim);
      
      gsz=1;
      
      for(int i=0; i<3; i++) gsz=gsz*griddim[i];
      grid_w.resize(gsz);
      fill(grid_w.begin(), grid_w.end(),1.); //already divided by the box volume
      
      //log check
      
      log.printf("gridsize:\t%d\n",gsz);

      //INSERTION MOLECULE INPUT
      
      string instring;
      parse("INSERT",instring);
      inmol.open(instring.c_str());
      if(!inmol.is_open()){
	log.printf("Insert file not found, ERROR.\n");
      }else{
	Nat=0; //number of atoms equals number of lines
	Mtot=0; //Total mass
	string buff;
	rcom.zero();
	while(!inmol.eof()){ //count lines
	  double Mass;
	  Vector r;
	  inmol >> Mass >> r[0] >> r[1] >> r[2];
	  getline(inmol,buff); 
	  if(!inmol.eof()){//fill arrays and compute rcom
	    ins_m.push_back(Mass);
	    ins_r.push_back(r);
	    rcom += Mass*delta(ins_r[0],ins_r[Nat]) ; 
	    Mtot+=Mass;
	    Nat++;
	  }
	}
	rcom = rcom/Mtot + ins_r[0];
	inmol.close();
      }
      
      //Refer input coordinates to CoM and print log
      log.printf("Inserting %d atom molecule.\n",Nat);
      log.printf("Com in %.4f,%.4f,%.4f and total mass %.4f (uma).\n",rcom[0],rcom[1],rcom[2],Mtot);
      Maxri=0.0;
      for(int i=0; i<Nat; i++){
	log.printf("Atom %d:  M = %.4f \tx = %.4f\ty = %.4f\tz = %.4f\n",ins_m[i],ins_r[i][0],ins_r[i][1],ins_r[i][2]);
	ins_r[i] = delta(rcom,ins_r[i]);
	log.printf("In COM reference: \tx = %.4f\ty = %.4f\tz = %.4f\n",ins_r[i][0],ins_r[i][1],ins_r[i][2]);
	//log.printf("Scaled: \tx = %.4f\ty = %.4f\tz = %.4f\n",ins_s[i][0],ins_s[i][1],ins_s[i][2]);
	if(ins_r[i].modulo()>Maxri) Maxri = ins_r[i].modulo(); // update max distance from COM
      }

      log.printf("Max d from COM = %.4f\n",Maxri);

      //Rotations: generate N ins_r[i] vectors, to be used for the insertion in each gridpoint
      //The ins_r[i] vectors can be generated via input, or via random direction generation
      string inang;
      parse("INANG",inang);
      inangles.open(inang.c_str());
      Nrot=0;
      if(!inangles.is_open()){
	log.printf("Insertion angles file not found, molecule will not be rotated.\n");
      }else{
	string buff;
	log.printf("The following rotations are tested:\n");
	while(!inangles.eof()){ //count lines
	  Vector theta;
	  inangles >> theta[0] >> theta[1] >> theta[2];
	  log.printf("%d alpha %.4f\t beta %.4f \t gamma %.4f\n",Nrot+1,theta[0], theta[1],theta[2]);
	  getline(inangles,buff); 
	  if(!inangles.eof()){ //generate Nrot-th rotation
	    //collect rotation angle arrays
	    ins_rot.push_back(theta);
	    Nrot++;
	  }
	}
	inangles.close();
      }
      
      //CHANGE automatic generation of rotations & random generation (for WIDOM)
      Nconf=Nrot+1; //original configuration plus rotations
      //ROTATE ins_r
      if(Nrot>0){
	for(int i=0; i<Nrot; i++){
	  for(int j=0; j<Nat; j++){
	    Vector vout;
	    rotate(ins_rot[i],ins_r[j],vout);
	    ins_r.push_back(vout);
	    log.printf("Rotation %d - %d x = %.4f\ty = %.4f\tz = %.4f\n",i+1,j,ins_r[i][0],ins_r[i][1],ins_r[i][2]);
	  }
	}
      }
      
      parseFlag("RIGID",rigid);
      if(rigid){
	log.printf("The inserted molecule geometry is rigid.\n");
      }else{
	log.printf("The inserted molecule geometry fluctuates with the box tensor.\n");
      }

      //CHANGE add Coulombic potential parameters and interface for the definition of the interactions

            
      //LJ potential  
      
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
      
      
      //Potential shift and soft-core keywords
      
      parseFlag("VSHIFT",vshift);


      
      parseFlag("QUAD",quad);
      if(quad){
	log.printf("Quadratic soft-core\n");
      }

      //init Verlet lists
      //CHANGE check plumed 2.0 lists and atoms -> molecules 
      r_cut=R_max + Maxri; //molecule size is included 
      parse("R_CUT",r_cut);
      r_skin=1.25*r_cut; 
      parse("R_SKIN",r_skin);
      
      //log check
      log.printf("Verlet list:\tR_cut = %.4f\tR_skin = %.4f\n",r_cut,r_skin);


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
      
      //exponential cutoff (default = 300.0)
      coff=300.0; 
      parse("COFF",coff);
      
      //log check
      log.printf("Cut-off = %.4f kT\n",coff);
      
      checkRead();
      addValueWithDerivatives(); 
      setNotPeriodic();
      log.flush();
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
      
      //fdbg.open("dbg.dat");
      //fdbg2.open("dbg2.dat");
      //fdbg2/fdbg.flush(); //DBG
    }
    
    void muinsw::rotate(Vector angle, Vector& vin, Vector& vout){
      //sin and cos of proper euler angles
      double c[3];
      double s[3];
      for(int i=0; i<3; i++){
	c[i]=cos(angle[i]);
	s[i]=sin(angle[i]);
      }
      Tensor R = Tensor(c[1],-c[2]*s[1],s[1]*s[2],c[0]*s[1],c[0]*c[1]*c[2]-s[0]*s[2],-c[2]*s[0]-c[0]*c[1]*s[2],s[0]*s[1],c[0]*s[2]+c[1]*c[2]*s[0],c[0]*c[2]-c[1]*s[0]*s[2]);
      vout = matmul(R,vin);
    }
    

    //fermi switchon (auto-width (rmin-r0)/20)
    double muinsw::swfon(double z, double Coff){
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
    double muinsw::swfoff(double z, double Coff){
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
    double muinsw::dswf(double z, double Coff){
      double dsig;
      if(fabs(z) > Coff){
	dsig=0.0;
      }else{
	dsig=0.5/(1.0+cosh(z));
      }
      return(dsig);
    }
    

    //LJ Potential
    double muinsw::phi(double rq, double rmin, double rmax, double sigma, double eps, double Vmax, double Vshift){
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
    double muinsw::dphi(double rq, double r, double rmin, double rmax, double sigma, double eps){
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
    
    //LJ potential w/ quadratic sof-core
    double muinsw::phiq(double rq, double r, double rmin, double rmax, double sigma, double eps, double Vmax, double Vshift, double dmin){
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
    
    
    //LJ derivative w/ quadratic sof-core
    double muinsw::dphiq(double rq, double r, double rmin, double rmax, double sigma, double eps, double dmin){
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
    //CHANGE cycle on molecule number (???)    

    void muinsw::newlist(vector<Vector> &grid, dulist_s &dulist)
    {
      int i,j,gridsize,listsize;
      Vector rij, test;
      double mod_rij;
      
      gridsize=grid.size(); //get size of the grid
      listsize=at_list.size(); //get the size of the list

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
	  mod_rij=rij.modulo();
	  if (mod_rij < dulist.rskin){ //if distance < rskin
	    dulist.ni[j][dulist.nn[j]]=i; //index in neighlist
	    dulist.nn[j]++; //increment nn 
	  }	  
	}
      }
    }
    
    //CHANGE cycle on molecule number (???)
    void muinsw::checklist(vector<Vector> &grid, dulist_s &dulist)
    {
      unsigned int j; 
      double dr=(dulist.rskin-dulist.rcut)*0.5;
      Vector rij;
        
      for (j=0; j<at_list.size(); ++j) { 
	rij = pbcDistance(getPosition(j),dulist.posit[j]); //check position variations
       	if( fabs(rij[0])>dr || fabs(rij[1])>dr || fabs(rij[2])>dr ) { 
	  id++; //rebuild list counter
	  newlist(grid,dulist); 
	  break; 
	}
      }
    } 
    
    // calculator
    void muinsw::calculate()    
    {
      int i,j,k,c;
      //Parallel parameters
      
      unsigned int stride;
      unsigned int rank; 
      
      stride=comm.Get_size();  //Number of processes
      rank=comm.Get_rank(); //Rank of pr

      //Box size
      
      vector<double> LBC(3);
      for(i=0;i<3;++i) LBC[i]=getBox()[i][i];
            
      double Vbox=getBox().determinant(); //box volume
      //Space grid
      
      Vector di;

      for(i=0;i<3;++i){
	di[i]=1./(1.*griddim[i]); //grid2 pace
      }
      
      if(dulist.step==0 && !rigid){ //pbc on insertion vectors
	for(i=0;i<Nconf*Nat;i++){
	  ins_s.push_back(getPbc().realToScaled(ins_r[i])); //scaled points
	  //fdbg<< i << "\t" << ins_s[i][0] <<"\t"<< ins_s[i][1]<<"\t"<< ins_s[i][2]<<endl;
	  //for(j=0;j<3;j++) ins_s[i][j]=Tools::pbc(ins_s[i][j]); //within PBC
	  //fdbg<< i << "\t" << ins_s[i][0] <<"\t"<< ins_s[i][1]<<"\t"<< ins_s[i][2]<<endl;
	  //ins_r[i]=getPbc().scaledToReal(ins_s[i]); //Unscale
	  //fdbg<< i << "\t" << ins_r[i][0] <<"\t"<< ins_r[i][1]<<"\t"<< ins_r[i][2]<<endl;
	}
      }

      vector<Vector> xgrid(gsz);
      Vector gs;
      int index;

      
      //derivative vector
      vector<Vector> deriv(getNumberOfAtoms());
      
      //virial contribution
      Tensor zerot;
      zerot.zero();
      Tensor virial;
      virial.zero();
      Vector ze;
      ze.zero();
      //init derivatives
      fill(deriv.begin(), deriv.end(), ze);
      
      //gridpoints

      for(i=0; i<griddim[0]; ++i){
	for(j=0; j<griddim[1]; ++j){
	  for(k=0; k<griddim[2]; ++k){
	    index=griddim[1]*griddim[2]*i+griddim[2]*j+k; //index
	    gs=Vector(di[0]*i,di[1]*j,di[2]*k);
	    xgrid[index]=matmul(transpose(getBox()),gs);
	    //fdbg<<setprecision(10);
	    //fdbg<< xgrid[index][0]<<"\t"<< xgrid[index][1]<<"\t"<< xgrid[index][2]<<endl;
	  }
	}
      }
     
      
      //initialize potential 
      double Vshift;
      if(vshift){
	Vshift=phi(R_max*R_max, 0.0, R_max+1., sigma_LJ, eps_LJ, 0.0, 0.0);
      }else{
	Vshift=0;
      }
      
      double Vmax;
      Vmax=phi(R_min*R_min, 0.0, R_max, sigma_LJ, eps_LJ, 0.0, Vshift);

      double dmin;
      if(quad) dmin=dphi(R_min*R_min, R_min, 0.0, R_max, sigma_LJ, eps_LJ);

      //CHANGE Coulomb initialization
      
      //create neigh lists  
      
      if(dulist.step==0) newlist(xgrid,dulist);
      ++dulist.step; 
      checklist(xgrid,dulist);
      
      double Su=0;
      double Du,modrkj,modq,zu,fDu,dfDu, ph;
      Vector rkj;
      int atomid;
      double beta=1./kT;
      
      double invM=1./gsz/Nconf;  //Number of insertions

      //cycle on space grid

      for(i=rank; i<gsz; i+=stride){
	//auxiliary array of vectors for derivatives (the contribution for each ins configuration is separated)
	vector<Vector> phiprime(dulist.nn[i]);
	//auxiliary array of tensors for derivatives  (the contribution for each ins configuration is separated)
	vector<Tensor> phit(dulist.nn[i]);
	
	for(c=0; c<Nconf; c++){ //loop on configurations (different Du values)      
	  
	  //initialize arrays and Du
	  fill(phiprime.begin(), phiprime.end(), ze);
	  fill(phit.begin(), phit.end(), zerot);
	  Du=0; 	  

	  
	  //loop on inserted atoms (contributing to the same Du)
	  for(k=0; k<Nat; ++k){
	    if(beta*Du>coff)   break; //speed up calculation (it might introduce an error)
	    
	    Vector insx;
	    if(rigid) {
	      insx=xgrid[i]+ins_r[k+c*Nat]; //insertion vector (rigid)
	    }else{
	      insx=xgrid[i]+getPbc().scaledToReal(ins_s[k+c*Nat]); //insertion vector (scaled)
	    }
	    
	    //fdbg2<<setprecision(10);
	    //fdbg2<< insx[0]<<"\t"<< insx[1]<<"\t"<< insx[2]<<endl;
	    
	    
	    //loop over neighbor atoms
	    
	    
	    for(j=0; j<dulist.nn[i]; ++j){
	      if(beta*Du>coff)   break;//speed up calculation (it might introduce an error)
	      
	      atomid = dulist.ni[i][j]; //atomindex
	      //Delta U

	      rkj = pbcDistance(insx,getPosition(atomid)); 
	      //fdbg setprecision(8) << modrij << "\t" << scientific << ph << "\t" << dph << endl; //DBG
	      modrkj=rkj.modulo();
	      modq=modrkj*modrkj;
	      
	      //CHANGE once k and j are associated to atoms, set the sigma,epsilon and Coulombic parameter

	      
	      //CHANGE ph and dph contain all the different interaction contributions

	      double ph,dph;
	      
	      if(quad){
		ph=phiq(modq, modrkj, R_min, R_max, sigma_LJ, eps_LJ, Vmax, Vshift, dmin);
		dph=dphiq(modq, modrkj, R_min, R_max, sigma_LJ, eps_LJ, dmin);
		//fdbg setprecision(8) << modrij << "\t" << scientific << ph << "\t" << dph << endl; //DBG
	      }else{
		ph=phi(modq, R_min, R_max, sigma_LJ, eps_LJ, Vmax, Vshift);
		dph=dphi(modq, modrkj, R_min, R_max, sigma_LJ, eps_LJ);
	      }
	      
	      Du+=ph; //update potential difference
	      Vector phiprimej=rkj*dph/modrkj;

	      phiprime[j] += phiprimej; //potential derivative (k-th atom contribution on j)

	      if(rigid){
		rkj = rkj+ins_r[k+c*Nat]; //refer to CoM
	      }
	      phit[j] += Tensor(rkj,phiprimej); //Virial component (k-th atom contribution on j)
	      //fdbg << i << "\t" << c << "\t" << k <<"\t"<<  j <<"\t"<< setprecision(10) << modrkj << "\t" << ph << "\t" << Du << endl;
	    }
	  }
	  
	  //exponential (with cut off at coff)
	  if(beta*Du>coff){
	    fDu=0.0;
	    dfDu=0.0;
	  }else{
	    fDu=exp(-beta*(Du));
	    dfDu=-beta*fDu;
	  }
	
	  Su+=fDu*invM*grid_w[i]; //evaluate integral
	  double dSoV=dfDu*invM*grid_w[i];
	
	  //derivatives -- cycle again over neighbors
	  //CHANGE cycle on all atoms of the neighborhood list
	  for(j=0; j < dulist.nn[i]; ++j){
	    atomid = dulist.ni[i][j]; //atomindex
	    deriv[atomid]+=dSoV*phiprime[j];
	    virial-=dSoV*phit[j]; //Tensor(phiprime[j],rkj); //Virial component
	  }
	}
	vector<Vector>().swap(phiprime);
	vector<Tensor>().swap(phit);
      }
      
      comm.Sum(Su);
      comm.Sum(deriv);
      comm.Sum(virial);
      
      //Logaritm of the sum
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
