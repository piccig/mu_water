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
#include <ctime>

using namespace std;


namespace PLMD{
  namespace colvar{

   
    //+PLUMEDOC COLVAR  
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
   
    class mudel : public Colvar {
      //      SwitchingFunction switchingFunction; //instance sw.f
      int Natoms, id, Ndel, Nsmall_av, Nlarge_av, Nhist_av, Nb, histprint, Nr;
      double sigma_LJ, eps_LJ, R_min, R_max, r_cut, r_skin, coff, kT, r0, Ubin,Du0;
      bool vshift, issig, shuff, rand, quad;
      vector<double> LJpar;
      vector<double> LJr;
      vector<double> avhist; //histogram time-averaged energy
      vector<double> histp; //histogram parameters
      vector<AtomNumber> at_list;
      vector<int> allind;
      vector<int> delind;
      
      // verlet list structure
      
      struct  ljlist_s{
	
	double rcut;
	double rskin;
	int  step;
	vector<int> nn;
	vector<vector<int> > ni;
	vector<Vector> pos;
      } ljlist;
      
    public:

      
      mudel(const ActionOptions&);
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
      void lj_newlist(vector<Vector> &del, ljlist_s &ljlist);
      void lj_checklist(vector<Vector> &del, ljlist_s &ljlist);
      static void registerKeywords(Keywords& keys);
      //debugfile DBG
      ofstream fdbg;
      ofstream rdbg;
      ofstream dhist;
      ofstream fhist;
    };

    PLUMED_REGISTER_ACTION(mudel,"MUDEL")

    void mudel::registerKeywords(Keywords& keys){
      Colvar::registerKeywords(keys);
      keys.add("atoms","GROUP","the group of atoms we are calculating the potential energy for");
      keys.add("optional","NDEL","Number of deletions");
      keys.add("optional","LJPAR","Lennard-Jones Parameters (default:  sigmaLJ = 0.3405 nm\t eps_LJ = 0.996 kJ/mol)");
      keys.add("optional","LJR","Lennard-Jones radii (default: R_min = 0.9*sigmaLJ \t R_max = 1.25 nm");
      keys.add("optional","R_0","For continuous shift to soft-core");
      keys.addFlag("SIG",false,"Set to TRUE if you want sigmoid switch for continuous derivatives (useless if R_0 is not set)");
      keys.addFlag("QUAD",false,"Set to TRUE if for the quadratic soft-core (R_0 will be ignored)");
      keys.addFlag("SHUFFLE",false,"Set to TRUE to get shuffle atom index for deletion");
      keys.addFlag("RAND",false,"Set to TRUE to compute the deletion of different atoms at any step");
      keys.addFlag("VSHIFT",false,"Set to TRUE if you want to shift the LJ potential to be 0 at the cut-off");
      keys.add("optional","R_SKIN","Verlet lists, default r_skin = 1.25 R_max");
      keys.add("optional","R_CUT","Verlet lists, default r_cut = R_max");
      keys.add("optional","DUZERO","Du0 value, to set the scale of the evaluated exponentials.");
      keys.add("optional","KT","kT, default 0.73167 kJ/mol =>  88 K");
      keys.add("optional","COFF","force sigma cutoff");
      keys.remove("NOPBC");
    }
    
    mudel::mudel(const ActionOptions&ao):
      PLUMED_COLVAR_INIT(ao)
      //init bool parameters
      //isnotmoving(false),
      //serial(false)
    {
      //Read atom groups
      
      //      log.printf("Before parse group\n");       log.flush(); //DBG
      
      parseAtomList("GROUP",at_list);
      
      Natoms=at_list.size();
      
      allind.resize(Natoms);
      
      //log check
      log.printf("Number of atoms:\t%d\n",Natoms);
      
      
      srand((unsigned)time(0)); //init rand()
      parse("NDEL",Ndel);
      log.printf("Average of %d deletions\n",Ndel);
      if(Ndel > Natoms){
	log.printf("Warning!!! Ndel > N atoms, Ndel set equal to N atoms.\n");
	Ndel=Natoms;
      }
      
      delind.resize(Ndel);

      

      for(int i=0; i<Natoms; i++){ 
	allind[i]=i; //atom indexes array
      }

      
      parseFlag("SHUFFLE",shuff);
      parseFlag("RAND",rand);
      if(shuff){
	log.printf("Shuffle indexes\n");
	std:: random_shuffle (allind.begin(), allind.end());
	if(rand) log.printf("Re choose at every step the deleted particles\n");
      }
      
      log.printf("1st deletion indexes:\n");
      for(int i=0; i<Ndel; i++){ // deletion indexes array
	delind[i]=allind[i];
	log.printf("%d\t",delind[i]);
      }
      log.printf("\n");
      


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
      
      
      //init Verlet lists
      
      r_cut=R_max; 
      parse("R_CUT",r_cut);
      r_skin=1.25*R_max; 
      parse("R_SKIN",r_skin);
      
      parseFlag("VSHIFT",vshift);
      
      //log check
      log.printf("Verlet list:\tR_cut = %.4f\tR_skin = %.4f\n",r_cut,r_skin);
      log.flush();
      ljlist.rcut=r_cut;  
      ljlist.rskin=r_skin;  
      ljlist.step=0;
      id=0;
      
      //create ljlist.pos, array containing atom positions X
      //vector<Vector> posvec(AtomNumber);
      Vector ze;
      ze.zero();
      
      ljlist.pos.resize(Natoms);
      fill(ljlist.pos.begin(), ljlist.pos.end(), ze);
      
      ljlist.nn.resize(Ndel);
      fill(ljlist.nn.begin(), ljlist.nn.end(), 0);
      //create ljlist.ni array containing the atom index of the neighbors X
      //vector<vector<double> > neighi;
      ljlist.ni.resize(Ndel);
      for (int i = 0; i < Ndel; ++i){
	ljlist.ni[i].resize(Natoms);
	fill(ljlist.ni[i].begin(), ljlist.ni[i].end(),0);
      }
      
      //log check
      log.printf("pos size:\t  %d\n",ljlist.pos.size());
      log.printf("neigh number size:\t  %d\n",ljlist.nn.size());
      log.printf("1st list size:\t  %d\n",ljlist.ni[0].size());
      log.printf("Last list size:\t  %d\n",ljlist.ni[ljlist.nn.size()-1].size());
      log.flush();
      
      //kT [KJ/mol] (default = 0.73167 => 88K)
      
      Du0=0.0; 
      parse("DUZERO",Du0);
      log.printf("DU_0 = %.4f\n",Du0);
         
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
      
      //open debug file DBG
      //fdbg.open("dbg.dat");
      
      
      unsigned int rank; 
      
      rank=comm.Get_rank(); //Rank of pr
      
      //stringstream A;
      //A<< "rdbg.dat."<< rank;
      //rdbg.open(A.str().c_str());
      //fdbg << "DEBUG FILE \n";  //DBG
      //fdbg.flush(); //DBG
      //log.printf("dbg.dat open  \n");
      
    }
    
    
    
    //fermi switchon (auto-width (rmin-r0)/20)
    double mudel::swfon(double z, double Coff){
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
    double mudel::swfoff(double z, double Coff){
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
    double mudel::dswf(double z, double Coff){
      double dsig;
      if(fabs(z) > Coff){
	dsig=0.0;
      }else{
	dsig=0.5/(1.0+cosh(z));
      }
      return(dsig);
    }
    


    //LJ Potential
    double mudel::phi(double rq, double rmin, double rmax, double sigma, double eps, double Vmax, double Vshift){
      double phi, sorq, sor6, sor12;
      if( rq <= rmin*rmin){
	phi=Vmax;
      }else if(rq > rmax*rmax){
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
    double mudel::dphi(double rq, double r, double rmin, double rmax, double sigma, double eps){
      double dphi, sorq, sor6, sor12;
      if( rq < rmin*rmin){
	dphi=0.0;
      }else if(rq > rmax*rmax){
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
    double mudel::phis(double rq, double r, double rmin, double rmax, double sigma, double eps, double Vmax, double Vshift, double r0, bool issig, double coff){
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
    double mudel::dphis(double rq, double r, double rmin, double rmax, double sigma, double eps, double Vmax, double Vshift, double r0, bool issig, double coff){
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
    

        
    double mudel::phiq(double rq, double r, double rmin, double rmax, double sigma, double eps, double Vmax, double Vshift, double dmin){
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
    double mudel::dphiq(double rq, double r, double rmin, double rmax, double sigma, double eps, double dmin){
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


    void mudel::lj_newlist(vector<Vector> &del, ljlist_s &ljlist)
    {
      int i,j;
      Vector rij, test;
      double mod_rij;
  

      for(j=0; j<Natoms; ++j){ //for the check
	//test=getPosition(j);
	ljlist.pos[j] = getPosition(j);
      }
      
      //fprintf(fdbg,"\n start build list cycle \n"); //DBG
      //fflush(fdbg);

      unsigned int stride;
      unsigned int rank; 
      
      stride=comm.Get_size();  //Number of processes
      rank=comm.Get_rank(); //Rank of pr
      
      for(j=rank; j<Ndel; j+=stride) {     // sum over deleted atoms
	ljlist.nn[j]=0; //reset nlistsize
	//fdbg << ind[j] << endl;
	for(i=0; i<Natoms; ++i) {                                           // sum over atoms
	  rij = pbcDistance(ljlist.pos[i],del[j]);
	  mod_rij=rij.modulo();
	  if (mod_rij < ljlist.rskin && i!=delind[j]){ //if distance < rskin and different atom
	    ljlist.ni[j][ljlist.nn[j]]=i; //index in neighlist
	    ljlist.nn[j]++; //increment nn 
	  }	  
	}
      }
      //comm.Sum(ljlist.nn);
      //comm.Sum(ljlist.ni[0]);
    }

    void mudel::lj_checklist(vector<Vector> &del, ljlist_s &ljlist)
    {
      unsigned int j; 
      double dr=(ljlist.rskin-ljlist.rcut)*0.5;
      Vector rij;
        
      for (j=0; j<Natoms; ++j) { 
	rij = pbcDistance(getPosition(j),ljlist.pos[j]); //check position variations
       	if( fabs(rij[0])>dr || fabs(rij[1])>dr || fabs(rij[2])>dr ) { 
	  id++; //rebuild list counter
	  lj_newlist(del,ljlist); 
	  break; 
	}
      }
    } 
    
    
    
    // calculator
    void mudel::calculate()    
    {
      int i,j,k;
      //Tensor virial;
      
    
      //Parallel parameters
      
      unsigned int stride;
      unsigned int rank; 
      
      stride=comm.Get_size();  //Number of processes
      rank=comm.Get_rank(); //Rank of pr
      
      double ledN = 1./Ndel;
      //parallel DBG

      //open rank debug file DBG
     
      //std::string fname = std::string("rdbg.dat.")+to_string(rank);
      //char rnk[10];
      //rnk=to_string(rank);
      //rdbg.open(fname.c_str());
		
      //Box size
      
      vector<double> LBC(3);
      for(i=0;i<3;++i) LBC[i]=getBox()[i][i];
      
      double Vbox=getBox().determinant(); //box volume
      //double dV=Vbox/gsz;
      //fdbg << "Number of atoms\t" << getNumberOfAtoms() << "\n"; //DBG
      //fdbg << "Box edges\t" << LBC[0] <<"\t"<< LBC[1]<<"\t"<< LBC[2] << "\n"; //DBG
      //fdbg << "Box volume "<<Vbox<<"\n";  //DBG
      //fdbg.flush(); //DBG

      //Space grid
      
      //ofstream fgrid; //DBG
      //fgrid.open("grid.dat"); //DBG
      //fgrid2.open("grid2.dat"); //DBG*/
      
      //fdbg << "gsz "<< gsz<<"\n";  //DBG
      //fdbg << "Grid dim "<<griddim[0]<<"\t"<<griddim[1]<<"\t"<<griddim[2]<<"\n";  //DBG
      //fdbg.flush(); //DBG
      
      //Vector dx;
      
      //derivative vector
      vector<Vector> deriv(Natoms);
      //deletion positions
      vector<Vector> delvec(Ndel);
      
      //Vector viriald;
      Tensor virial;
      virial.zero();
      Vector ze;
      ze.zero();
      //init derivatives
      fill(deriv.begin(), deriv.end(), ze);
      
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

      
      double dmin;
      if(quad) dmin=dphi(R_min*R_min, R_min, 0.0, R_max, sigma_LJ, eps_LJ);
      
      
      if(rand){
	std:: random_shuffle (allind.begin(), allind.end()); //reshuffle indexes
	for(i=0; i<Ndel; i++) delind[i]=allind[i];  // rechose deletion indexes
      }
      
      
      for(i=0; i<Ndel; i++){ //find deletion positions
	delvec[i]=getPosition(delind[i]);
      }
      
      //create neigh lists 
      if(ljlist.step==0){ //first step, create
	lj_newlist(delvec,ljlist);
      }else{ // further steps 
	if(rand){ // new deletion list, recreate neigh lists
	  lj_newlist(delvec,ljlist);
	}else{ // might be the same, check
	  lj_checklist(delvec,ljlist);
	}
      }
      ++ljlist.step;
      
      //init CV
      
      //fdbg << ljlist.step <<"\t"<<ljlist.nn[500]<< endl;
      /*fdbg << "Number of atoms\t" << getNumberOfAtoms() << "\t"<< at_list.size()<<endl; //get size of the grid
      fdbg << "Box edges\t" << LBC[0] <<"\t"<< LBC[1]<<"\t"<< LBC[2] <<endl; //DBG
      fdbg << "Box volume "<<Vbox<<endl;  //DBG
      fdbg << "gridsize\t" << gsz << "\t"<< xgrid.size()<<endl; //get size of the grid
      fdbg.flush();*/

      double Su=0;
      double Du,modrij,modq,zu,fDu,dfDu, MDu, ph;
      Vector rij;
      int atomid;
      double beta=1./kT;

      //cycle on space grid
      
      for(i=rank; i<Ndel; i+=stride){
	
	//auxiliary array of vectors for derivatives
	vector<Vector> phiprime(ljlist.nn[i]);
	Vector phiprime_i;
	phiprime_i.zero(); //derivative on deleted molecule
	//auxiliary array of tesors for derivatives
	
	//auxiliary 
	Du=0;
	MDu=0;
	//fprintf(fdbg,"%d \n",i); //DBG
	//fflush(fdbg);
	//fdbg << ljlist.nn[i] << "\t"; //DBG
	//cycle over neighbor particles
	
	for(j=0; j<ljlist.nn[i]; ++j){
	  //nlist.ni[i][j] i: grid point, j: neighbor molecule
	  atomid = ljlist.ni[i][j]; //atomindex
	  //Delta U
	  rij = pbcDistance(delvec[i],getPosition(atomid));
	  modrij=rij.modulo();
	  modq=modrij*modrij;
	  
	  //phi(double rq, double rmin, double rmax, double sigma, double eps, double Vmax)
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
	  //fdbg << setprecision(8) << delind[i] << "\t" << atomid << "\t" << ph << endl;
	  Du+=ph; //update potential difference
	  
	  //dphi(double rq, double r, double rmin, double rmax, double sigma, double eps)
	  phiprime[j] += rij*dph/modrij; //potential derivative
	  phiprime_i -=  rij*dph/modrij; //potential derivative
	  
	  
	  //if(i==10 && atomid%1000 == 0){
	}
	
	//fdbg <<endl;
	//fdbg << setprecision(8) << delind[i] << "\t" << Du << endl;

	//exponential (with cutoff at coff)
	
	if(beta*(Du-Du0)<-coff){
	  fDu=0.0;
	  dfDu=0.0;
	}else{
	  fDu=exp(beta*(Du-Du0));
	  dfDu=beta*fDu;
	}

	Su+=fDu*ledN; //evaluate integral
	double dSoV=dfDu*ledN;//Vbox;
	//fdbg << delind[i] << "\t" << setprecision(8) << scientific << Du << "\t" << fDu << "\t" << Su <<endl;
	//derivatives -- cycle again over neighbors
	deriv[delind[i]]+=dSoV*phiprime_i; //contribution on deleted particle

	for(j=0; j < ljlist.nn[i]; ++j){
	  atomid = ljlist.ni[i][j]; //atomindex
	  deriv[atomid]+=dSoV*phiprime[j];
	  rij = pbcDistance(delvec[i],getPosition(atomid));
	  virial-=dSoV*Tensor(phiprime[j],rij); //Virial component
	  
	}
	vector<Vector>().swap(phiprime);
      }
      
      comm.Sum(Su);
      comm.Sum(deriv);
      comm.Sum(virial);
      
      if(Su==0){
	setValue(-kT*coff+Du0); //When the integral is 0 mu=coff* kT (very unlikely!!!)
      }else{
	setValue(Du0+kT*(std::log(Su)));
      }
      for(j=0; j < Natoms; ++j)	{
	setAtomsDerivatives(j, kT/Su*deriv[j]);
      }
      
      //Add total volume derivative
      
      setBoxDerivatives(kT/Su*virial);//+SuT);
      
    }
  }
}
