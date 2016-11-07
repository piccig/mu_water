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
   
    class mudelab : public Colvar {
      //      SwitchingFunction switchingFunction; //instance sw.f
      int gsz,func, Natoms, Na, Nb, Ndel, id, rolls, cube_dx, Nsigma;
      double l_func, sigma_A, eps_A,  sigma_B, eps_B, R_min, R_min_A, R_min_B,R_max, u,  r_cut, r_skin, coff, kT, r0, fraction, Du0;
      bool vshift, issig, quad, onecutoff, nolist, shuff, delb, rand;
      vector<double> LJpar_A;
      vector<double> LJpar_B;
      vector<double> LJr;
      vector<int> griddim;
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
	vector<Vector> posit;
      } ljlist;
      
    public:

      
      mudelab(const ActionOptions&);
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
    };

    PLUMED_REGISTER_ACTION(mudelab,"MUDELAB")

    void mudelab::registerKeywords(Keywords& keys){
      Colvar::registerKeywords(keys);
      keys.add("atoms","GROUP","the group of atoms we are calculating the dipole moment for");
      keys.add("optional","LJA","Lennard-Jones Parameters interaction with A specie (default:  sigma = 0.3405 nm\t eps = 0.996 kJ/mol)");
      keys.add("optional","LJB","Lennard-Jones Parameters interaction with B specie (default:  sigma = 0.88 sigma_A\t eps = 0.5 eps_A)");
      keys.add("optional","LJR","Lennard-Jones radii (default: R_min = 0.881*sigma \t R_max = 4.4053 sigma");
      keys.add("optional","R_CUT","Verlet lists, default r_cut = R_max");
      keys.add("optional","R_0","For continuous shift to flat soft-core");
      keys.addFlag("VSHIFT",false,"Set to TRUE if you want to shift the LJ potential to be 0 at the cut-off");
      keys.addFlag("SIG",false,"Set to TRUE if you want sigmoid switch for continuous derivatives (useless if R_0 is not set)");
      keys.addFlag("QUAD",false,"Set to TRUE if for the quadratic soft-core (R_0 will be ignored)");
      keys.addFlag("ONECUT",false,"Set to TRUE for equal cut-off lengths for A and B");
      keys.addFlag("SHUFFLE",false,"Set to TRUE to get shuffle atom index for deletion");
      keys.addFlag("RAND",false,"Set to TRUE to compute the deletion of different atoms at any step");
      keys.addFlag("NOLIST",false,"Set to TRUE to remove neighlist construction");
      keys.add("optional","R_SKIN","Verlet lists, default r_skin = 1.25 R_max");
      keys.add("optional","NDEL","Number of deletions");
      keys.addFlag("DELB",false,"Set to TRUE if the deletion of B kind is considered");
      keys.add("compulsory","NB","number of NB atoms (they MUST BE at the bottom of the configuration!!!)");
      keys.add("optional","KT","kT, default 0.73167 kJ/mol =>  88 K");
      keys.add("optional","COFF","force sigma cutoff");
      keys.add("optional","DUZERO","Du0 value, to set the scale of the evaluated exponentials.");
      keys.remove("NOPBC");
    }
    
    mudelab::mudelab(const ActionOptions&ao):
      PLUMED_COLVAR_INIT(ao)
    {
      //Read atom groups
      
      parseAtomList("GROUP",at_list);
      
      Natoms=at_list.size();

      parse("NB",Nb);
      Na=Natoms-Nb;
      
      
      //log check
      log.printf("Number of atoms:\t%d\n",Natoms);
      log.printf("Number of A:\t%d\n",Na);
      log.printf("Number of B:\t%d\n",Nb);
      
          
      srand((unsigned)time(0)); //init rand()
      parse("NDEL",Ndel);
      log.printf("Average of %d deletions\n",Ndel);
      
      
      parseFlag("DELB",delb);
      if(delb) {
	log.printf("B kind deletion\n");
	if(Ndel > Nb){
	  log.printf("Warning!!! Ndel > N B atoms, Ndel set equal to N B.\n");
	  Ndel=Nb;
	}
      }
      
      if(!delb){
	if(Ndel > Na){
	  log.printf("Warning!!! Ndel > N A atoms, Ndel set equal to N A.\n");
	  Ndel=Na;
	}
      }
      delind.resize(Ndel);
      

      if(delb){
	allind.resize(Nb);
	for(int i=0; i<Nb; i++){ 
	  allind[i]=i+Na; //atom indexes array (only B indexes)
	}
      }else{
	allind.resize(Na);
	for(int i=0; i<Na; i++){ 
	  allind[i]=i; //atom indexes array (only A indexes)
	}
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

      LJpar_A.resize(2);
      LJpar_B.resize(2);
      //default: binary LJ with Ar as A
      
      LJpar_A[0] = 0.996; //kJ/mol
      LJpar_A[1] = 0.3405; //LJ parameters 
      LJpar_B[0] = 1.5*eps_A; 
      LJpar_B[1] = 0.8*sigma_A;

      parseVector("LJA",LJpar_A);
      eps_A = LJpar_A[0];
      sigma_A = LJpar_A[1];
      
      parseVector("LJB",LJpar_B);
      eps_B = LJpar_B[0];
      sigma_B = LJpar_B[1];
      
            
      parseFlag("ONECUT",onecutoff);
      
      if(!onecutoff){
	LJr.resize(2);
	R_min = 0.88106; // sigma units
	R_max = 4.4053; //
	parseVector("LJR",LJr);
	R_min = LJr[0];
	R_max = LJr[1];
      }else{ //one R_max cut-off (different manual R_min)
	LJr.resize(3);
	R_min_A = 0.3; // nm
	R_min_B = 0.8*R_min_A;
	R_max = 2.5*sigma_A; //
	parseVector("LJR",LJr);
	R_min_A = LJr[0];
	R_min_B = LJr[1];
	R_max = LJr[2];
      }


      //log check
      log.printf("LJ parameters:\teps_A = %.4f\tsigma_A = %.4f\n",eps_A,sigma_A);
      log.printf("LJ parameters:\teps_B = %.4f\tsigma_B = %.4f\n",eps_B,sigma_B);
      if(!onecutoff){
	log.printf("LJ radii:\tR_min = %.4f\tR_max = %.4f\t in sigma units\n",R_min,R_max);
      }else{
	log.printf("LJ radii:\tR_min_A = %.4f\tR_min_B = %.4f\tR_max = %.4f\t in nm\n",R_min_A,R_min_B,R_max);
      }
      
      r0 = 0.0;
      parse("R_0",r0);
      if(r0>0.0 && r0<1.0){
	log.printf("derivable shift to soft-core: r0=%.4f R_min\n",r0);
      }else if(r0>=1.0){
	log.printf("Error, r0>Rmin, derivable shift to soft-core inactive");
	r0=0.0;
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

      parseFlag("NOLIST",nolist);

      //init Verlet lists
      if(!onecutoff){
	r_cut=R_max*sigma_A;
	if(r_cut<R_max*sigma_B) r_cut=R_max*sigma_B; //set the largest R_max
      }else{
	r_cut=R_max;
      }
      parse("R_CUT",r_cut);
      r_skin=1.25*r_cut; 
      parse("R_SKIN",r_skin);
      
      parseFlag("VSHIFT",vshift);
  
      double Lbox=0.0;
      //log check
      if(!nolist){
	log.printf("Verlet list:\tR_cut = %.4f\tR_skin = %.4f\n",r_cut,r_skin);
      }else{
	log.printf("Verlet list inactive (r_skin set to box-length)\n");
      }
      
      log.flush();

      ljlist.rcut=r_cut;  
      ljlist.rskin=r_skin;  
      ljlist.step=0;
      id=0;
      //create ljlist.pos, array containing atom positions X
      Vector ze;
      ze.zero();

      //init neigh lists
      
      ljlist.posit.resize(Natoms);
      fill(ljlist.posit.begin(), ljlist.posit.end(), ze);
      //create ljlist.nn array containing the number of neighbors per gridpoint X
      ljlist.nn.resize(Ndel);
      fill(ljlist.nn.begin(), ljlist.nn.end(), 0);
      //create ljlist.ni array containing the atom index of the neighbors X
      ljlist.ni.resize(Ndel);
      for (int i = 0; i < Ndel; ++i){
	ljlist.ni[i].resize(Natoms);
	fill(ljlist.ni[i].begin(), ljlist.ni[i].end(),0);
      }
      
      //log check
      log.printf("pos size:\t  %d\n",ljlist.posit.size());
      log.printf("neigh number size:\t  %d\n",ljlist.nn.size());
      log.printf("1st list size:\t  %d\n",ljlist.ni[0].size());
      log.printf("Last list size:\t  %d\n",ljlist.ni[ljlist.nn.size()-1].size());
      log.flush();


      //kT [KJ/mol] (default = 0.73167 => 88K)
      
      Du0=0.0; 
      parse("DUZERO",Du0);
      log.printf("DU_0 = %.4f\n",Du0);
      
      kT=0.73167; 
      parse("KT",kT);
      
      //log check
      log.printf("kT = %.4f\tkJ/mol\n",kT);
      
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
      //fdbg.open("dbg.dat");
      //fdbg.flush(); //DBG
    }

    //fermi switchon (auto-width (rmin-r0)/20)
    double mudelab::swfon(double z, double Coff){
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
    double mudelab::swfoff(double z, double Coff){
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
    double mudelab::dswf(double z, double Coff){
      double dsig;
      if(fabs(z) > Coff){
	dsig=0.0;
      }else{
	dsig=0.5/(1.0+cosh(z));
      }
      return(dsig);
    }
    

    //LJ Potential
    double mudelab::phi(double rq, double rmin, double rmax, double sigma, double eps, double Vmax, double Vshift){
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
    double mudelab::dphi(double rq, double r, double rmin, double rmax, double sigma, double eps){
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
    double mudelab::phis(double rq, double r, double rmin, double rmax, double sigma, double eps, double Vmax, double Vshift, double r0, bool issig, double coff){
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
    double mudelab::dphis(double rq, double r, double rmin, double rmax, double sigma, double eps, double Vmax, double Vshift, double r0, bool issig, double coff){
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
    
    double mudelab::phiq(double rq, double r, double rmin, double rmax, double sigma, double eps, double Vmax, double Vshift, double dmin){
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
    double mudelab::dphiq(double rq, double r, double rmin, double rmax, double sigma, double eps, double dmin){
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
    

    void mudelab::lj_newlist(vector<Vector> &del, ljlist_s &ljlist)
    {
      int i,j;
      Vector rij, test;
      double mod_rij;
      
      //fill ljlist.pos with new list positions

      for(j=0; j<Natoms; ++j){
	ljlist.posit[j] = getPosition(j);
      }
      
      unsigned int stride;
      unsigned int rank; 
      
      stride=comm.Get_size();  //Number of processes
      rank=comm.Get_rank(); //Rank of pr

       
      for(j=rank; j<Ndel; j+=stride) {     // sum over grid
	ljlist.nn[j]=0; //reset nlistsize
	for(i=0; i<Natoms; ++i) {                                           // sum over atoms
	  rij = pbcDistance(ljlist.posit[i],del[j]);
	  mod_rij=rij.modulo();
	  if (mod_rij < ljlist.rskin && i!=delind[j]){ //if distance < rskin and different atom
	    ljlist.ni[j][ljlist.nn[j]]=i; //index in neighlist
	    ljlist.nn[j]++; //increment nn 
	  }	  
	}
      }
    }

    void mudelab::lj_checklist(vector<Vector> &del, ljlist_s &ljlist)
    {
      unsigned int j; 
      double dr=(ljlist.rskin-ljlist.rcut)*0.5;
      Vector rij;
        
      for (j=0; j<Natoms; ++j) { 
	rij = pbcDistance(getPosition(j),ljlist.posit[j]); //check position variations
       	if( fabs(rij[0])>dr || fabs(rij[1])>dr || fabs(rij[2])>dr ) { 
	  id++; //rebuild list counter
	  lj_newlist(del,ljlist); 
	  break; 
	}
      }
    } 
    
    // calculator
    void mudelab::calculate()    
    {
      int i,j,k;
      //Parallel parameters
      
      unsigned int stride;
      unsigned int rank; 
     
      stride=comm.Get_size();  //Number of processes
      rank=comm.Get_rank(); //Rank of pr
      
  
      //Box size
      double Lmax=0.0;
      vector<double> LBC(3);
      for(i=0;i<3;++i){
	LBC[i]=getBox()[i][i];
	if(LBC[i]>Lmax) Lmax=LBC[i];
      }
      if(nolist) ljlist.rskin=Lmax;
      //fdbg << ljlist.rskin << endl;
            
      double Vbox=getBox().determinant(); //box volume
      double dVi=0.;  //volume element, in [Vbox] units
      double ledN = 1./Ndel;
      //Space grid
      
      //derivative vector
      vector<Vector> deriv(Natoms);
      //deletion positions
      vector<Vector> delvec(Ndel);
      
      //virial contribution
      Tensor virial;
      virial.zero();
      Vector ze;
      ze.zero();
      //init derivatives
      fill(deriv.begin(), deriv.end(), ze);
         
      //scale potential cut-off parameters
      
      double R_max_A,R_max_B;
      if(!onecutoff){
	R_max_A=R_max*sigma_A;
	R_max_B=R_max*sigma_B;
	
	R_min_A=R_min*sigma_A;
	R_min_B=R_min*sigma_B;
      }else{
	R_max_A=R_max;
	R_max_B=R_max;
      }
	
      double r0_A=r0*R_min_A;
      double r0_B=r0*R_min_B;
      
      
      //initialize phi (BIN: values for each interaction)
      double Vshift_A, Vshift_B;
      if(vshift){
	Vshift_A=phi(R_max_A*R_max_A, 0.0, R_max_A+1., sigma_A, eps_A, 0.0, 0.0);
	Vshift_B=phi(R_max_B*R_max_B, 0.0, R_max_B+1., sigma_B, eps_B, 0.0, 0.0);
      }else{
	Vshift_A=0;
	Vshift_B=0;
      }
      
      double Vmax_A, Vmax_B;
      if(r0>0.0){
	Vmax_A=phi(r0_A*r0_A, 0.0, R_max_A, sigma_A, eps_A, 0.0, Vshift_A);
	Vmax_B=phi(r0_B*r0_B, 0.0, R_max_B, sigma_B, eps_B, 0.0, Vshift_B);
      }else{
	Vmax_A=phi(R_min_A*R_min_A, 0.0, R_max_A, sigma_A, eps_A, 0.0, Vshift_A);
	Vmax_B=phi(R_min_B*R_min_B, 0.0, R_max_B, sigma_B, eps_B, 0.0, Vshift_B);
      }

      double dmin_A,dmin_B;
      if(quad){
	dmin_A=dphi(R_min_A*R_min_A, R_min_A, 0.0, R_max_A, sigma_A, eps_A);
	dmin_B=dphi(R_min_B*R_min_B, R_min_B, 0.0, R_max_B, sigma_B, eps_B);
      }

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
	}else if(!nolist){ // might be the same, check only if it doesn't contain all atoms (nolist active)
	  lj_checklist(delvec,ljlist);
	}
      }
      ++ljlist.step;

     
      double Su=0;
      double Du,modrij,modq,zu,fDu,dfDu, MDu, ph;
      Vector rij;
      int atomid;
      double beta=1./kT;
      double eps, sigma, Vshift, Vmax, R_min, R_max, r0,dmin;
      
      //cycle on space grid
      
      for(i=rank; i<Ndel; i+=stride){

	//fdbg.flush();
	
	//auxiliary array of vectors for derivatives
	vector<Vector> phiprime(ljlist.nn[i]);
	Vector phiprime_i;
	phiprime_i.zero(); //derivative on deleted molecule
		
	//auxiliary 
	Du=0;
	MDu=0;

	//cycle over neighbor particles
	
	//fdbg << i << "\t" ; //DBG
	
	//fdbg.flush();
	
	for(j=0; j<ljlist.nn[i]; ++j){
	  
	  atomid = ljlist.ni[i][j]; //atomindex
	  
	  if(atomid<Na){ //select interaction
	    eps=eps_A;
	    sigma=sigma_A;
	    Vshift=Vshift_A;
	    Vmax=Vmax_A;
	    R_min=R_min_A;
	    R_max=R_max_A;
	    r0=r0_A;
	    dmin=dmin_A;
	  }else{
	    eps=eps_B;
	    sigma=sigma_B;
	    Vshift=Vshift_B;
	    Vmax=Vmax_B;
	    R_min=R_min_B;
	    R_max=R_max_B;
	    r0=r0_B;
	    dmin=dmin_B;
	  }

	  //Delta U
	  rij = pbcDistance(delvec[i],getPosition(atomid));
	  modrij=rij.modulo();
	  modq=modrij*modrij;
	  double ph,dph;
	  
	  if(r0>0.0){
	    ph=phis(modq, modrij, R_min, R_max, sigma, eps, Vmax, Vshift, r0, issig, coff); //potential difference
	    dph= dphis(modq, modrij, R_min, R_max, sigma, eps, Vmax, Vshift, r0, issig, coff); //potential derivative
	  }else if(quad){
	    ph=phiq(modq, modrij, R_min, R_max, sigma, eps, Vmax, Vshift, dmin);
	    dph=dphiq(modq, modrij, R_min, R_max, sigma, eps, dmin);
	    //fdbg setprecision(8) << modrij << "\t" << scientific << ph << "\t" << dph << endl; //DBG
	  }else{
	    ph=phi(modq, R_min, R_max, sigma, eps, Vmax, Vshift);
	    dph=dphi(modq, modrij, R_min, R_max, sigma, eps);
	  }
	  Du+=ph; //update potential difference
	  phiprime[j] += rij*dph/modrij; //potential derivative
	  phiprime_i -= rij*dph/modrij;
	}
	
	//fdbg << Du << endl;

	//exponential (with cut off at coff)
	if(beta*(Du-Du0)<-coff){
	  fDu=0.0;
	  dfDu=0.0;
	}else{
	  fDu=exp(beta*(Du-Du0));
	  dfDu=beta*fDu;
	}
	
	Su+=fDu*ledN; //evaluate integral
	double dSoV=dfDu*ledN;//Vbox;
	//fdbg << delind[i] << "\t" << setprecision(8) << scientific << Du << "\t" << fDu << "\t" << Su <<"\t"<< ledN<< "\t" <<endl;
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
