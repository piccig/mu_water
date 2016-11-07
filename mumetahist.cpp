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
   
    class mumetahist : public Colvar {
      //      SwitchingFunction switchingFunction; //instance sw.f
      int gsz,func, Natoms, id, gsz2; 
      double l_func, sigma_LJ, eps_LJ, R_min, R_max, u,  r_cut, r_skin, coff, kT, r0, zdist, eps_Wall, sigma_Wall, r_cut_Wall, wallv;
      bool vshift, issig, quad, nopbcz,widom,gridfile;
      vector<double> LJpar;
      vector<double> Wallpar;
      vector<double> LJr;
      vector<Vector> grid_in;
      vector<double> grid_w;
      vector<int> griddim;
      vector<AtomNumber> at_list;
      
      //histogram variables
      bool rwhist, xir ,dufile, nofirst;
      int Nbin, Nsmall_av, Nlarge_av, Nhist_av,histprint, Ngdr_av,Nbing;
      double Ubin, Ztot, Wmax, bing;
      vector<double> avhist; //histogram time-averaged energy
      vector<double> avhist2; //histogram time-averaged energy
      vector<double> histp; //histogram parameters
      vector<double> shist; //histogram limits
      vector<double> gdrp; //gdr parameters
      vector<double> avgdr; //time-averaged gdr
      vector<double> avxir; //time-averaged xir
      vector<double> xirav; //time-averaged xir

      ofstream dhist;
      ofstream fhist;
      ifstream wfile;
      ofstream fdu;

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

      
      mumetahist(const ActionOptions&);
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

    PLUMED_REGISTER_ACTION(mumetahist,"MUMETAHIST")
    
    void mumetahist::registerKeywords(Keywords& keys){
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
      keys.add("optional","KT","kT, default 0.73167 kJ/mol =>  88 K");
      keys.add("optional","COFF","force sigma cutoff");
      keys.add("optional","XY","xy plane grid");
      keys.add("optional","WALLJ","Wall potential parameters (for LJ walls)");
      keys.add("optional","WALLVAL","Wall potential value");
      keys.addFlag("WIDOM",false,"generates a random uniform grid. ONLY FOR MONITORING.");
      keys.addFlag("NOPBCZ",false,"remove pbc along z");
      keys.remove("NOPBC");
      keys.add("optional","HISTPRINT","how often the energy/gdr histogram is printed (plumed steps)");
      keys.add("optional","HIST","Collect energy histogram");
      keys.addFlag("WEIGHTS",false,"set to TRUE if histogram or gdr has to be reweighed");
      keys.add("optional","GDR","Collect the GDR");
      keys.addFlag("DUFILE",false,"Plot Insertion Energy File");
      keys.addFlag("NOFIRST",false,"first frame is not considered in the histogram and deltaU collection (for restarts)");
      keys.add("optional","SHIST","Collect the histograms only up to the given s value");      
      keys.addFlag("XIR",false,"Collect the local fugacity xi(R). Only with GDR!!!");

    }
    
    mumetahist::mumetahist(const ActionOptions&ao):
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
      
      
      //2D grid option

      zdist=-1; //default, 3d grid
      parse("XY",zdist); 
      if(zdist!=-1){
	if(zdist<0){
	  log.printf("warning xy plane cannot be negative!!! Back to normal grid");
	  zdist = -1;
	}else if(zdist>=0.5){
	  log.printf("warning xy plane cannot be beyond 0.5*Z!!! Set to 0.5 Z");
	  zdist=0.5;
	  griddim[2]=1;
	}else{
	  log.printf("warning xy plane at %f Z",zdist);
	  griddim[2]=2;
	}
      }
      
      
      for(int i=0; i<3; i++) gsz=gsz*griddim[i];
      //log check
      log.printf("gridsize:\t%d\n",gsz);
      
      parseFlag("WIDOM",widom);
      if(widom){
	srand((unsigned)time(0)); 
	log.printf("Random grid generation at every step. No Forces!\n The grid size is %d, 2d grid option neglected.",gsz);
      }
      
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
	log.printf("Error, r0>Rmin, derivable shift to soft-core inactive\n");
	r0=0.0;
      }else{
	log.printf("r0 = 0.0, non-derivable shift to soft-core\n");
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

      dulist.rcut=r_cut;  
      dulist.rskin=r_skin;  
      dulist.step=-1;
      
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
  
      //LJ WALLS parameters
      
      //check for wall parameters
      Wallpar.resize(3);
      Wallpar[0] = 0.0; //If not found set to 0
      Wallpar[1] = 0.0; //LJ parameters 
      Wallpar[2] = 0.0; //LJ parameters 
      parseVector("WALLJ",Wallpar);
      eps_Wall = Wallpar[0];
      sigma_Wall = Wallpar[1];
      r_cut_Wall = Wallpar[2];
      if(eps_Wall!=0)  log.printf("Wall parameters:\teps_w = %.4f\tsigma_w = %.4f\tcut_w = %.4f\n",eps_Wall,sigma_Wall,r_cut_Wall);
      wallv=0;
      parse("WALLVAL",wallv);
      if(wallv!=0){
	eps_Wall = 0; //skip LJ wall calculation
	log.printf("Wall value at the plane:\tphi_W = %.4f\n",wallv);
      }
      
      parseFlag("NOPBCZ",nopbcz);
      if(nopbcz){
	log.printf("No pbcs along z\n");
      }
          
      //kT [KJ/mol] (default = 0.73167 => 88K)
      
      kT=0.73167; 
      parse("KT",kT);
      
      //log check
      log.printf("kT = %.4f\tkJ/mol\n",kT);
      
      //exponential cutoff (default = 20.0)
      coff=20.0; 
      parse("COFF",coff);
      
      //log check
      log.printf("Cut-off = %.4f kT\n",coff);
     

      //open debug file DBG
      //fdbg.open("dbg.dat");
     

      //HISTOGRAM INPUT
      
      parseFlag("NOFIRST",nofirst);
      if(nofirst){
	log.printf("Nofirst option active, neglecting first frame in data collection.\n");
	dulist.step=-2;
      }

      //histogram s filter
      
      

      shist.resize(3);
      shist[0] = -kT*coff;
      shist[1] = kT*coff;
      shist[2] = -5.; //default value for Dumin (in epsilon)
      parseVector("SHIST",shist);
      log.printf("If any histogram is collected, only for %.4f <= s <= %.4f\n",shist[0],shist[1]);
      log.printf("Exp weights normalized with Du_min = %.4f epsilon\n",shist[2]);

      histprint=500; //default (500)
      parse("HISTPRINT",histprint);      
      histp.resize(3);
      histp[2]=-1.; //default hist non active
      

      //Energy histogram
      parseFlag("WEIGHTS",rwhist);
      parseVector("HIST",histp);  //get histogram parameters, min, max, N bins
      if(histp[2]>0.0){
	Nbin=(int)histp[2];
	Ubin=(histp[1]-histp[0])/(1.0*Nbin);
	//log check
	log.printf("energy histogram activated:\tmin:\t%.5f\tmax:\t%.5f\tbins:\t%d\tdU:\t%.5f\n",histp[0],histp[1],Nbin,Ubin);
	fhist.open("pdUins.dat");
	//header fhist
	fhist << setprecision(5);
	for(int i=0; i<Nbin; i++){
	  double duval=histp[0]+(i+0.5)*Ubin;
	  fhist << duval << "\t";
	}
	fhist<< endl;
	dhist.open("insdata.dat"); //file containing the characteristic numbers for the collected histogram
	//header dhist
	dhist << "dulist.step\tNsmall\t<Nsmall>\tNlarge\t<Nlarge>\tNhist\t<Nhist>\tNtot\t" << endl;
	//fdbg<< "files prepared" << endl; fdbg.flush();
	//init average histogram
	avhist.resize(Nbin);
	avhist2.resize(Nbin);
	fill(avhist.begin(), avhist.end(), 0);
	fill(avhist2.begin(), avhist2.end(), 0);
	//fdbg<< "av hist initialized to 0" << endl; fdbg.flush();
	//init numbers
	Nhist_av=0;
	Nsmall_av=0;
	Nlarge_av=0;
	Ztot=0.0;
	
	if(rwhist){
	  wfile.open("weights.dat");
	  Wmax=0.0;
	  while(!wfile.eof()){
	    int step;
	    double bias,coft;
	    wfile >> step >> bias >> coft;
	    if(bias-coft>Wmax) Wmax=bias-coft;
	  }
	  wfile.close();
	  wfile.open("weights.dat");
	  log.printf("weight file read and open. Wmax = %.5f\n",Wmax);
	}
	
      }
      
      gdrp.resize(3);
      gdrp[0]=0.;
      gdrp[1]=0.;
      gdrp[2]=-1.; //default gdr non active
    
      //look for flags
      parseFlag("XIR",xir);

      //GDR
      parseVector("GDR",gdrp);  //get gdr parameters, min, max, N bins
      if(gdrp[2]>0.0){
	Nbing=(int)gdrp[2];
	bing=(gdrp[1]-gdrp[0])/(1.0*Nbing);
	//log check
	log.printf("gdr activated:\tmin:\t%.5f\tmax:\t%.5f\tbins:\t%d\tdr:\t%.5f\n",gdrp[0],gdrp[1],Nbing,bing);
	//init average histogram
	avgdr.resize(Nbing);
	fill(avgdr.begin(), avgdr.end(), 0);
	
	if(xir){
	  avxir.resize(Nbing);
	  fill(avxir.begin(), avxir.end(), 0);
	  xirav.resize(Nbing);
	  fill(xirav.begin(), xirav.end(), 0);
	}
	Ngdr_av=0;
	Ztot=0.0;
	if(histp[2]<=0.0){
	  if(rwhist){
	    wfile.open("weights.dat");
	    Wmax=0.0;
	    while(!wfile.eof()){
	      int step;
	      double bias,coft;
	      wfile >> step >> bias >> coft;
	      if(bias-coft>Wmax) Wmax=bias-coft;
	    }
	    wfile.close();
	    wfile.open("weights.dat");
	    log.printf("weight file read and open. Wmax = %.5f\n",Wmax);
	  }
	}	
      }
      
      //fdbg<< "gdr???" << gdrp[2] << endl; fdbg.flush();
      
      //print DU file?
      parseFlag("DUFILE",dufile);
      if(dufile){
	fdu.open("duins.dat");
	//fdu << "N\t(U(N+1)-U(N))/kT\tstep" <<endl; 
      }


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
      
      //rank=comm.Get_rank(); //Rank of pr
      //stringstream A;
      //A<< "rdbg.dat."<< rank;
      //rdbg.open(A.str().c_str());


      //fdbg.flush(); //DBG
    }

    //fermi switchon (auto-width (rmin-r0)/20)
    double mumetahist::swfon(double z, double Coff){
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
    double mumetahist::swfoff(double z, double Coff){
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
    double mumetahist::dswf(double z, double Coff){
      double dsig;
      if(fabs(z) > Coff){
	dsig=0.0;
      }else{
	dsig=0.5/(1.0+cosh(z));
      }
      return(dsig);
    }
    

    //LJ Potential
    double mumetahist::phi(double rq, double rmin, double rmax, double sigma, double eps, double Vmax, double Vshift){
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
    double mumetahist::dphi(double rq, double r, double rmin, double rmax, double sigma, double eps){
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
    double mumetahist::phis(double rq, double r, double rmin, double rmax, double sigma, double eps, double Vmax, double Vshift, double r0, bool issig, double coff){
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
    double mumetahist::dphis(double rq, double r, double rmin, double rmax, double sigma, double eps, double Vmax, double Vshift, double r0, bool issig, double coff){
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
	  double inw=20./(rmin-r0);
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
    
    double mumetahist::phiq(double rq, double r, double rmin, double rmax, double sigma, double eps, double Vmax, double Vshift, double dmin){
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
    double mumetahist::dphiq(double rq, double r, double rmin, double rmax, double sigma, double eps, double dmin){
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
    

    void mumetahist::du_newlist(vector<AtomNumber> &list, vector<Vector> &grid, dulist_s &dulist)
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

    void mumetahist::du_checklist(vector<AtomNumber> &list, vector<Vector> &grid, dulist_s &dulist)
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
    void mumetahist::calculate()    
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
      
            
      double Vbox=getBox().determinant(); //box volume
      double dVi=0.;  //volume element, in [Vbox] units

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
   
      dVi=1./gsz;
      
      if(gridfile){
	for(i=0; i<gsz; ++i){
	  xgrid[i]=matmul(transpose(getBox()),grid_in[i]);
	}
      }else{ //generate grid
	

	if(widom){ // random grid
	  double r1,r2,r3;
	  for(i=0; i<gsz; i++){
	    r1=static_cast <float> (rand()) / static_cast <float> (RAND_MAX+1.);
	    r2=static_cast <float> (rand()) / static_cast <float> (RAND_MAX+1.);
	    r3=static_cast <float> (rand()) / static_cast <float> (RAND_MAX+1.);
	    gs=Vector(r1,r2,r3);
	    xgrid[i]=matmul(transpose(getBox()),gs);
	  }
	}else{ //regular grid
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
      }
            
      //rdbg << "Grid buildt" << endl; rdbg.flush();
      
      double xvec=0.0;
      
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
	  phiwall=phi(zdist*zdist*LBC[2]*LBC[2], rwallmax, r_cut_Wall, sigma_Wall, eps_Wall, wallmax, wallshift)+phi((1.-zdist)*(1.-zdist)*LBC[2]*LBC[2], rwallmax, r_cut_Wall, sigma_Wall, eps_Wall, wallmax, wallshift);
	}
      }

      //if fixed wallval is given in input
      if(wallv!=0)  phiwall=wallv; //LJ computation inactive

      
      double dmin;
      if(quad) dmin=dphi(R_min*R_min, R_min, 0.0, R_max, sigma_LJ, eps_LJ);

      
      //create neigh lists  
      if(nofirst){
	if(dulist.step==-2){ //first step
	  du_newlist(at_list,xgrid,dulist);
	}else{
	  if(widom){
	    du_newlist(at_list,xgrid,dulist); //recompute lists
	  }else{
	    du_checklist(at_list,xgrid,dulist); //check lists
	  }
	}
      }else{
	if(dulist.step==-1){ //first step
	  du_newlist(at_list,xgrid,dulist);
	}else{
	  if(widom){
	    du_newlist(at_list,xgrid,dulist); //recompute lists
	  }else{
	    du_checklist(at_list,xgrid,dulist); //check lists
	  }
	}
      }

      ++dulist.step; //0 first step, or -1 if nofirst
      //rdbg << "List buildt" << endl; rdbg.flush();

      double Su=0;
      double Du,modrij,modq,zu,fDu, fDuw, dfDu, MDu, ph;
      Vector rij;
      int atomid;
      double beta=1./kT;
      
      //energy array for dufile
      vector<double> bDuvec;
      if(dufile){  
	bDuvec.resize(gsz);
	fill(bDuvec.begin(), bDuvec.end(), 0.);
      }
      


      //rdbg << "Start" << endl; rdbg.flush();
      
      //HISTOGRAM INITIALIZATION
      vector<int> hist;
      if(histp[2]>0.0){  
	hist.resize(Nbin);
	fill(hist.begin(), hist.end(), 0);
      }

      //rdbg << "Inst histogram initialized" << endl; rdbg.flush();
      vector<int> gdr;
      vector<int> gdri;
      vector<double> xirhist;

       
      if(gdrp[2]>0.0){  
	gdr.resize(Nbing);
	fill(gdr.begin(), gdr.end(), 0);
	if(xir) {
	  gdri.resize(Nbing); //for the ith gridpoint
	  xirhist.resize(Nbing);
	  fill(xirhist.begin(), xirhist.end(), 0.0);
	}
      }
            
      int Ntot,Nsmall,Nlarge,Nhist,Ntotgdr,Ngdr;
      double weight, wcav;
      Ntot=0;
      Nsmall=0;
      Nlarge=0;
      Nhist=0;
      Ntotgdr=0;
      Ngdr=0;
      weight=1.0; //if not WEIGHTS
      wcav=0.0;
      
      
      //cycle on space grid
      //rdbg << "Start cycle" << endl; rdbg.flush();
      for(i=rank; i<gsz; i+=stride){

	//fdbg.flush();
	
	//auxiliary array of vectors for derivatives
	vector<Vector> phiprime(dulist.nn[i]);
	//auxiliary array of tesors for derivatives
	vector<Tensor> phit(dulist.nn[i]);

	//z dependent wall 
	if(zdist==-1 && eps_Wall !=0.0) phiwall=phi(xgrid[i][2]*xgrid[i][2], rwallmax, r_cut_Wall, sigma_Wall, eps_Wall, wallmax, wallshift)+phi((LBC[2]-xgrid[i][2])*(LBC[2]-xgrid[i][2]), rwallmax, r_cut_Wall, sigma_Wall, eps_Wall, wallmax, wallshift);
	
	//auxiliary 
	Du=phiwall;
	MDu=0;

	
	
	//xir counters to zero (grid point gdr)
	if(xir) fill(gdri.begin(), gdri.end(), 0);
	
	//fdbg << i << "\t" << xgrid[i][0] << "\t" << xgrid[i][1] << "\t" << xgrid[i][2] << endl; fdbg.flush();
	//cycle over neighbor particles
	
	for(j=0; j<dulist.nn[i]; ++j){
	  atomid = dulist.ni[i][j]; //atomindex
	  //Delta U
	  rij = pbcDistance(xgrid[i],getPosition(atomid));

	  modrij=rij.modulo();
	  modq=modrij*modrij;
	  
	  
	  //Fill gdr histogram
	  if(gdrp[2]>0.0 && dulist.step>=0){
	    if(modrij>=gdrp[0] && modrij<gdrp[1]){
	      int bin= (int)((modrij-gdrp[0])/bing);
	      gdr[bin]++;
	      Ngdr++;
	      if(xir) gdri[bin]++;
	    }
	    Ntotgdr++;
	  }
	  
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
	  
	  
	  //fdbg << i << "\t" << atomid << "\t" << rij.modulo() << "\t" << ph << endl; fdbg.flush();
	}
	

	//rdbg << i << "Du" << Du << endl; rdbg.flush();
	//Fill insertion energy histogram
	if(dulist.step>=0){
	  //print insertion energy in array
	  if(dufile) bDuvec[i] = beta*Du;
	  
	  if(histp[2]>0.0){
	    double bDu = beta*Du;
	    if(bDu<histp[0]){
	      Nsmall++;
	    }else if(bDu>=histp[1]){
	      Nlarge++;
	    }else{
	      int bin= (int)((bDu-histp[0])/Ubin);
	      //fdbg << i << "\t" << bDu << "\t" << Du << "\t" << histp[0] <<  "\t" << bin << endl; fdbg.flush();
	      hist[bin]++;
	      Nhist++;
	    }
	  }
	  Ntot++;
	}
	//rdbg << i << "Histogram filled" << endl; rdbg.flush();
	//exponential (with cut off at coff)
	if(beta*Du>coff){
	  fDu=0.0;
	  dfDu=0.0;
	}else{
	  fDu=exp(-beta*(Du));
	  fDuw=exp(-beta*(Du-shist[2]));
	  dfDu=-beta*fDu;
	}
	
	//rdbg << i << " fDu " << fDu << "\tfDuw " << fDuw << "\t" << xir << endl; rdbg.flush();
	//fill xir histogram
	if(xir && dulist.step>=0){
	  for(int nb=0; nb<Nbing; nb++){
	    if(gdri[nb]>0){
	      xirhist[nb]+=fDuw*gdri[nb];
	    }
	  }
	  wcav+=fDuw;
	}

	
	Su+=fDu*dVi*grid_w[i]; //evaluate integral
	double dSoV=dfDu*dVi*grid_w[i];//Vbox;
	
	//derivatives -- cycle again over neighbors
	//rdbg << i << "Before derivatives" << endl; rdbg.flush();
	for(j=0; j < dulist.nn[i]; ++j){
	  atomid = dulist.ni[i][j]; //atomindex
	  deriv[atomid]+=dSoV*phiprime[j];
	  
	  virial-=dSoV*phit[j];
	}
	vector<Vector>().swap(phiprime);
	vector<Tensor>().swap(phit);
	
      }
      //rdbg << "Compute loop closed" << endl; rdbg.flush();
      comm.Sum(Su);
      comm.Sum(deriv);
      comm.Sum(virial);
      
      if(dufile && dulist.step>=0){
	comm.Sum(bDuvec);
	if(rank == 0){
	  for(int j=0 ; j<gsz; j++){ 
	    fdu << setprecision(8);
	    fdu << (dulist.step)*gsz+j << "\t" << bDuvec[j] << "\t" << dulist.step << endl;
	  }
	}
      }
      




      //rdbg << "Results communicated" << endl; rdbg.flush();
      
      double S;
      if(Su==0){
	S=kT*coff; //When the integral is 0 mu=coff* kT
      }else{
	S=-kT*(std::log(Su));
      }
      setValue(S);
      for(j=0; j < Natoms; ++j)	{
	setAtomsDerivatives(j, -kT/Su*deriv[j]);
      }
      
      //Add total volume derivative
      
      setBoxDerivatives(-kT/Su*virial);//+SuT);
      
      //WEIGHTED DISTRIBUTION
      //read configuration weight
      if(dulist.step>=0){
	if(rwhist){
	  int step;
	  double coft, bias;
	  wfile >> step >> bias >> coft;
	  weight=exp(beta*(bias-coft-Wmax));
	}

	//rdbg << "Configuration weight:  " << weight << endl; rdbg.flush();
	//increase partition function (only if S is less than the threshold)
	if(S>=shist[0] && S<=shist[1]) Ztot+=weight;
      
	//INSERTION ENERGY HISTOGRAM: print and update
	
	if(histp[2]>0.0){
	  
	  //rdbg << "Ready to communicate data" << endl; rdbg.flush();
	  //communicate sum
	  comm.Sum(hist);
	  comm.Sum(Ntot);
	  comm.Sum(Nsmall);
	  comm.Sum(Nlarge);
	  comm.Sum(Nhist);
	  
	  //rdbg << "Test on S:" << S << "\t"<<shist[0]<<"\t"<<shist[1]<< endl; rdbg.flush();
	  if(S>=shist[0] && S<=shist[1]){
	    //Increase averages
	    Nsmall_av+=Nsmall;
	    Nlarge_av+=Nlarge;
	    Nhist_av+=Nhist;
	    for(j=0; j<Nbin; ++j){
	      avhist[j]=(avhist[j]*(Ztot-weight)+(weight*hist[j])/(Ntot*Ubin))/Ztot; //update total histogram
	      avhist2[j]+=hist[j]; //update total histogram
	    }
	  }
	  //rdbg << "Averages updated" << endl; rdbg.flush();
	  //print on files
	  if((dulist.step+1)%histprint == 0 && rank == 0){
	    ofstream favhist;
	    favhist.open("avpdUins.dat"); //file t-averaged histogram
	    fhist << setprecision(5);
	    favhist << setprecision(5);
	    for(j=0; j<Nbin; ++j){
	      fhist << (1.0*hist[j])/(Ntot*Ubin) <<"\t";   // current configuration histogram
	      favhist << (histp[0]+(j+0.5)*Ubin) <<"\t"<< avhist[j] << "\t" << avhist2[j] << endl; //t-averaged histogram
	    }
	    fhist << endl;
	    dhist << setprecision(5);
	    dhist << dulist.step << "\t" << Nsmall << "\t"  << 1.0*Nsmall_av/(dulist.step+1) << "\t"<< Nlarge << "\t"  << 1.0*Nlarge_av/(dulist.step+1)<< "\t" << Nhist << "\t"  << 1.0*Nhist_av/(dulist.step+1) << "\t" << Ntot << "\t" << fixed << setprecision(8) << weight << endl; //hist N data
	    //rdbg << "Files updated" << endl; rdbg.flush();
	    fhist.flush();
	    dhist.flush();
	    favhist.close();
	    //rdbg << "Files Closed" << endl; rdbg.flush();
	    
	  }
	  
	}
      
	
	
	//PAIR DISTRIBUTION HISTOGRAM: print and update
	if(gdrp[2]>0.0){
	  
	  //communicate sum
	  comm.Sum(gdr);
	  comm.Sum(Ntotgdr);
	  comm.Sum(Ngdr);
	  if(xir){
	    comm.Sum(xirhist);
	    comm.Sum(wcav);
	  }
	  
	  //Increase averages (depending on s value)
	  if(S>=shist[0] && S<=shist[1]){
	    
	    Ngdr_av+=Ngdr;
	    for(j=0; j<Nbing; ++j){
	      avgdr[j]=(avgdr[j]*(Ztot-weight)+(weight*gdr[j]/gsz))/Ztot; //update average gdrw
	      if(xir){ // && gdr[j]>0){ 	  //Fill local fugacity array
		xirav[j]=(xirav[j]*(Ztot-weight)+(weight*xirhist[j]/wcav))/Ztot; //update average xir
		if(gdr[j]>0){
		  avxir[j]=(avxir[j]*(Ztot-weight)+(weight*xirhist[j]/gdr[j]))/Ztot; //update average xir
		}else{
		  avxir[j]=(avxir[j]*(Ztot-weight))/Ztot; //update average xir
		}
	      }
	    }
	  }
	  
	  //print on files
	  if((dulist.step+1)%histprint == 0 && rank == 0){
	    ofstream favgdr, favxir;
	    favgdr.open("avgdr.dat"); //file t-averaged histogram
	    favgdr << setprecision(5);
	    if(xir){
	      favxir.open("avxir.dat"); //file t-averaged histogram
	      favxir << setprecision(5);
	    }
	    double integral=0.0;
	    double rhomol=(1.*getNumberOfAtoms())/Vbox;
	    double rhogrid=(1.*gsz)/Vbox;
	    for(j=0; j<Nbing; ++j){
	      double rj0=(gdrp[0]+(j)*bing);
	      double rj1=(gdrp[0]+(j+1)*bing);
	      double rjh=(gdrp[0]+(j+0.5)*bing);
	      double Vj=(4./3.*M_PI*(rj1*rj1*rj1-rj0*rj0*rj0));
	      integral+=avgdr[j];
	      favgdr << rjh <<"\t"<< avgdr[j] << "\t" << avgdr[j]/Vj/rhomol << "\t" << integral*bing/Vj/rhomol << endl; //t-averaged histogram
	      if(xir) favxir << rjh <<"\t"<< avxir[j]/Vj/rhomol << "\t" << xirav[j]/Vj/rhomol << "\t" << wcav << endl; //t-averaged histogram
	    }
	    favgdr.close();
	    if(xir) favxir.close();
	    
	  }
	  
	}
      }
	

    }
  }
}
