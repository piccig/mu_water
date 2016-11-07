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
   
    class mumetaswhist : public Colvar {
      //      SwitchingFunction switchingFunction; //instance sw.f
      int gsz,func, Natoms, id, gsz2;
      double l_func, sigma_LJ, eps_LJ, R_min, R_max, u,  r_cut, r_skin, vcoff, scoff, kT, r0, zdist, eps_Wall, sigma_Wall, r_cut_Wall, wallv;
      bool vshift, issig, quad, nopbcz,gridfile,widom;
      vector<double> LJpar;
      vector<double> Wallpar;
      vector<double> LJr;
      vector<Vector> grid_in;
      vector<double> grid_w;
      vector<int> griddim;
      vector<AtomNumber> at_list;
      double Z0;

      
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


      //SW parameters
      double A, B, gamma, lambda, epsilon, sigma, costh0, a;
      int p,q;
      
      // verlet list structure
      
      struct  swlist_s{
	
	double rcut;
	double rskin;
	int  step;
	vector<int> nn;
	vector<vector<int> > ni;
	vector<Vector> posit;
      } swlist;
      
    public:

      
      mumetaswhist(const ActionOptions&);
      virtual void calculate();
      string to_string(int s);
      double theta_ab(Vector av, Vector bv);
      double dtheta_ab(Vector av, Vector bv, Vector& dtheta_a, Vector& dtheta_b);
      double h_ab(Vector av, Vector bv, Vector& dh_a, Vector& dh_b);
      double v2(double r, double rmin, double rmax, double Vmax);
      double dv2(double r, double rmin, double rmax, double vr);
      double v2q(double r, double rmin, double rmax, double Vmax, double dv0);
      double dv2q(double r, double rmin, double rmax, double vr, double dv0);
      double v3(Vector rA, Vector ri, Vector rj, double rmax, Vector& dv3_dri, Vector& dv3_drj);
      double phi(double rq, double rmin, double rmax, double sigma, double eps, double Vmax, double Vshift);
      double phiq(double rq, double r, double rmin, double rmax, double sigma, double eps, double Vmax, double Vshift, double dmin);
      double dphi(double rq, double r, double rmin, double rmax, double sigma, double eps);
      double dphiq(double rq, double r, double rmin, double rmax, double sigma, double eps, double dmin);
      void sw_newlist(vector<Vector> &grid, swlist_s &swlist);
      void sw_checklist(vector<Vector> &grid, swlist_s &swlist);
      
      static void registerKeywords(Keywords& keys);
      //debugfile DBG
      ofstream fdbg;
      ofstream rdbg;
      ifstream ingrid;
    };

    PLUMED_REGISTER_ACTION(mumetaswhist,"MUMETASWHIST")

    void mumetaswhist::registerKeywords(Keywords& keys){
      Colvar::registerKeywords(keys);
      keys.add("atoms","GROUP","the group of atoms of the system");
      keys.add("optional","ASW","A parameter of SW potential");
      keys.add("optional","BSW","B parameter of SW potential");
      keys.add("optional","PSW","P parameter of SW potential");
      keys.add("optional","QSW","Q parameter of SW potential");
      keys.add("optional","GAMMASW","gamma parameter of SW potential");
      keys.add("optional","ACUTSW","a parameter of SW potential");
      keys.add("optional","CTHETASW","cos(theta0) parameter of SW potential");
      keys.add("optional","LAMBDASW","lambda parameter of SW potential");
      keys.add("optional","EPSSW","epsilon parameter of SW potential");
      keys.add("optional","SIGMASW","sigma parameter of SW potential");
      keys.addFlag("QUAD",false,"Set to TRUE if for the quadratic soft-core (at R_0)");
      keys.add("optional","R_0","2-body term potential inner cut-off");
      keys.add("optional","R_CUT","Verlet lists, default r_cut = 2*a*sigma");
      keys.add("optional","R_SKIN","Verlet lists, default r_skin = 2.5 a*sigma");
      keys.add("compulsory","GRID","space grid");
      keys.addFlag("GRIDFILE",false,"the space grid is read from an input file grid.dat");
      keys.addFlag("WIDOM",false,"generates a random uniform grid. NO FORCES COMPUTED, only for monitoring.");
      keys.add("optional","XY","xy plane grid");
      keys.addFlag("NOPBCZ",false,"remove pbc along z");
      keys.add("optional","WALLJ","Wall potential parameters (for LJ walls)");
      keys.add("optional","WALLVAL","Wall potential value");
      keys.add("optional","KT","kT, default 2.49434 kJ/mol =>  300 K");
      keys.add("optional","SCOFF","general exponential cutoff (default=300)");
      keys.add("optional","VCOFF","potential exponentials cutoff (default=32)");
      
      //histogram part
      keys.add("optional","HISTPRINT","how often the energy/gdr histogram is printed (plumed steps)");
      keys.add("optional","HIST","Collect energy histogram");
      keys.add("optional","GDR","Collect the GDR");
      keys.addFlag("XIR",false,"Collect the local fugacity xi(R). Only with GDR!!!");
      keys.addFlag("WEIGHTS",false,"set to TRUE if histogram or gdr has to be reweighed");
      keys.addFlag("DUFILE",false,"Plot Insertion Energy File");
      keys.addFlag("NOFIRST",false,"first frame is not considered in the histogram and deltaU collection (for restarts)");
      keys.add("optional","SHIST","Collect the histograms only up to the given s value");      
      keys.remove("NOPBC");
    }
    
    mumetaswhist::mumetaswhist(const ActionOptions&ao):
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
      
      parseFlag("WIDOM",widom);
      if(widom){
	srand((unsigned)time(0)); 
	log.printf("Random grid generation at every step. No Forces!\n The grid size is %d, 2d grid option neglected.",gsz);
      }


      parseFlag("GRIDFILE",gridfile);
      if(gridfile){
	//parseFlag("W2",w2instead);
	ingrid.open("grid.dat"); 
	log.printf("Space grid provided from a file, defined grid, xy, or widom options not used\n");
	if(!ingrid.is_open()){
	  log.printf("grid.dat not found, using defined grid.\n");
	}else{//read file
	  zdist = -1;
	  gsz=0;
	  string buff;

	  while(!ingrid.eof()){ //count lines
	    Vector v;
	    int ind,count;
	    double w1;
	    ingrid >> ind >> v[0] >> v[1] >> v[2] >> w1 >> count;
	    getline(ingrid,buff);
	    if(!ingrid.eof()){//fill arrays
	      grid_in.push_back(v);
	      grid_w.push_back(count/w1); //already divided by the box volume
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
      
      //log check
      log.printf("gridsize:\t%d\n",gsz);

      //SW potential  parameters (default MW model for water)
      
      
      
      A = 7.049556377;
      parse("ASW",A);
      B = 0.6022245584;
      parse("BSW",B);
      gamma = 1.2;
      parse("GAMMASW",gamma);
      costh0 = -1./3.;
      parse("CTHETASW",costh0);
      p=4;
      parse("PSW",p);
      q=0;
      parse("QSW",q);
      a = 1.8;
      parse("ACUTSW",a);
      epsilon=25.894776;
      parse("EPSSW",epsilon);
      sigma=0.23925;
      parse("SIGMASW",sigma);
      lambda=23.15;
      parse("LAMBDASW",lambda);
      r0=sigma/100.0;
      parse("R_0",r0);


      if(q != 0) {
	log.printf("Warning SW parameter q set automatically to q, q != 0 is not available for now.\n");
	q=0;
      }
      
      log.printf("SW parameters:\nA = %.8f\nB = %.8f\ngamma = %.2f\ncos(theta0) = %.8fp = %d\nq = 0\na = %.2f\nepsilon = %.8f kJ/mol\nsigma = %.8f nm\nlambda = %.4f\nr0 = %.4f\n",A,B,gamma,costh0,p,a,epsilon,sigma,lambda,r0);
      
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
      
      parseFlag("QUAD",quad);
      if(quad){
	log.printf("Quadratic soft-core\n");
      }



      parseFlag("NOPBCZ",nopbcz);
      if(nopbcz){
	log.printf("No pbcs along z\n");
      }

      //init Verlet lists
      
      r_cut=2*a*sigma; 
      parse("R_CUT",r_cut);
      r_skin=2.5*a*sigma; 
      parse("R_SKIN",r_skin);
      
      //log check
      log.printf("Verlet list:\tR_cut = %.4f\tR_skin = %.4f\n",r_cut,r_skin);
      log.flush();

      swlist.rcut=r_cut;  
      swlist.rskin=r_skin;  
      swlist.step=0;
      id=0;
      //create swlist.pos, array containing atom positions X
      Vector ze;
      ze.zero();

      //init neigh lists
      
      swlist.posit.resize(Natoms);
      fill(swlist.posit.begin(), swlist.posit.end(), ze);
      //create swlist.nn array containing the number of neighbors per gridpoint X
      swlist.nn.resize(gsz);
      fill(swlist.nn.begin(), swlist.nn.end(), 0);
      //create swlist.ni array containing the atom index of the neighbors X
      swlist.ni.resize(gsz);
      for (int i = 0; i < gsz; ++i){
	swlist.ni[i].resize(Natoms);
	fill(swlist.ni[i].begin(), swlist.ni[i].end(),0);
      }
      
      //kT [kJ/mol] (default = 2.49434 => 300K)
      
      kT=2.49434; 
      parse("KT",kT);
      
      //log check
      log.printf("kT = %.4f\tkJ/mol",kT);

      //general s exponential cutoff (default = 300.0)
      scoff=300.; 
      parse("SCOFF",scoff);
      log.printf("Small potential cutoff exponent: %.4f\n",scoff);


      //V exponential cutoff (default = 32.0)
      vcoff=32.; 
      parse("VCOFF",vcoff);
      log.printf("Small potential cutoff exponent: %.4f\n",vcoff);


      //log check
      
      
      //HISTOGRAM INPUT
      
      parseFlag("NOFIRST",nofirst);
      if(nofirst){
	log.printf("Nofirst option active, neglecting first frame in data collection.\n");
	swlist.step=-2;
      }

      //histogram s filter
      
      

      shist.resize(3);
      shist[0] = -kT*scoff;
      shist[1] = kT*scoff;
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
	dhist << "swlist.step\tNsmall\t<Nsmall>\tNlarge\t<Nlarge>\tNhist\t<Nhist>\tNtot\t" << endl;
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
      
      rank=comm.Get_rank(); //Rank of pr
      
      

      
      //stringstream A;
      //A<< "rdbg.dat."<< rank;
      //rdbg.open(A.str().c_str());

      //open debug file DBG
      //fdbg.open("dbg.dat");
      
      //fdbg.flush(); //DBG
    }

    
    //get the angle between two vectors a b
    double mumetaswhist::theta_ab(Vector av, Vector bv){
      double cdot = av[0]*bv[0]+av[1]*bv[1]+av[2]*bv[2];
      double amod = av.modulo();
      double bmod = bv.modulo();
      double C_ab = cdot/(amod*bmod);
      return(acos(C_ab));
    }    
    
    //get the angle between two vectors a b, plus the gradients dtheta/da dtheta/db
    double mumetaswhist::dtheta_ab(Vector av, Vector bv, Vector& dtheta_a, Vector& dtheta_b)
    {
      double cdot = av[0]*bv[0]+av[1]*bv[1]+av[2]*bv[2];
      double amod = av.modulo();
      double bmod = bv.modulo();
      double C_ab = cdot/(amod*bmod);
      double eps = 1.e-8;
      /*if(C_ab > 1.0-eps){
	C_ab = 1.0-eps;
      }else if(C_ab < -1.0+eps){
	C_ab=-1.0+eps;
	}*/
      double theta = acos(C_ab);
      if(fabs(C_ab*C_ab-1.) < eps){ //singularity, derivatives = 0
	dtheta_a[0] = dtheta_a[1] = dtheta_a[2] = 0.0;
	dtheta_b[0] = dtheta_b[1] = dtheta_b[2] = 0.0;
	return(theta);
      }
      double dtheta_dC = -(1./(sqrt(1.-C_ab*C_ab)));
      double den = 1./(amod*amod*bmod*bmod);
      double coeff=dtheta_dC*den;
      
      dtheta_a[0]=coeff * (bv[0]*amod*bmod-cdot*av[0]*bmod/amod);
      dtheta_a[1]=coeff * (bv[1]*amod*bmod-cdot*av[1]*bmod/amod);
      dtheta_a[2]=coeff * (bv[2]*amod*bmod-cdot*av[2]*bmod/amod);
      
      dtheta_b[0]=coeff * (av[0]*amod*bmod-cdot*bv[0]*amod/bmod);
      dtheta_b[1]=coeff * (av[1]*amod*bmod-cdot*bv[1]*amod/bmod);
      dtheta_b[2]=coeff * (av[2]*amod*bmod-cdot*bv[2]*amod/bmod);
      
      return(theta);
    }
    
    //get the angle penalty term between two vectors a b, plus the gradients dh/da dh/db
    double mumetaswhist::h_ab(Vector av, Vector bv, Vector& dh_a, Vector& dh_b)
    {

       //SW parameters
      //double A, B, gamma, lambda, epsilon, sigma, costh0, a;
      //int p,q;

      Vector dtheta_a, dtheta_b;
      double theta = dtheta_ab(av,bv,dtheta_a,dtheta_b);
      double amod = av.modulo();
      double bmod = bv.modulo();
      double rmax = a*sigma;
      double expexp = exp(gamma*sigma/(amod-rmax))*exp(gamma*sigma/(bmod-rmax));     
      double hab = lambda*epsilon*(cos(theta)-costh0)*(cos(theta)-costh0)*expexp;
      //for(int i=0; i<3; i++){
      //dh_a[i] = -lambda*epsilon*2.*(cos(theta)-costh0)*sin(theta)*expexp*dtheta_a[i]-hab*gamma*sigma/(amod-rmax)/(amod-rmax)*av[i]/amod;
      //dh_b[i] = -lambda*epsilon*2.*(cos(theta)-costh0)*sin(theta)*expexp*dtheta_b[i]-hab*gamma*sigma/(bmod-rmax)/(bmod-rmax)*bv[i]/bmod;
      //}
      
      dh_a = -lambda*epsilon*2.*(cos(theta)-costh0)*sin(theta)*expexp*dtheta_a-hab*gamma*sigma/(amod-rmax)/(amod-rmax)*av/amod;
      dh_b = -lambda*epsilon*2.*(cos(theta)-costh0)*sin(theta)*expexp*dtheta_b-hab*gamma*sigma/(bmod-rmax)/(bmod-rmax)*bv/bmod;
      
      //fdbg << "Inside" << setprecision(8) << dh_a[0] << "  " << dh_a[1]<< "  " << dh_a[2] << "  " << dh_b[0] << "  " << dh_b[1]<< "  " << dh_b[2] << endl; //DBG
      //fdbg.flush();
      return(hab);
    }
    
    //2 body potential
    double mumetaswhist::v2(double r, double rmin, double rmax, double Vmax)
    {
      
      double phi,sor,eps;
      eps = sigma/vcoff;
      if( r <= rmin){
	phi=Vmax;
      }else if(r > rmax-eps){
	phi=0.0;
      }else{
	sor=sigma/r;
	phi=A*epsilon*(B*pow(sor,1.*p)-1.)*exp(sigma/(r-rmax));
      }
      return(phi);
    }
    
    //2 body potential derivative
    double mumetaswhist::dv2(double r, double rmin, double rmax, double vr){
      double dphi,sor,eps;
      eps = sigma/vcoff;
      if( r <= rmin){
	dphi=0.0;
      }else if(r > rmax-eps){
	dphi=0.0;
      }else{
	sor=sigma/r;
	dphi=-A*epsilon*B*p*pow(sigma,1.*p)*pow(r,-1.*(p+1))*exp(sigma/(r-rmax))-vr*sigma/(r-rmax)/(r-rmax);
      }
      return(dphi);
    }


    //2 body potential (soft-core)
    double mumetaswhist::v2q(double r, double rmin, double rmax, double Vmax, double dv0){
      double phi,sor,eps;
      eps = sigma/vcoff;
      if( r <= rmin){
	phi=0.5*dv0*(r*r/rmin-rmin)+Vmax;
      }else if(r > rmax-eps){
	phi=0.0;
      }else{
	sor=sigma/r;
	phi=A*epsilon*(B*pow(sor,1.*p)-1.)*exp(sigma/(r-rmax));
      }
      return(phi);
    }
    
    //2 body potential derivative (soft-core)
    double mumetaswhist::dv2q(double r, double rmin, double rmax, double vr, double dv0){
      double dphi,sor,eps;
      eps = sigma/vcoff;
      if( r <= rmin){
	dphi=r*dv0/rmin;
      }else if(r > rmax-eps){
	dphi=0.0;
      }else{
	sor=sigma/r;
	dphi=-A*epsilon*B*p*pow(sigma,1.*p)*pow(r,-1.*(p+1))*exp(sigma/(r-rmax))-vr*sigma/(r-rmax)/(r-rmax);
      }
      return(dphi);
    }



    //3 body potential term + gradients, evaluate the 3 angular contributions
    double mumetaswhist::v3(Vector rA, Vector ri, Vector rj, double rmax, Vector& dv3_dri, Vector& dv3_drj){
      
      double theta;
      //h1 term, jAi
      Vector av = pbcDistance(rA,rj); //rjA
      double amod = av.modulo();
      Vector bv = pbcDistance(rA,ri); //riA
      double bmod = bv.modulo();
      Vector dh1_a,dh1_b;
      double eps = gamma*sigma/vcoff;
      double h1;
      if( amod > rmax-eps || bmod > rmax-eps){
	h1=0.0;
	dh1_a.zero();
	dh1_b.zero();
      }else{
	//fdbg << setprecision(12) << h1 << "  " << dh1_a[0] << "  " << dh1_a[1]<< "  " << dh1_a[2] << "  " << dh1_b[0] << "  " << dh1_b[1]<< "  " << dh1_b[2] << endl; //DBG
	h1 = h_ab(av,bv,dh1_a,dh1_b);
	//fdbg << setprecision(12) << h1 << "  " << dh1_a[0] << "  " << dh1_a[1]<< "  " << dh1_a[2] << "  " << dh1_b[0] << "  " << dh1_b[1]<< "  " << dh1_b[2] << endl; //DBG
	//fdbg.flush();
      }
 
      //h2 term, Aij
      av = pbcDistance(ri,rA); //rAi
      amod = av.modulo();
      bv = pbcDistance(ri,rj); //rji
      bmod = bv.modulo();
      Vector dh2_a,dh2_b;
      double h2;
      if( amod > rmax-eps || bmod > rmax-eps){
	h2=0.0;
	dh2_a.zero();
	dh2_b.zero();
      }else{
	h2 = h_ab(av,bv,dh2_a,dh2_b);
      }

      //h3 term, ijA
      av = pbcDistance(rj,ri); //rij
      amod = av.modulo();
      bv = pbcDistance(rj,rA); //rAj
      bmod = bv.modulo();
      Vector dh3_a,dh3_b;
      double h3;
      if( amod > rmax-eps || bmod > rmax-eps){
	h3=0.0;
	dh3_a.zero();
	dh3_b.zero();
      }else{
	h3 = h_ab(av,bv,dh3_a,dh3_b);
      }

      //evaluate derivatives

      //for(int i=0; i<3; i++){
      //dv3_dri[i] = dh1_b[i] - dh2_a[i] - dh2_b[i] + dh3_a[i];
      //dv3_drj[i] = dh1_a[i] + dh2_b[i] - dh3_a[i] - dh3_b[i];
      //}
      dv3_dri = dh1_b - dh2_a - dh2_b + dh3_a;
      dv3_drj = dh1_a + dh2_b - dh3_a - dh3_b;
      return(h1+h2+h3);
    }

    //LJ Potential
    double mumetaswhist::phi(double rq, double rmin, double rmax, double sig, double eps, double Vmax, double Vshift){
      double phi, sorq, sor6, sor12;
      if( rq <= rmin*rmin){
	phi=Vmax;
      }else if(rq >= rmax*rmax){
	phi=0.0;
      }else{
	sorq = sig*sig/rq;
	sor6 = sorq * sorq * sorq;
	sor12 = sor6 * sor6;
	phi=4*eps*(sor12-sor6)-Vshift;
      }
      return(phi);
    }

    
    //LJ derivative
    double mumetaswhist::dphi(double rq, double r, double rmin, double rmax, double sig, double eps){
      double dphi,sorq,sor6,sor12;
      if( rq <= rmin*rmin){
	dphi=0.0;
      }else if(rq >= rmax*rmax){
	dphi=0.0;
      }else{
	sorq = sig*sig/rq;
	sor6 = sorq * sorq * sorq;
	sor12 = sor6 * sor6;
	dphi=-24*eps*(2*sor12-sor6)/r;
      }
      return(dphi);
    }
    
    double mumetaswhist::phiq(double rq, double r, double rmin, double rmax, double sig, double eps, double Vmax, double Vshift, double dmin){
      double phi, sorq, sor6, sor12;
      if(rq <= rmin*rmin){
	phi=0.5*dmin*(rq/rmin-rmin)+Vmax;
      }else if(rq >= rmax*rmax){
	phi=0.0;
      }else{
	sorq = sig*sig/rq;
	sor6 = sorq * sorq * sorq;
	sor12 = sor6 * sor6;
	phi=4*eps*(sor12-sor6)-Vshift;
      }
      return(phi);
    }
    
    
    //LJ derivative
    double mumetaswhist::dphiq(double rq, double r, double rmin, double rmax, double sig, double eps, double dmin){
      double dphi, sorq, sor6, sor12;
      if( rq <= rmin*rmin){
	dphi=r*dmin/rmin;
      }else if(rq >= rmax*rmax){
	dphi=0.0;
      }else{
	sorq = sig*sig/rq;
	sor6 = sorq * sorq * sorq;
	sor12 = sor6 * sor6;
	dphi=-24*eps*(2*sor12-sor6)/r;
      }
      return(dphi);
    }
    
    

    // verlet list formation
    

    void mumetaswhist::sw_newlist(vector<Vector> &grid, swlist_s &swlist)
    {
      int i,j,gridsize,listsize;
      Vector rij, test;
      double mod_rij;
      
      gridsize=grid.size(); //get size of the grid
      listsize=at_list.size(); //get the size of the list

      //fill swlist.pos with new list positions

      for(j=0; j<listsize; ++j){
	swlist.posit[j] = getPosition(j);
      }
      
      unsigned int stride;
      unsigned int rank; 
      
      stride=comm.Get_size();  //Number of processes
      rank=comm.Get_rank(); //Rank of pr

      for(j=rank; j<gridsize; j+=stride) {     // sum over grid
	swlist.nn[j]=0; //reset nlistsize
	for(i=0; i<listsize; ++i) {                                           // sum over atoms
	  rij = pbcDistance(swlist.posit[i],grid[j]);
	  if(nopbcz) rij[2] = grid[j][2]-swlist.posit[i][2]; //nopbc z
	  mod_rij=rij.modulo();
	  if (mod_rij < swlist.rskin){ //if distance < rskin
	    swlist.ni[j][swlist.nn[j]]=i; //index in neighlist
	    swlist.nn[j]++; //increment nn 
	  }	  
	}
      }
    }

    void mumetaswhist::sw_checklist(vector<Vector> &grid, swlist_s &swlist)
    {
      unsigned int j; 
      double dr=(swlist.rskin-swlist.rcut)*0.5;
      Vector rij;
        
      for (j=0; j<at_list.size(); ++j) { 
	rij = pbcDistance(getPosition(j),swlist.posit[j]); //check position variations
	if(nopbcz) rij[2] = swlist.posit[j][2]-getPosition(j)[2]; //nopbc z
       	if( fabs(rij[0])>dr || fabs(rij[1])>dr || fabs(rij[2])>dr ) { 
	  id++; //rebuild list counter
	  sw_newlist(grid,swlist); 
	  break; 
	}
      }
    } 
    
    // calculator
    void mumetaswhist::calculate()    
    {
      int A,i,j,k;
      
      //fdbg<< swlist.step <<"\t start"<<endl;
      //fdbg.flush(); //DBG
      
      //Parallel parameters
      
      unsigned int stride;
      unsigned int rank; 
     
      stride=comm.Get_size();  //Number of processes
      rank=comm.Get_rank(); //Rank of pr

      //Box size
      Pbc pbc=getPbc();
      vector<double> LBC(3);
      for(i=0;i<3;++i) LBC[i]=pbc.getBox()[i][i];
      if(swlist.step==0) Z0=LBC[2];
      //if(swlist.step>744500) rdbg << "started check" << swlist.step << endl;

      double Vbox=getBox().determinant(); //box volume
      double dVA=1./gsz;  //volume element, in [Vbox] units
      if(gridfile) dVA=1./gsz2; //divided by box volume for imported grid
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
      
      //if(swlist.step>744500) rdbg << "grid ready " << gsz << endl;

      //fdbg<< swlist.step <<"\t grid built"<<endl;
      //fdbg.flush(); //DBG

      
      //init SW potential
      double Vmax;
      Vmax=v2(r0, 0.0, a*sigma, 0.0);
      //fdbg<< Vmax <<endl;
      
      double dv0;      
      if(quad) dv0=dv2(r0, 0.0, a*sigma, Vmax);
      
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
      //fdbg<< swlist.step << "\t" << zdist*zdist*LBC[2]*LBC[2] << "\t" << phiwall <<endl;
      
      //fdbg<< swlist.step <<"\t V initialized"<<endl;
      //fdbg.flush(); //DBG

      
      //create neigh lists  
      
      //create neigh lists  
      if(nofirst){
	if(swlist.step==-2){ //first step
	  sw_newlist(xgrid,swlist);
	}else{
	  if(widom){
	    sw_newlist(xgrid,swlist); //recompute lists
	  }else{
	    sw_checklist(xgrid,swlist); //check lists
	  }
	}
      }else{
	if(swlist.step==-1){ //first step
	  sw_newlist(xgrid,swlist);
	}else{
	  if(widom){
	    sw_newlist(xgrid,swlist); //recompute lists
	  }else{
	    sw_checklist(xgrid,swlist); //check lists
	  }
	}
      }

      ++swlist.step; //0 first step, or -1 if nofirst



      //fdbg<< swlist.step <<"\t V initialized"<<endl;
      //fdbg.flush(); //DBG

      //(swlist.step>744501) rdbg << "nl ready " << endl;
      double Su=0;
      double Du,modrij,modrAi,modrAj,modq,zu,fDu,dfDu, fDuw, MDu, ph;
      Vector rij,rAi,rAj;
      int atomi, atomj;
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
      
      for(int A=rank; A<gsz; A+=stride){

	//fdbg.flush();

	//fdbg<< swlist.step <<"\t grid loop "<< A << endl;
	//fdbg.flush(); //DBG
      
	//auxiliary array of vectors for derivatives
	vector<Vector> phiprime(swlist.nn[A]);
	fill(phiprime.begin(), phiprime.end(), ze);
	//auxiliary array of tesors for derivatives
	//vector<Tensor> phit(swlist.nn[A]);
	
	//z dependent wall 
	if(zdist == -1 && eps_Wall !=0.0) phiwall=phi(xgrid[A][2]*xgrid[A][2], rwallmax, r_cut_Wall, sigma_Wall, eps_Wall, wallmax, wallshift)+phi((LBC[2]-xgrid[A][2])*(LBC[2]-xgrid[A][2]), rwallmax, r_cut_Wall, sigma_Wall, eps_Wall, wallmax, wallshift);
	//if(swlist.step == 1){
	//  fdbg<< setprecision(8);
	//  fdbg << xgrid[i][2] << "\t" << phiwall <<endl;
	//}
	//auxiliary 
	Du=phiwall; //zero if no wall is present
	
	MDu=0;
	
	
	//xir counters to zero (grid point gdr)
	if(xir) fill(gdri.begin(), gdri.end(), 0);
	
	//cycle over neighbor particles
	
	for(i=0; i<swlist.nn[A]; ++i){
	  
	  //Skip further computations if threshold insertion energy is exceeded
	  if(beta*Du>scoff)   break;
	  
	  atomi = swlist.ni[A][i]; //i atomindex
	  //if(i==1){
	  //  fdbg<< atomid <<"  ";
	  //  fdbg.flush();
	  //}
	  //Delta U
	  rAi = pbcDistance(xgrid[A],getPosition(atomi));
	  modrAi=rAi.modulo();
	  
	  //Fill gdr histogram
	  if(gdrp[2]>0.0 && swlist.step>=0){
	    if(modrAi>=gdrp[0] && modrAi<gdrp[1]){
	      int bin= (int)((modrAi-gdrp[0])/bing);
	      gdr[bin]++;
	      Ngdr++;
	      if(xir) gdri[bin]++;
	    }
	    Ntotgdr++;
	  }



	  //2 BODY TERM

	  double ph,dph;
	  
	  if(quad){
	    //double mumetaswhist::v2q(double r, double rmin, double rmax, double Vmax, double dv0){
	    ph=v2q(modrAi, r0, a*sigma, Vmax, dv0);
	    
	    //double mumetaswhist::dv2q(double r, double rmin, double rmax, double Vmax, double vr, double dv0){
	    dph=dv2q(modrAi, r0, a*sigma, ph, dv0);

	    
	    
	  }else{
	    //double mumetaswhist::v2(double r, double rmin, double rmax, double Vmax)
	    ph=v2(modrAi, r0, a*sigma, Vmax);
	    
	    //double mumetaswhist::dv2(double r, double rmin, double rmax, double Vmax, double vr){
	    dph=dv2(modrAi, r0, a*sigma, ph);
	  }
	  
	  Du+=ph; //update potential difference
	  phiprime[i] = phiprime[i] + rAi*dph/modrAi; //potential derivative
	  //phit[i] = Tensor(phiprime[i],rAi); //Virial component

	  //3 BODY TERM
	  
	  //find triplets (j>i)
	  
	  for(j=i+1; j<swlist.nn[A]; ++j){
	    
	    atomj = swlist.ni[A][j]; //j atomindex
	    
	    Vector dv3i,dv3j;
	    
	    //compute 3 body for Aij triplet

	    //double mumetaswhist::v3(Vector rA, Vector ri, Vector rj, double rmax, Vector dv3_dri, Vector dv3_drj){
	    
	    double ph3 = v3(xgrid[A], getPosition(atomi), getPosition(atomj), a*sigma, dv3i, dv3j);
	    
	    //Increment Du with V3
	    
	    Du+=ph3; //update potential difference
	    
	    //Increment dDu/dr_i with 3 body contributions

	    phiprime[i] = phiprime[i] + dv3i; 
	    //Increment dDu/dr_j

	    phiprime[j] = phiprime[j] + dv3j; 
	    //fdbg << setprecision(8) << i << "\t" << j << "\t" << dv3i[0] << "\t" << dv3i[1]<< "\t" << dv3i[2] << "\t" << dv3j[0] << "\t" << dv3j[1] << "\t" << dv3j[2] << endl; //DBG
	    //fdbg.flush();
	    
	    }
	  
	  
	}
	//Fill insertion energy histogram
	if(swlist.step>=0){
	  //print insertion energy in array
	  if(dufile) bDuvec[A] = beta*Du;
	  
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
	
	//exponential (with cut off at coff)
	if(beta*Du>scoff){
	  fDu=0.0;
	  dfDu=0.0;
	}else{
	  fDu=exp(-beta*(Du));
	  fDuw=exp(-beta*(Du-shist[2]));
	  dfDu=-beta*fDu;
	}


	if(xir && swlist.step>=0){
	  for(int nb=0; nb<Nbing; nb++){
	    if(gdri[nb]>0){
	      xirhist[nb]+=fDuw*gdri[nb];
	    }
	  }
	  wcav+=fDuw;
	}
	
	Su+=fDu*dVA*grid_w[A]; //evaluate integral
	double dSoV=dfDu*dVA*grid_w[A];//Vbox;
	
	//if(swlist.step > 744501) rdbg << setprecision(8) << A << "\t" << scientific << Su << endl; //DBG
	//derivatives -- cycle again over neighbors
	
	//fdbg<< swlist.step <<"\t Re-loop over i for derivatives." << endl;
	//fdbg.flush(); //DBG
	if(!widom){
	  for(i=0; i < swlist.nn[A]; ++i){
	    atomi = swlist.ni[A][i]; //atomindex
	    deriv[atomi]+=dSoV*phiprime[i];
	    rAi = pbcDistance(xgrid[A],getPosition(atomi));
	    virial-=dSoV*Tensor(phiprime[i],rAi);
	  }
	}
	//if(swlist.step>744501)  rdbg << setprecision(8) << A << "\tDerivatives ready\t" << scientific << virial[0][0] << endl; //DBG
	vector<Vector>().swap(phiprime);
	//vector<Tensor>().swap(phit);
      }
      
      
      //fdbg<< swlist.step <<"\t Loop closed, Du = "<< Du << endl;
      //fdbg.flush(); //DBG
      

      comm.Sum(Su);
      comm.Sum(deriv);
      comm.Sum(virial);
      
      //fdbg<< swlist.step <<"\t Communication done" << endl;
      //fdbg.flush(); //DBG

      if(dufile && swlist.step>=0){
	comm.Sum(bDuvec);
	if(rank == 0){
	  for(int j=0 ; j<gsz; j++){ 
	    fdu << setprecision(8);
	    fdu << (swlist.step)*gsz+j << "\t" << bDuvec[j] << "\t" << swlist.step << endl;
	  }
	}
      }


      double S;
      if(Su==0){
	S=kT*scoff; //When the integral is 0 mu=coff* kT
	//if(swlist.step>744501) rdbg << "s value = " << kT*scoff << endl;
      }else{
	S=-kT*(std::log(Su));
	//if(swlist.step>744501) rdbg << "s value = " << -kT*(std::log(Su)) << endl;
      }
      setValue(S);
      for(j=0; j < Natoms; ++j)	{
	setAtomsDerivatives(j, -kT/Su*deriv[j]);
      }
      
      //Add total volume derivative
      
      setBoxDerivatives(-kT/Su*virial);//+SuT);

      //fdbg<< swlist.step <<"\t Everything set" << endl;
      //fdbg.flush(); //DBG
      vector<Vector>().swap(xgrid);
      //if(swlist.step>744501)  rdbg << "All set, next cycle" << endl; //DBG

     //WEIGHTED DISTRIBUTION
      //read configuration weight
      if(swlist.step>=0){
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
	  if((swlist.step+1)%histprint == 0 && rank == 0){
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
	    dhist << swlist.step << "\t" << Nsmall << "\t"  << 1.0*Nsmall_av/(swlist.step+1) << "\t"<< Nlarge << "\t"  << 1.0*Nlarge_av/(swlist.step+1)<< "\t" << Nhist << "\t"  << 1.0*Nhist_av/(swlist.step+1) << "\t" << Ntot << "\t" << fixed << setprecision(8) << weight << endl; //hist N data
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
	  if((swlist.step+1)%histprint == 0 && rank == 0){
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
