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
   
    class mumetasw : public Colvar {
      //      SwitchingFunction switchingFunction; //instance sw.f
      int gsz,func, Natoms, id, gsz2;
      double l_func, sigma_LJ, eps_LJ, R_min, R_max, u,  r_cut, r_skin, vcoff, scoff, kT, r0, zdist, eps_Wall, sigma_Wall, r_cut_Wall, wallv;
      bool vshift, issig, quad, nopbcz,gridfile,w2instead;
      vector<double> LJpar;
      vector<double> Wallpar;
      vector<double> LJr;
      vector<Vector> grid_in;
      vector<double> grid_w;
      vector<int> griddim;
      vector<AtomNumber> at_list;
      double Z0;
      
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

      
      mumetasw(const ActionOptions&);
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

    PLUMED_REGISTER_ACTION(mumetasw,"MUMETASW")

    void mumetasw::registerKeywords(Keywords& keys){
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
      keys.add("optional","XY","xy plane grid");
      keys.addFlag("NOPBCZ",false,"remove pbc along z");
      keys.add("optional","WALLJ","Wall potential parameters (for LJ walls)");
      keys.add("optional","WALLVAL","Wall potential value");
      keys.add("optional","KT","kT, default 2.49434 kJ/mol =>  300 K");
      keys.add("optional","SCOFF","general exponential cutoff (default=300)");
      keys.add("optional","VCOFF","potential exponentials cutoff (default=32)");
      keys.remove("NOPBC");
    }
    
    mumetasw::mumetasw(const ActionOptions&ao):
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
      
      parseFlag("GRIDFILE",gridfile);
      if(gridfile){
	//parseFlag("W2",w2instead);
	ingrid.open("grid.dat"); 
	log.printf("Space grid provided from a file, defined grid, xy not used\n");
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
	log.printf("Warning SW parameter q set automatically to 0, q != 0 is not available for now.\n");
	q=0;
      }
      
      log.printf("SW parameters:\nA = %.8f\nB = %.8f\ngamma = %.2f\ncos(theta0) = %.8f\np = %d\nq = 0\na = %.2f\nepsilon = %.8f kJ/mol\nsigma = %.8f nm\nlambda = %.4f\nr0 = %.4f\n",A,B,gamma,costh0,p,a,epsilon,sigma,lambda,r0);
      
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
      log.printf("kT = %.4f\tkJ/mol\n",kT);

      //general s exponential cutoff (default = 300.0)
      scoff=300.; 
      parse("SCOFF",scoff);
      log.printf("S cutoff exponent: %.4f\n",scoff);


      //V exponential cutoff (default = 32.0)
      vcoff=32.; 
      parse("VCOFF",vcoff);
      log.printf("Small potential cutoff exponent: %.4f\n",vcoff);
      //log check
      
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
    double mumetasw::theta_ab(Vector av, Vector bv){
      double cdot = av[0]*bv[0]+av[1]*bv[1]+av[2]*bv[2];
      double amod = av.modulo();
      double bmod = bv.modulo();
      double C_ab = cdot/(amod*bmod);
      return(acos(C_ab));
    }    
    
    //get the angle between two vectors a b, plus the gradients dtheta/da dtheta/db
    double mumetasw::dtheta_ab(Vector av, Vector bv, Vector& dtheta_a, Vector& dtheta_b)
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
    double mumetasw::h_ab(Vector av, Vector bv, Vector& dh_a, Vector& dh_b)
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
    double mumetasw::v2(double r, double rmin, double rmax, double Vmax)
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
    double mumetasw::dv2(double r, double rmin, double rmax, double vr){
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
    double mumetasw::v2q(double r, double rmin, double rmax, double Vmax, double dv0){
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
    double mumetasw::dv2q(double r, double rmin, double rmax, double vr, double dv0){
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
    double mumetasw::v3(Vector rA, Vector ri, Vector rj, double rmax, Vector& dv3_dri, Vector& dv3_drj){
      
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
    double mumetasw::phi(double rq, double rmin, double rmax, double sig, double eps, double Vmax, double Vshift){
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
    double mumetasw::dphi(double rq, double r, double rmin, double rmax, double sig, double eps){
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
    
    double mumetasw::phiq(double rq, double r, double rmin, double rmax, double sig, double eps, double Vmax, double Vshift, double dmin){
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
    double mumetasw::dphiq(double rq, double r, double rmin, double rmax, double sig, double eps, double dmin){
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
    

    void mumetasw::sw_newlist(vector<Vector> &grid, swlist_s &swlist)
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

    void mumetasw::sw_checklist(vector<Vector> &grid, swlist_s &swlist)
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
    void mumetasw::calculate()    
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
      
      if(swlist.step==0) sw_newlist(xgrid,swlist);
      ++swlist.step; 
      sw_checklist(xgrid,swlist);
      
      //fdbg<< swlist.step <<"\t V initialized"<<endl;
      //fdbg.flush(); //DBG
      
      //(swlist.step>744501) rdbg << "nl ready " << endl;
      double Su=0;
      double Du,modrij,modrAi,modrAj,modq,zu,fDu,dfDu, MDu, ph;
      Vector rij,rAi,rAj;
      int atomi, atomj;
      double beta=1./kT;
      
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
	
	//cycle over neighbor particles
	
	for(i=0; i<swlist.nn[A]; ++i){
	  
	  //Skip further computations if threshold insertion energy is exceeded
	  if(beta*Du>scoff)   break;
	  
	  atomi = swlist.ni[A][i]; //i atomindex
	  //if(i==1){
	  //fdbg << atomi <<"  " << i << endl;
	  fdbg.flush();
	  //}
	  //Delta U
	  rAi = pbcDistance(xgrid[A],getPosition(atomi));
	  modrAi=rAi.modulo();
	  
	  //2 BODY TERM

	  double ph,dph;
	  
	  if(quad){
	    //double mumetasw::v2q(double r, double rmin, double rmax, double Vmax, double dv0){
	    ph=v2q(modrAi, r0, a*sigma, Vmax, dv0);
	    
	    //double mumetasw::dv2q(double r, double rmin, double rmax, double Vmax, double vr, double dv0){
	    dph=dv2q(modrAi, r0, a*sigma, ph, dv0);

	    
	    
	  }else{
	    //double mumetasw::v2(double r, double rmin, double rmax, double Vmax)
	    ph=v2(modrAi, r0, a*sigma, Vmax);
	    
	    //double mumetasw::dv2(double r, double rmin, double rmax, double Vmax, double vr){
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

	    //double mumetasw::v3(Vector rA, Vector ri, Vector rj, double rmax, Vector dv3_dri, Vector dv3_drj){
	    
	    double ph3 = v3(xgrid[A], getPosition(atomi), getPosition(atomj), a*sigma, dv3i, dv3j);
	    
	    //Increment Du with V3
	    
	    Du+=ph3; //update potential difference
	    
	    //Increment dDu/dr_i with 3 body contributions
	    //fdbg << atomj <<"  " << j << "\t\t";
	    //fdbg.flush();
	    phiprime[i] = phiprime[i] + dv3i; 
	    //Increment dDu/dr_j

	    phiprime[j] = phiprime[j] + dv3j; 
	    //fdbg << setprecision(8) << i << "\t" << j << "\t" << dv3i[0] << "\t" << dv3i[1]<< "\t" << dv3i[2] << "\t" << dv3j[0] << "\t" << dv3j[1] << "\t" << dv3j[2] << endl; //DBG
	    //fdbg.flush();
	    
	    }
	  
	  
	}
	
	//exponential (with cut off at coff)
	if(beta*Du>scoff){
	  fDu=0.0;
	  dfDu=0.0;
	}else{
	  fDu=exp(-beta*(Du));
	  dfDu=-beta*fDu;
	}
	
	Su+=fDu*dVA*grid_w[A]; //evaluate integral
	double dSoV=dfDu*dVA*grid_w[A];//Vbox;
	
	//if(swlist.step > 744501) rdbg << setprecision(8) << A << "\t" << scientific << Su << endl; //DBG
	//derivatives -- cycle again over neighbors
	
	//fdbg<< swlist.step <<"\t Re-loop over i for derivatives." << endl;
	//fdbg.flush(); //DBG
	for(i=0; i < swlist.nn[A]; ++i){
	  atomi = swlist.ni[A][i]; //atomindex
	  deriv[atomi]+=dSoV*phiprime[i];
	  rAi = pbcDistance(xgrid[A],getPosition(atomi));
	  virial-=dSoV*Tensor(phiprime[i],rAi);
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
      
      if(Su==0){
	setValue(kT*scoff); //When the integral is 0 mu=coff* kT
	//if(swlist.step>744501) rdbg << "s value = " << kT*scoff << endl;
      }else{
	setValue(-kT*(std::log(Su)));
	//if(swlist.step>744501) rdbg << "s value = " << -kT*(std::log(Su)) << endl;
      }
      for(j=0; j < Natoms; ++j)	{
	setAtomsDerivatives(j, -kT/Su*deriv[j]);
      }
      
      //Add total volume derivative
      
      setBoxDerivatives(-kT/Su*virial);//+SuT);

      //fdbg<< swlist.step <<"\t Everything set" << endl;
      //fdbg.flush(); //DBG
      vector<Vector>().swap(xgrid);
      //if(swlist.step>744501)  rdbg << "All set, next cycle" << endl; //DBG

    }
  }
}
