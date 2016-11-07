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
   
    class mudelsw : public Colvar {
      //      SwitchingFunction switchingFunction; //instance sw.f
      int Natoms, id, Ndel, Nr;
      double R_min, R_max, r_cut, r_skin, scoff, vcoff, r0, Ubin, Wmax, eps_Wall, Ztot, kT, Du0;
      bool vshift, issig, shuff, rand, quad, nopbcz;
      vector<double> Swr;
      vector<AtomNumber> at_list;
      vector<int> allind;
      vector<int> delind;
      
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
	vector<Vector> pos;
      } swlist;

    public:

      
      mudelsw(const ActionOptions&);
      virtual void calculate();
      string to_string(int s);
      double theta_ab(Vector av, Vector bv);
      double dtheta_ab(Vector av, Vector bv, Vector& dtheta_a, Vector& dtheta_b);
      double h_ab(Vector av, Vector bv, Vector& dh_a, Vector& dh_b);
      double v2(double r, double rmin, double rmax, double Vmax);
      double dv2(double r, double rmin, double rmax, double Vmax, double vr);
      double v2q(double r, double rmin, double rmax, double Vmax, double dv0);
      double dv2q(double r, double rmin, double rmax, double Vmax, double vr, double dv0);
      double v3(Vector ri, Vector rj, Vector rk, double rmax, Vector& dv3_dri, Vector& dv3_drj, Vector& dv3_drk);
      double phi(double rq, double rmin, double rmax, double sigma, double eps, double Vmax, double Vshift);

      double phiq(double rq, double r, double rmin, double rmax, double sigma, double eps, double Vmax, double Vshift, double dmin);
      double dphi(double rq, double r, double rmin, double rmax, double sigma, double eps);
      double dphiq(double rq, double r, double rmin, double rmax, double sigma, double eps, double dmin);
      void sw_newlist(vector<Vector> &del, swlist_s &swlist);
      void sw_checklist(vector<Vector> &del, swlist_s &swlist);
      static void registerKeywords(Keywords& keys);
      //debugfile DBG
      ofstream fdbg;
      ofstream rdbg;
    };

    PLUMED_REGISTER_ACTION(mudelsw,"MUDELSW") 

    void mudelsw::registerKeywords(Keywords& keys){
      Colvar::registerKeywords(keys);
      keys.add("atoms","GROUP","the group of atoms we are calculating the potential energy for");
      keys.add("optional","NDEL","Number of deletions");
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
      keys.add("optional","HISTPRINT","how often print the histogram (steps)");
      keys.add("optional","R_0","2-body term potential inner cut-off");
      keys.add("optional","R_SKIN","Verlet lists, default r_skin = 2.5 a*sigma");
      keys.add("optional","R_CUT","Verlet lists, default r_cut = 2*a*sigma");
      keys.add("optional","SCOFF","general exponential cutoff (default=300)");
      keys.add("optional","VCOFF","potential exponentials cutoff (default=32)");
      keys.add("optional","DUZERO","Du0 value, to set the scale of the evaluated exponentials.");
      keys.addFlag("SHUFFLE",false,"Set to TRUE to get shuffle atom index for deletion");
      keys.addFlag("RAND",false,"Set to TRUE to compute the deletion of different atoms at any step");
      keys.addFlag("QUAD",false,"Set to TRUE if for the quadratic soft-core");
      keys.add("optional","KT","kT, default 0.73167 kJ/mol =>  88 K");
      keys.remove("NOPBC");
    }
    
    mudelsw::mudelsw(const ActionOptions&ao):
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
      
      //SW potential  parameters
      
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

      log.printf("SW parameters:\nA = %.8f\nB = %.8f\ngamma = %.2f\ncos(theta0) = %.8f\np = %d\nq = 0\na = %.2f\nepsilon = %.8f kcal/mol\nsigma = %.8f A\nlambda = %.4f\nr0 = %.4f\n",A,B,gamma,costh0,p,a,epsilon,sigma,lambda,r0);
              
      
      parseFlag("QUAD",quad);
      if(quad){
	log.printf("Quadratic soft-core\n");
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
      //vector<Vector> posvec(AtomNumber);
      Vector ze;
      ze.zero();
      
      swlist.pos.resize(Natoms);
      fill(swlist.pos.begin(), swlist.pos.end(), ze);
      
      swlist.nn.resize(Ndel);
      fill(swlist.nn.begin(), swlist.nn.end(), 0);
      //create swlist.ni array containing the atom index of the neighbors X
      //vector<vector<double> > neighi;
      swlist.ni.resize(Ndel);
      for (int i = 0; i < Ndel; ++i){
	swlist.ni[i].resize(Natoms);
	fill(swlist.ni[i].begin(), swlist.ni[i].end(),0);
      }
      
      //log check
      log.printf("pos size:\t  %d\n",swlist.pos.size());
      log.printf("neigh number size:\t  %d\n",swlist.nn.size());
      log.printf("1st list size:\t  %d\n",swlist.ni[0].size());
      log.printf("Last list size:\t  %d\n",swlist.ni[swlist.nn.size()-1].size());
      log.flush();
      

      //kT [kJ/mol] (default = 2.49434 => 300K)
      
      kT=2.49434; 
      parse("KT",kT);
      
      //log check
      log.printf("kT = %.4f\tkJ/mol\n",kT);

      Du0=30.0; 
      parse("DUZERO",Du0);
      log.printf("DU_0 = %.4f\n",Du0);
      
      //general s exponential cutoff (default = 300.0)
      scoff=300.; 
      parse("SCOFF",scoff);
      log.printf("S cutoff exponent: %.4f\n",scoff);
      
      
      //V exponential cutoff (default = 32.0)
      vcoff=32.; 
      parse("VCOFF",vcoff);
      log.printf("Small potential cutoff exponent: %.4f\n",vcoff);
      
      
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
      
      //      stringstream A;
      //A<< "rdbg.dat."<< rank;
      //rdbg.open(A.str().c_str());
      //fdbg << "DEBUG FILE \n";  //DBG
      //fdbg.open("dbg.dat");
      //fdbg.flush(); //DBG
      //log.printf("dbg.dat open  \n");
      
    }
    
    
    
    //get the angle between two vectors a b
    double mudelsw::theta_ab(Vector av, Vector bv){
      double cdot = av[0]*bv[0]+av[1]*bv[1]+av[2]*bv[2];
      double amod = av.modulo();
      double bmod = bv.modulo();
      double C_ab = cdot/(amod*bmod);
      return(acos(C_ab));
    }    
    
    //get the angle between two vectors a b, plus the gradients dtheta/da dtheta/db
    double mudelsw::dtheta_ab(Vector av, Vector bv, Vector& dtheta_a, Vector& dtheta_b)
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
    double mudelsw::h_ab(Vector av, Vector bv, Vector& dh_a, Vector& dh_b)
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
    double mudelsw::v2(double r, double rmin, double rmax, double Vmax)
    {
      
      double phi,sor,eps;
      eps = sigma/vcoff;
      if( r <= rmin){
	phi=Vmax;
      }else if(r > rmax-eps){
	phi=0.0;
      }else{
	sor=sigma/r;
	phi=A*epsilon*(B*pow(sor,1.0*p)-1.)*exp(sigma/(r-rmax));
      }
      return(phi);
    }
    
    //2 body potential derivative
    double mudelsw::dv2(double r, double rmin, double rmax, double Vmax, double vr){
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
    double mudelsw::v2q(double r, double rmin, double rmax, double Vmax, double dv0){
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
    double mudelsw::dv2q(double r, double rmin, double rmax, double Vmax, double vr, double dv0){
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
    double mudelsw::v3(Vector ri, Vector rj, Vector rk, double rmax, Vector& dv3_dri, Vector& dv3_drj, Vector& dv3_drk){
      
      //h1 term, kij
      Vector av = pbcDistance(ri,rk); //rki
      double amod = av.modulo();
      Vector bv = pbcDistance(ri,rj); //rji
      double bmod = bv.modulo();
      Vector dh1_a,dh1_b;
      double eps = gamma*sigma/vcoff;
      double h1;
      if( amod > rmax-eps || bmod > rmax-eps){
	h1=0.0;
	dh1_a.zero();
	dh1_b.zero();
      }else{
	h1 = h_ab(av,bv,dh1_a,dh1_b);
      }

      //h2 term, ijk
      av = pbcDistance(rj,ri); //rij
      amod = av.modulo();
      bv = pbcDistance(rj,rk); //rkj
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

      //h3 term, jki
      av = pbcDistance(rk,rj); //rjk
      amod = av.modulo();
      bv = pbcDistance(rk,ri); //rik
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

      dv3_dri = - dh1_a - dh1_b + dh2_a + dh3_b;
      dv3_drj = + dh1_b - dh2_a - dh2_b + dh3_a;
      dv3_drk = + dh1_a + dh2_b - dh3_a - dh3_b;
      
      return(h1+h2+h3);
    }

    

    //LJ Potential
    double mudelsw::phi(double rq, double rmin, double rmax, double sigma, double eps, double Vmax, double Vshift){
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
    double mudelsw::dphi(double rq, double r, double rmin, double rmax, double sigma, double eps){
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

     
    double mudelsw::phiq(double rq, double r, double rmin, double rmax, double sigma, double eps, double Vmax, double Vshift, double dmin){
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
    double mudelsw::dphiq(double rq, double r, double rmin, double rmax, double sigma, double eps, double dmin){
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


    void mudelsw::sw_newlist(vector<Vector> &del, swlist_s &swlist)
    {
      int i,j,ndel,listsize;
      Vector rij, test;
      double mod_rij;
  
      for(j=0; j<Natoms; ++j){
	swlist.pos[j] = getPosition(j);
      }
      
      unsigned int stride;
      unsigned int rank; 
      
      stride=comm.Get_size();  //Number of processes
      rank=comm.Get_rank(); //Rank of pr
      
      for(j=rank; j<Ndel; j+=stride) {     // sum over deleted atoms
	swlist.nn[j]=0; //reset nlistsize
	for(i=0; i<Natoms; ++i) {                                           // sum over atoms
	  rij = pbcDistance(swlist.pos[i],del[j]);
	  mod_rij=rij.modulo();
	  if (mod_rij < swlist.rskin && i!=delind[j]){ //if distance < rskin and different atom
	    swlist.ni[j][swlist.nn[j]]=i; //index in neighlist
	    swlist.nn[j]++; //increment nn 
	  }	  
	}
      }
    }
    
    void mudelsw::sw_checklist(vector<Vector> &del, swlist_s &swlist)
    {
      unsigned int j; 
      double dr=(swlist.rskin-swlist.rcut)*0.5;
      Vector rij;
      
      for (j=0; j<Natoms; ++j) { 
	rij = pbcDistance(getPosition(j),swlist.pos[j]); //check position variations
       	if( fabs(rij[0])>dr || fabs(rij[1])>dr || fabs(rij[2])>dr ) { 
	  id++; //rebuild list counter
	  sw_newlist(del,swlist); 
	  break; 
	}
      }
    } 
    
    
    
    
    
    // calculator
    void mudelsw::calculate()    
    {
      //fdbg<< swlist2.step <<"\t start"<<endl;
      //fdbg.flush(); //DBG
      int i,j,k;
      //Tensor virial;
      
      //Parallel parameters
      
      unsigned int stride;
      unsigned int rank; 
      
      stride=comm.Get_size();  //Number of processes
      rank=comm.Get_rank(); //Rank of pr
      
      double ledN = 1./Ndel;
      //Box size
      
      vector<double> LBC(3);
      Pbc pbc=getPbc();
      for(i=0;i<3;++i) LBC[i]=pbc.getBox()[i][i];

      double Vbox=getBox().determinant(); //box volume
      
      
      
      int index;
      
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
   
      //init SW potential
      double Vmax;
      Vmax=v2(r0, 0.0, a*sigma, 0.0);
      //fdbg << "Vmax = "<< Vmax << endl;
      double dv0;      
      if(quad) dv0=dv2(r0, 0.0, a*sigma, Vmax, Vmax);
      
      if(rand){
	std:: random_shuffle (allind.begin(), allind.end()); //reshuffle indexes
	for(i=0; i<Ndel; i++) delind[i]=allind[i];  // rechose deletion indexes
      }
      
      for(i=0; i<Ndel; i++){ //find deletion positions
	delvec[i]=getPosition(delind[i]);
      }
      
      
      //create neigh lists 
      if(swlist.step==0){ //first step, create
	sw_newlist(delvec,swlist);
      }else{ // further steps 
	if(rand){ // new deletion list, recreate neigh lists
	  sw_newlist(delvec,swlist);
	}else{ // might be the same, check
	  sw_checklist(delvec,swlist);
	}
      }
      ++swlist.step;
      
      //init CV
      
      double Su=0.0;
      double Du,modrij,modrik,modkj,modq,zu,fDu,dfDu, MDu;
      Vector rij,rik,rkj;
      int atomj,atomk;
      double beta=1./kT;
      Vector dv3i,dv3j,dv3k;
      
      for(i=rank; i<Ndel; i+=stride){

	vector<Vector> phiprime(swlist.nn[i]);
	Vector phiprime_i;
	phiprime_i.zero(); //derivative on deleted molecule
	//cycle over neighbor particles
	
	Du=0;
	MDu=0;


	for(j=0; j<swlist.nn[i]; ++j){ //loop on all ij-couples 
	  //nlist.ni[i][j] i: grid point, j: neighbor molecule
	  atomj = swlist.ni[i][j]; //atomindex
	  Vector rj = getPosition(atomj);
	  //if(swlist.step == 1){
	  //  rdbg << "\t\t" << atomj << endl;
	  //  rdbg.flush();
	  //}
	  //fdbg << atomj << "  " << j << endl;
	  //fdbg.flush();

	  rij = pbcDistance(delvec[i],rj);
	  modrij=rij.modulo();
	  
	  //2 BODY TERM
	  double ph,dph;
	  if(quad){
	    //double mumetasw::v2q(double r, double rmin, double rmax, double Vmax, double dv0){
	    ph=v2q(modrij, r0, a*sigma, Vmax, dv0);
	    
	    //double mumetasw::dv2q(double r, double rmin, double rmax, double Vmax, double vr, double dv0){
	    dph=dv2q(modrij, r0, a*sigma, Vmax, ph, dv0);

	  }else{
	    //double mumetasw::v2(double r, double rmin, double rmax, double Vmax)
	    ph=v2(modrij, r0, a*sigma, Vmax);
	    
	    //double mumetasw::dv2(double r, double rmin, double rmax, double Vmax, double vr){
	    dph=dv2(modrij, r0, a*sigma, Vmax, ph);
	  }
	  
	  //if(swlist2.step == 1) fdbg << j << "\t" << atomj << "\t" << swlist2.nn[i] << "\t" << modrij << "\t" << ph <<  endl;
	  
	  Du+=ph;
	  phiprime_i  = phiprime_i - rij*dph/modrij; //v2 derivative on deleted part
	  phiprime[j] = phiprime[j] + rij*dph/modrij; //v2 derivative other part

	  //3 BODY TERM
	  
	  //find triplets (k>j)
	  
	  for(k=j+1; k<swlist.nn[i]; ++k){ 
	    
	    atomk = swlist.ni[i][k]; //k atomindex
	    
	    //compute 3 body for ijk triplet
	    
	    //double mumetasw::v3(Vector rA, Vector ri, Vector rj, double rmax, Vector dv3_dri, Vector dv3_drj){
	    
	    double ph3 = v3(delvec[i], rj, getPosition(atomk), a*sigma, dv3i, dv3j, dv3k);
	    
	    //Increment Du with V3
	    
	    Du+=ph3; //update potential difference

	    //Increment dDu/dr_i with 3 body contributions

	    //fdbg << atomk <<"  " << k << "\t\t";
	    //fdbg.flush();
	  
	    
	    phiprime_i = phiprime_i + dv3i; 
	    
	    //Increment dDu/dr_j
	    
	    phiprime[j] = phiprime[j] + dv3j; 
	    
	    //Increment dDu/dr_k
	    
	    phiprime[k] = phiprime[k] + dv3k; 
	    
	  }
	  //fdbg << swlist.step << "\t" << delind[i] << "\t" << atomj << "\t" << setprecision(8) << scientific << Du << "\t" <<endl;
	       
	}
	if(beta*(Du-Du0)<-scoff){
	  fDu=0.0;
	  dfDu=0.0;
	}else{
	  fDu=exp(beta*(Du-Du0));
	  dfDu=beta*fDu;
	}
	
	Su+=fDu*ledN; //evaluate integral
	double dSoV=dfDu*ledN;//Vbox;
	
	//derivatives -- cycle again over neighbors
	deriv[delind[i]]+=dSoV*phiprime_i; //contribution on deleted particle
	
	for(j=0; j < swlist.nn[i]; ++j){
	  atomj = swlist.ni[i][j]; //atomindex
	  deriv[atomj]+=dSoV*phiprime[j];
	  rij = pbcDistance(delvec[i],getPosition(atomj));
	  virial-=dSoV*Tensor(phiprime[j],rij); //Virial component
	  
	}
	vector<Vector>().swap(phiprime);
	  
	fdbg << swlist.step << "\t" << delind[i]+1 << "\t" << setprecision(8) << scientific << Du << "\t" << fDu << "\t" << Su <<endl;

      }

    
      comm.Sum(Su);
      comm.Sum(deriv);
      comm.Sum(virial);


      if(Su==0){
	setValue(-kT*scoff+Du0); //When the integral is 0 mu=coff* kT (very unlikely!!!)
      }else{
	setValue(Du0+kT*(std::log(Su)));
      }
      for(j=0; j < Natoms; ++j)	{
	setAtomsDerivatives(j, kT/Su*deriv[j]);
      }
      
      for(j=0; j < Natoms; ++j)	setAtomsDerivatives(j, deriv[j]);
      setBoxDerivatives(kT/Su*virial);
      
      
    }
  }
}

