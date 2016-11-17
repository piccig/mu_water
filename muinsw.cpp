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
      int gsz,func, Natoms, id, gsz2;
      double l_func, u,  coff, kT, zdist, eps_Wall, sigma_Wall, r_cut_Wall, wallv;
      double qOO, qHH, qOH, C0, sigmaOO, sigmaHH, sigmaOH, epsilonOO, epsilonHH, epsilonOH, AOO, AHH, AOH, BOO, BHH, BOH;  //Potential parameters
      double R_min_OO, R_min_HH, R_min_OH, R_max, r_cut, r_skin, r0;                                                       //Structural parameters
      double V_minOO, V_minHH, V_minOH, dVmin_OO, dV_minHH, dV_minOH;                                                      //Potential parameters
      double rtest, rstep, test_ph, test_ph_qs, test_dph, test_dph_qs                                                      //Potential testing tools
      double Z0;
      bool vshift, issig, quad, nopbcz,gridfile;
      vector<double> TIP3Pcharge;
      vector<double> TIP3Psigma;
      vector<double> TIP3Pepsilon;
      vector<double> Wallpar;
      vector<double> LJr;
      vector<Vector> grid_in;
      vector<double> grid_w;
      vector<int> griddim;
      vector<AtomNumber> at_list;
      
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
      double phis(double rq, double r, double rmin, double rmax, double sigma, double eps, double Vmax, double Vshift, double r0, bool issig, double coff);
      double phiq(double rq, double r, double rmin, double rmax, double sigma, double eps, double Vmax, double Vshift, double dmin);
      double phi_tip3p_qs(double r, double rmin, double rmax, double qq, double Apar, double Bpar, double coff, double Vmin, double dmin){
      double dphi_tip3p(double r, double qq, double Apar, double Bpar){
      double dphi(double rq, double r, double rmin, double rmax, double sigma, double eps);
      double dphis(double rq, double r, double rmin, double rmax, double sigma, double eps, double Vmax, double Vshift, double r0, bool issig, double coff);
      double dphiq(double rq, double r, double rmin, double rmax, double sigma, double eps, double dmin);
      double dphi_tip3p(double r, double qq, double Apar, double Bpar){
      double dphi_tip3p_qs(double r, double rmin, double rmax, double qq, double Apar, double Bpar, double coff, double Vmax, double dmin){
      void du_newlist(vector<AtomNumber> &list, vector<Vector> &grid, dulist_s &dulist);
      void du_checklist(vector<AtomNumber> &list, vector<Vector> &grid, dulist_s &dulist);
      
      static void registerKeywords(Keywords& keys);
      //debugfile DBG
      ofstream fdbg;
      ofstream rdbg;
      ifstream ingrid;
    };

    PLUMED_REGISTER_ACTION(muinsw,"MUINSW")

    void muinsw::registerKeywords(Keywords& keys){
      Colvar::registerKeywords(keys);
      keys.add("atoms","GROUP","the group of atoms we are calculating the dipole moment for");
      keys.add("optional","CHARGES","TIP3P charges (default:   qO = -0.834 C\t qH = 0.417 C)");
      keys.add("optional","LJSIGMA","TIP3P LJ sigma parameters (default:   sigmaOO = 0.31507 nm\t sigmaHH = 0.04000 nm\t sigmaOH = 0.17753 nm)");
      keys.add("optional","LJEPSILON","TIP3P LJ epsilon parameters (default:   epsilonOO = 0.63681228 kJ/mol\t epsilonHH = 0.19259280 kJ/mol\t epsilonOH = 0.35001648 kJ/mol)");
      keys.add("optional","LJR","Lennard-Jones radii (default: R_min = 0.9*sigmaLJ \t R_max =  0.18 nm");
      keys.add("optional","R_CUT","Verlet lists, default r_cut = R_max");
      keys.add("optional","R_0","For continuous shift to soft-core");
      keys.add("optional","R_SWITCH","Distance at which the long-range switching function starts to act up to R_MAX (default: 0.9*R_MAX)");
      keys.add("optional","MAX_DER","Max gradient value determining soft-core threshold (default: max_der = -400 kJ/mol AA");
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
      
      //CHANGE add angles to each insertion point (ANGLE GRID)


      //log check
      
      log.printf("gridsize:\t%d\n",gsz);


      //TIP3P parameters
      
      //Coulomb parameters
      TIP3Pcharge.resize(2);
      TIP3Pcharge[0] = -0.834;   //Oxygen charge (C)
      TIP3Pcharge[1] =  0.417;   //Hydrogen charge
      paerseVector("CHARGES",TIP3Pcharge);
      C0  = 139.028429916;                  //Electrostatic constant x squared electron charge in kJ/mol*nm --> 332.0637 (kcal/mol*AA) * 4.1868 J->cal * 0.1 AA->nm
      qOO = C0*TIP3Pcharge[0]*TIP3Pcharge[0];  //Product of charges 
      qHH = C0*TIP3Pcharge[1]*TIP3Pcharge[1];  //
      qOH = C0*TIP3Pcharge[0]*TIP3Pcharge[1];  //

      //LJ parameters  
      
      //sigma
      TIP3Psigma.resize(3);
      TIP3Psigma[0] = 0.31507;    //O-O LJ sigma in nm
      TIP3Psigma[1] = 0.04000;    //H-H LJ sigma in nm
      TIP3Psigma[2] = 0.17753;    //O-H LJ sigma in nm
      parseVector("LJSIGMA",TIP3Psigma);
      sigmaOO = TIP3Psigma[0];
      sigmaHH = TIP3Psigma[1];
      sigmaOH = TIP3Psigma[2];

      //epsilon
      TIP3Pepsilon.resize(3);
      TIP3Pepsilon[0] = 0.63681228;  //O-O LJ epsilon in kJ/mol 
      TIP3Pepsilon[1] = 0.19259280;  //H-H LJ epsilon in kJ/mol
      TIP3Pepsilon[2] = 0.35001648;  //O-H LJ epsilon in kJ/mol
      parseVector("LJESPILON",TIP3Pepsilon);
      epsilonOO = TIP3Pepsilon[0];
      epsilonHH = TIP3Pepsilon[1];
      epsilonOH = TIP3Pepsilon[2];

      //LJ AB parameters
      AOO = 4.0*epsilonOO*pow(sigmaOO,12.0);  // A = 4*eps*sig^12
      AHH = 4.0*epsilonHH*pow(sigmaHH,12.0);
      AOH = 4.0*epsilonOH*pow(sigmaOH,12.0);
   
      BOO = 4.0*epsilonOO*pow(sigmaOO,6.0);   // B = 4*eps*sig^6 
      BHH = 4.0*epsilonHH*pow(sigmaHH,6.0);
      BOH = 4.0*epsilonOH*pow(sigmaOH,6.0);

      //LJr.resize(2);
      //LJr[0] = 0.9*sigma_LJ; //nm
      //LJr[1] = 1.25; //
      //parseVector("LJR",LJr);
      R_max = 0.18;   //R_max defined the region for switching to zero the potential tail, default 18 AA

      //Calculate single atom-atom interaction gradient cut-off
      //Can be better done, this is just for testing
      max_der = -400.0;
      parse("MAX_DER",max_der);

      rtest = R_max;
      rstep = 0.01;    //Just for test purpose
      for(;;){
         rtest -= rstep;
         test_dph = dphi_tip3p(rtest,qOO,AOO,BOO);
         test_ph  =  phi_tip3p(rtest,qOO,AOO,BOO);
         if(test_dph > max_der){
           continue;
         }else{
          R_minOO  = rtest;
          V_minOO  = test_ph;
          dV_minOO = test_dph;
          break;
         }
      }

      rtest = R_max;
      rstep = 0.01;    //Just for test purpose
      for(;;){
         rtest -= rstep;
         test_dph = dphi_tip3p(rtest,qHH,AHH,BHH);
         test_ph  =  phi_tip3p(rtest,qHH,AHH,BHH);
         if(test_dph > max_der){
           continue;
         }else{
          R_minHH  = rtest;
          V_minHH  = test_ph;
          dV_minHH = test_dph;
          break;
         }
      }

      rtest = R_max;
      rstep = 0.01;    //Just for test purpose
      for(;;){
         rtest -= rstep;
         test_dph = dphi_tip3p(rtest,qOH,AOH,BOH);
         test_ph  =  phi_tip3p(rtest,qOH,AOH,BOH);
         if(test_dph > max_der){
           continue;
         }else{
          R_minOH  = rtest;
          V_minOH  = test_ph;
          dV_minOH = test_dph;
          break;
         }
      }

      //exponential cutoff (default = 300.0)
      coff=300.0;
      parse("COFF",coff);

      //log check
      log.printf("Cut-off = %.4f kT\n",coff);


      //////////////// CHECK POTENTIALS //////////////
      rtest = 0.1;
      rstep = 0.01;
      ofstream potoo;
      potoo.open("PotOO.dat"); //file 
      ofstream potoo;
      pothh.open("PotHH.dat"); //file 
      ofstream pothh;
      potoh.open("PotOH.dat"); //file 
      for(;;){
         rtest += rstep;
         test_ph     = phi_tip3p(rtest,qOO,AOO,BOO);
         test_ph_qs  = phi_tip3p_qs(rtest,R_minOO,R_max,qOO,AOO,BOO,coff,V_minOO,dV_minOO);
         test_dph    = dphi_tip3p(rtest,qOO,AOO,BOO);
         test_dph_qs = dphi_tip3p_qs(rtest,R_minOO,R_max,qOO,AOO,BOO,coff,V_minOO,dV_minOO);
         potoo << rtest << "\t" << test_ph << "\t" << test_dph << "\t" << test_ph_qs << "\t" << test_dph_qs << endl; //O-O potential

         test_ph     = phi_tip3p(rtest,qHH,AHH,BHH);
         test_ph_qs  = phi_tip3p_qs(rtest,R_minHH,R_max,qHH,AHH,BHH,coff,V_minHH,dV_minHH);
         test_dph    = dphi_tip3p(rtest,qHH,AHH,BHH);
         test_dph_qs = dphi_tip3p_qs(rtest,R_minHH,R_max,qHH,AHH,BHH,coff,V_minHH,dV_minHH);
         pothh << rtest << "\t" << test_ph << "\t" << test_dph << "\t" << test_ph_qs << "\t" << test_dph_qs << endl; //H-H potential

         test_ph     = phi_tip3p(rtest,qOH,AOH,BOH);
         test_ph_qs  = phi_tip3p_qs(rtest,R_minOH,R_max,qOH,AOH,BOH,coff,V_minOH,dV_minOH);
         test_dph    = dphi_tip3p(rtest,qHH,AHH,BHH);
         test_dph_qs = dphi_tip3p_qs(rtest,R_minOH,R_max,qOH,OHH,BOH,coff,V_minOH,dV_minOH);
         potoh << rtest << "\t" << test_ph << "\t" << test_dph << "\t" << test_ph_qs << "\t" << test_dph_qs << endl; //O-H potential

         if(rstep > R_max){
           break;
         }else{
          continue;
         }
      }
      potoo.close();
      pothh.close();
      potoh.close();

      /////////////// END CHECK POTENTIALS ////////////
       
      //log check
      log.printf("TIP3P Coulomb parameters:\tqO = %.4f\tqH = %.4f\n",TIP3Pcharge[0],TIP3Pcharge[1]);
      log.printf("TIP3P LJ sigma parameters:\tsigmaOO = %.4f\tsigmaHH = %.4f\tsigmaOH = %.4f\n",sigmaOO,sigmaHH,sigmaOH);
      log.printf("TIP3P LJ epsilon parameters:\tepsilonOO = %.4f\tepsilonHH = %.4f\tepsilonOH = %.4f\n",epsilonOO,epsilonHH,epsilonOH);
      log.printf("SOFT-CORE radii:\tR_min O-O = %.4f\tR_min H-H = %.4f\tR_min O-H = %.4f\n",R_minOO,R_minHH,R_minOH);
      log.printg("LONG_RANGE switching:\tR_max = %.4f\n",R_max);
      
      //CHANGE add Coulombic potential parameters

	     
      //Potential shift and soft-core keywords

      parseFlag("VSHIFT",vshift);
      
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
      //CHANGE check plumed 2.0 lists and atoms -> molecules 
      r_cut=R_max; 
      parse("R_CUT",r_cut);
      r_skin=1.25*R_max; 
      parse("R_SKIN",r_skin);
      r_switch=0.9*R_max;
      parse("R_SWITCH",r_switch);
      if(r_switch>=R_max){
        log.printf("Error, R_SWITCH>R_CUT, this cannot be!");
      }
      
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
    
    //LJ Potential w/ sigmoidal sof-core
    double muinsw::phis(double rq, double r, double rmin, double rmax, double sigma, double eps, double Vmax, double Vshift, double r0, bool issig, double coff){
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

    //LJ derivative w/ sigmoidal sof-core
    double muinsw::dphis(double rq, double r, double rmin, double rmax, double sigma, double eps, double Vmax, double Vshift, double r0, bool issig, double coff){
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

    //TIP3P Potential w/ quadratic soft-core and sigmoidla long-range switching
    double muinsw::phi_tip3p_qs(double r, double rmin, double rmax, double qq, double Apar, double Bpar, double coff, double Vmin, double dmin){
      double phi, ph, r6, r12;
      if(r <= rmin){
        phi=0.5*dmin*(r*r/rmin-rmin)+Vmin; 
      }else if(r >= rmax){
        r6  = pow(r,6.0);
        r12 = pow(r,12.0);
        ph  = qq/r + Apar/r12 - Bpar/r6;
        double Z0=0.5*(r_switch+rmax);
        double w=0.05*(rmax-r_switch);
        double z=(r-Z0)/w;
        double soff=swfoff(z,coff);
        phi=ph*soff;
      }else{
        r6  = pow(r,6.0);
        r12 = pow(r,12.0);
        phi = qq/r + Apar/r12 - Bpar/r6;
      }
      return(phi);
    }

    //TIP3P Potential derivative w/ quadratic soft-core and sigmoidla long-range switching
    double muinsw::dphi_tip3p_qs(double r, double rmin, double rmax, double qq, double Apar, double Bpar, double coff, double Vmax, double dmin){
      double dphi, dph, ph, r6, r7, r12, r13;
      if(r <= rmin){
        dphi=dmin*r/rmin;
      }else if(r >= rmax){
        r2  = pow(r,2.0);
        r6  = pow(r,6.0);
        r7  = pow(r,7.0);
        r12 = pow(r,12.0);
        r13 = pow(r,13.0);
        ph  = qq/r + Apar/r12 - Bpar/r6;
        dph = -qq/r2 - 12.0*Apar/r13 + 6.0*Bpar/r7
        double Z0=0.5*(r_switch+rmax);
        double w=0.05*(rmax-r_switch);
        double z=(r-Z0)/w;
        double soff=swfoff(z,coff);
        double ds=dswf(z,coff)*w;
        dphi=ph*ds+dph*soff;
      }else{
        r2  = pow(r,2.0);
        r7  = pow(r,7.0);
        r13 = pow(r,13.0);
        dphi = -qq/r2 - 12.0*Apar/r13 + 6.0*Bpar/r7
      }
      return(dphi);
    }


    //TIP3P test for cut-off calculation
    double muinsw::phi_tip3p(double r, double qq, double Apar, double Bpar){
      double phi, r6, r12;
      r6  = pow(r,6.0);
      r12 = pow(r,12.0);
      phi = qq/r + Apar/r12 - Bpar/r6;
      return(phi);
    }
                                      

    //TIP3P derivative test for cut-off calculation
    double muinsw::dphi_tip3p(double r, double qq, double Apar, double Bpar){
      double dphi, r2, r7, r13;
      r2  = pow(r,2.0);
      r7  = pow(r,7.0);
      r13 = pow(r,13.0);
      dphi = -qq/r2 - 12.0*Apar/r13 + 6.0*Bpar/r7
      return(dphi);
    }
    
    

    // verlet list formation
    //CHANGE cycle on molecule number (???)    

    void muinsw::du_newlist(vector<AtomNumber> &list, vector<Vector> &grid, dulist_s &dulist)
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
	  mod_rij=rij.modulo();
	  if (mod_rij < dulist.rskin){ //if distance < rskin
	    dulist.ni[j][dulist.nn[j]]=i; //index in neighlist
	    dulist.nn[j]++; //increment nn 
	  }	  
	}
      }
    }
    
    //CHANGE cycle on molecule number (???)
    void muinsw::du_checklist(vector<AtomNumber> &list, vector<Vector> &grid, dulist_s &dulist)
    {
      unsigned int j; 
      double dr=(dulist.rskin-dulist.rcut)*0.5;
      Vector rij;
        
      for (j=0; j<list.size(); ++j) { 
	rij = pbcDistance(getPosition(j),dulist.posit[j]); //check position variations
       	if( fabs(rij[0])>dr || fabs(rij[1])>dr || fabs(rij[2])>dr ) { 
	  id++; //rebuild list counter
	  du_newlist(list,grid,dulist); 
	  break; 
	}
      }
    } 
    
    // calculator
    void muinsw::calculate()    
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
      double dVi=1./gsz;  //volume element, in [Vbox] units
      //Space grid
      
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
      
      
      for(i=0; i<griddim[0]; ++i){
	  for(j=0; j<griddim[1]; ++j){
	      for(k=0; k<griddim[2]; ++k){
		index=griddim[1]*griddim[2]*i+griddim[2]*j+k; //index
		gs=Vector(di[0]*i,di[1]*j,di[2]*k);
		xgrid[index]=matmul(transpose(getBox()),gs);
	      }
	    }
	  }
     
      //CHANGE define possible insertion configurations (per grid point) 


//Vshift is not used (fermi-swtiching funct instead, Vmax defined at the beginning
      
//      //initialize potential 
//      double Vshift;
//      if(vshift){
//       Vshift=phi(R_max*R_max, 0.0, R_max+1., sigma_LJ, eps_LJ, 0.0, 0.0);
//      }else{
//	Vshift=0;
//      }
//      
//      double Vmax;
//      if(r0>0.0){
//        Vmax=phi(r0*r0, 0.0, R_max, sigma_LJ, eps_LJ, 0.0, Vshift);
//      }else{
//        Vmax=phi(R_min*R_min, 0.0, R_max, sigma_LJ, eps_LJ, 0.0, Vshift);
//      }

      double dmin;
      if(quad) dmin=dphi(R_min*R_min, R_min, 0.0, R_max, sigma_LJ, eps_LJ);

      //CHANGE Coulomb initialization
      
      //create neigh lists  
      
      if(dulist.step==0) du_newlist(at_list,xgrid,dulist);
      ++dulist.step; 
      du_checklist(at_list,xgrid,dulist);
      
      double Su=0;
      double Du,modrij,modq,zu,fDu,dfDu, ph;
      Vector rij;
      int atomid;
      double beta=1./kT;
      
      //cycle on space grid

      for(i=rank; i<gsz; i+=stride){
      //CHANGE loop on configurations      
         
	//CHANGE loop on test palticle atoms 
	
	//auxiliary array of vectors for derivatives
	//CHANGE if neighbor list contains molecules this has to be changed accordingly
	vector<Vector> phiprime(dulist.nn[i]);
	//auxiliary array of tensors for derivatives
	vector<Tensor> phit(dulist.nn[i]);
	
	//auxiliary 
	Du=0; //zero if no wall is present
	

	//cycle over neighbor particles
	//CHANGE (molecules or atoms)
	
	for(j=0; j<dulist.nn[i]; ++j){
	  if(beta*Du>coff)   break;

	  //CHANGE watch out if the loop is over atoms or molecules
	  atomid = dulist.ni[i][j]; //atomindex
	  //Delta U
	  rij = pbcDistance(xgrid[i],getPosition(atomid));
	  modrij=rij.modulo();
	  modq=modrij*modrij;
	  
	  //CHANGE once i and j are associated to atoms, set the sigma,epsilon and Coulombic parameters
	  //CHANGE ph and dph contain all the contributions
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

	  //CHANGE j is now referred to molecules, change phiprime and phit to the proper atom index
	  phiprime[ ] = rij*dph/modrij; //potential derivative
	  phit[ ] = Tensor(phiprime[ ],rij); //Virial component
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
	double dSoV=dfDu*dVi*grid_w[i];
	
	//derivatives -- cycle again over neighbors
	//CHANGE cycle on all atoms of the neighborhood list
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
