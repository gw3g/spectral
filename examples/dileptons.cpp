#include "core.hh"
#include "quad.hh"
#include "map.hh"
#include <fstream>
#include <string>

double k0,k;

using namespace std;

ofstream fout;
int Print_D(double); int elapsed; float percentage;

int main(int argc, char *argv[]) {
  /* for lattice comparisons: */
  // nf=0
  //if (argc>1) { k=atof(argv[1]); Print_D(k); }
  //Print_D(3.*2.*M_PI/3.); // T=1.1Tc
  //Print_D(3.*7.*M_PI/12.); // T=1.3Tc
  // nf=2
  //Print_D(sqrt(1.)*M_PI/2.); // T=1.2Tc
  //Print_D(1.5*M_PI); // T=1.2Tc
  Print_D(.5*M_PI*sqrt(14)); // T=1.2Tc

  //Print_D(.5);
  //Print_D(1.);
  //Print_D(1.5);

  //Print_D(.3);
  //Print_D(1.5);
  //Print_D(3.);
  //Print_D(6.);
  //Print_D(9.);
}

/*--------------------------------------------------------------------*/

struct Rho_V
{
  double lo, nlo;
  int S[3];
  double Nc=3., cF=(Nc*Nc-1.)/(2.*Nc); // group factors

  // notation à la 1310.0164
  Master
    *rho_b, *rho_bb, *rho_d, *rho_db, *rho_g, *rho_h, *rho_hp, *rho_j;
  double
    _b, _bb, _d, _db, _g, _h, _hp, _j;

  Rho_V() {

    S[0] = +1; S[1] = - 1; S[2] = -1; // statistics

    rho_b =  _11020(0,0,S); // 'vector' channel
    rho_bb=  _10120(0,0,S);
    rho_d =  _11010(0,0,S);
    rho_db=  _10110(0,0,S);
    rho_g =  _11011(0,0,S);
    rho_h =  _11110(0,0,S);
    rho_hp=  _Star( 0,0,S);
    rho_j =  _11111(0,0,S);

  };
  void operator ()() {

    lo = +Nc*K2*psi0(-1,-1)*OOFP; // note sign convention
    //lo = Nc*K2*OOFP; // large-K2
    nlo= lo*3.*cF*SQR(OOFP);

    double res, err;
    ////////////////////////////////////
  double epsabs = 1e-3, epsrel = 0;
  double rr = 1.;

  gsl_set_error_handler_off();
  size_t limit = 5e2;

  quad wsp1(limit);
  quad wsp2(limit);

  cout << 
          2.*(rho_h->integrand)(.1,.3)+
          2.*(rho_hp->integrand)(.1,.3)*(kp/km)+
          -(rho_j->integrand)(.1,.3)*(kp/km)
          << endl;


  auto outer = make_gsl_function( [&](double x) 
  {
    double inner_result, inner_abserr;
    auto inner = make_gsl_function( [&](double y) {
          return
          2.*(rho_h->integrand)(x,y)+
          2.*(rho_hp->integrand)(x,y)*(kp/km)+
          -(rho_j->integrand)(x,y)*(kp/km)
          ;
        } );
    gsl_integration_qag( inner, .0+1e-10,1., epsabs, epsrel,
                         limit, 6, wsp1, &inner_result, &inner_abserr );
    return inner_result;
  } );
  gsl_integration_qag( outer, .0+1e-10,1., epsabs, epsrel*2,
                       limit, 6, wsp2, &res, &err  );//*/
    ////////////////////////////////////

    _b = (*rho_b )(k0,k)*K2;
    _bb= (*rho_bb)(k0,k)*K2;
    _d = (*rho_d )(k0,k);
    _db= (*rho_db)(k0,k);
    _g = (*rho_g )(k0,k)*K2;
    //_h = (*rho_h )(k0,k)*K2;
    //_hp= (*rho_hp)(k0,k);
    //_j = (*rho_j )(k0,k)*SQR(K2);

    nlo -=
    //8.*Nc*cF*( 2.*(_b-_bb+_d-_db) - 1.5*_g + 2.*(_h+_hp) - _j );
    8.*Nc*cF*( 2.*(_b-_bb+_d-_db) - 1.5*_g + res*CUBE(OOFP) );

  };
};

struct Rho_00
{
  double lo, nlo;
  int S[3];
  double Nc=3., cF=(Nc*Nc-1.)/(2.*Nc);

  Master
    *rho_b_0, *rho_bb_0,
    *rho_b_1, *rho_bb_1,
    *rho_b_2, *rho_bb_2,
    *rho_g,
    *rho_h_0, *rho_h_1,
    *rho_hp,
    *rho_j_0, *rho_j_2;
  double
    _b_0, _bb_0,
    _b_1, _bb_1,
    _b_2, _bb_2,
    _g,
    _h_0, _h_1,
    _hp,
    _j_0, _j_2;

  Rho_00() {

    S[0] = +1; S[1] = - 1; S[2] = -1; // statistics

    //rho_b_0 =  _11020(0,0,S); // notation: rho_<tag>_<power of ..>
    //rho_bb_0=  _10120(0,0,S);
    //rho_b_1 =  _11020(1,0,S);
    //rho_bb_1=  _10120(1,0,S);
    //rho_b_2 =  _11020(2,0,S);
    //rho_bb_2=  _10120(2,0,S);
    rho_g   =  _11011(0,0,S);
    rho_h_0 =  _11110(0,0,S);
    rho_h_1 =  _11110(0,1,S);
    rho_hp  =  _Star( 0,0,S);
    rho_j_0 =  _11111(0,0,S);
    rho_j_2 =  _11111(2,0,S);

  };
  void operator ()() {

    lo = -2.*Nc*( k0*k0*(psi1(-1,-1)-psi2(-1,-1))-.25*K2*psi0(-1,-1) )*OOFP;
    nlo = -cF*Nc*( k*k*psi0(-1,-1) )*CUBE(OOFP);
    //nlo = -cF*Nc*( k*k )*CUBE(OOFP); // large-K2
    double res, err;
    ////////////////////////////////////
  double epsabs = 1e-3, epsrel = 0;
  double rr = 1.;

  gsl_set_error_handler_off();
  size_t limit = 5e2;

  quad wsp1(limit);
  quad wsp2(limit);

  cout <<
          2*(rho_h_0->integrand)(.1,.3)+
          2*(rho_hp->integrand)(.1,.3)*(kp/km)+
          -8*(rho_h_1->integrand)(.1,.3)*k0*k0/K2
          +  (rho_j_0->integrand)(.1,.3)*(k0*k0+k*k)*(kp/km)/K2
          -4*(rho_j_2->integrand)(.1,.3)*(kp/km) << endl;

  auto outer = make_gsl_function( [&](double x) 
  {
    double inner_result, inner_abserr;
    auto inner = make_gsl_function( [&](double y) {
          return
          2*(rho_h_0->integrand)(x,y)+
          2*(rho_hp->integrand)(x,y)*(kp/km)+
          -8*(rho_h_1->integrand)(x,y)*k0*k0/K2
          +  (rho_j_0->integrand)(x,y)*(k0*k0+k*k)*(kp/km)/K2
          -4*(rho_j_2->integrand)(x,y)*(kp/km)
          ;
        } );
    gsl_integration_qag( inner, .0+1e-10,1., epsabs, epsrel,
                         limit, 6, wsp1, &inner_result, &inner_abserr );
    return inner_result;
  } );
  gsl_integration_qag( outer, .0+1e-10,1., epsabs, epsrel*2,
                       limit, 6, wsp2, &res, &err  );//*/
    ////////////////////////////////////

    //_b_0 = (*rho_b_0 )(k0,k)*K2;
    //_bb_0= (*rho_bb_0)(k0,k)*K2;
    //_b_1 = (*rho_b_1 )(k0,k)*k0;
    //_bb_1= (*rho_bb_1)(k0,k)*k0;
    //_b_2 = (*rho_b_2 )(k0,k);
    //_bb_2= (*rho_bb_2)(k0,k);
    _g   = (*rho_g   )(k0,k)*k*k;
    //_h_0 = (*rho_h_0 )(k0,k)*K2;
    //_h_1 = (*rho_h_1 )(k0,k)*k0;
    //_hp  = (*rho_hp  )(k0,k);
    //_j_0 = (*rho_j_0 )(k0,k)*(k0*k0+k*k)*K2;
    //_j_2 = (*rho_j_2 )(k0,k)*K2;

    nlo -=
    4.*Nc*cF*( //2.*(_b_0-_bb_0-4.*(_b_1-_bb_1)+4.*(_b_2-_bb_2)) // =0 
             //+ _g + 2.*(_h_0+_hp) - 8.*_h_1 + _j_0 - 4.*_j_2 );
             + _g + res*CUBE(OOFP) );

  };
};


int Print_D(double k_curr) {
  int N_k0;
  double res, s, k0_min, k0_max;
  k=k_curr;

  // filename
  char k_name[20];
  sprintf(k_name,"{k=%.2f}",k);
  string fname = "NLO_rho_"
               + string(k_name)
               + ".dat";

  cout << "\n:: Creating table for k = " << k <<  " ..." << endl << endl;
  Rho_V rV;
  Rho_00 r00;
  fout.open(fname);
  fout << 
  "# Columns: k0/T, rhoV_LO/T2, rho00_LO/T2, rhoV_NLO/(g2*T2), rho00_NLO/(g2*T2)" 
       << endl
       << "# ( k=" << k << " )" << endl;

  signal( SIGALRM, sigalrm_handler );
  elapsed=0; alarm(1);

  // Here are some parameters that can be changed:
  N_k0=10; 

  k0_min=15.+1e-1;
  k0_max=30.;
  //k0_min=40.;
  //k0_max=50.;
  // don't change anything after that.

  //s=pow(k0_max/k0_min,1./(N_k0-1));
  //s=(k0_max-k0_min)/((double)N_k0-1.);
  s = 1e-1;
  k0=k0_min;

  int i=0;
  //for (int i=0;i<N_k0;i++) { 
  while (k0<k0_max) {
    //percentage=(float)i/((float)N_k0);
    percentage = k0/k0_max;

    rV();
    r00();

    fout << scientific << k0        // k0/T
         <<     "    " << rV.lo     // leading-order: rho_V ,
         <<     "    " << r00.lo    //                rho_00
         <<     "    " << rV.nlo    // next-to-LO   : rho_V ,
         <<     "    " << r00.nlo   //                rho_00
         << endl;

    k0+=s; 
  }
  cout << endl << ":: Saved to file [" << fname << "]" << endl;
  fout.close();

  return 0;
}

