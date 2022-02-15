/*--------------------------------------------------------------------*/
// for output while running:
#include <signal.h>
#include <unistd.h>
#include <iomanip>
extern int elapsed;
extern float percentage;

static void sigalrm_handler( int sig ) {
  elapsed+=1; // monitor progress every 10 sec
  cout.flush();
  cout << "  "<< setw(2) << setfill('0') << elapsed/60 ;
  cout << ":" ;
  cout << setw(2) << setfill('0') << elapsed%60 <<" "; //<< right << 
  //cout << setw(10) << "K = ("<< setprecision(2) << k0 << "," << k << ")" << "  ";
  cout << '[';
  for (int i=0;i<50;i++) {
    if ((double)i<percentage*50.) {cout << '#';}
    else {cout << '-';};
  }
  cout << "] " << setw(2) << setfill(' ') << (int)(percentage*100.);
  cout << "%" << '\r';
  alarm(1);
}

/*--------------------------------------------------------------------*/
