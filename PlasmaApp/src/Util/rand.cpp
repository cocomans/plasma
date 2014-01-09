#include <cmath>
extern "C"
{
#define MODULUS    2147483647 /* DON'T CHANGE THIS VALUE                  */
#define MULTIPLIER 48271      /* DON'T CHANGE THIS VALUE                  */
#define CHECK      399268537  /* DON'T CHANGE THIS VALUE                  */
#define STREAMS    256        /* # of streams, DON'T CHANGE THIS VALUE    */
#define A256       22925      /* jump multiplier, DON'T CHANGE THIS VALUE */
#define DEFAULT    123456789  /* initial seed, use 0 < DEFAULT < MODULUS  */

int seedt = 1;

void   Randoms(double &rx, int &seed)
/* ----------------------------------------------------------------
 * Random returns a pseudo-random real number uniformly distributed 
 * between 0.0 and 1.0. 
 * ----------------------------------------------------------------
 */
{
  const long Q = MODULUS / MULTIPLIER;
  const long R = MODULUS % MULTIPLIER;
        long t;

  t = MULTIPLIER * (seed % Q) - R * (seed / Q);
  if (t > 0) 
    seed = t;
  else 
    seed = t + MODULUS;
  rx = ((double) seed / MODULUS);
}

  //Lehmer random number generator
  //x(t+1) = (a*x(t)) mod m
  void lcg_rand(double &rx, int &seed)
  {
    const double a =  279470273.0;
    const double m = 4294967291.0;
    double temp = a * seed;
    seed = (int)(fmod(temp,m));
    rx = double(seed)/m;
  }


  //Carbon linear congruential generator
  //X[n+1] = (a*X[n] + c) mod m 
  //with a = 16807, c = 0, and m = 2^31 - 1
  //  double LCG_rand( int &seed ) {
  void LCG_rand(double &rx, int &seed ) {
    // Input
    //   seed    Integer seed (DO NOT USE A SEED OF ZERO)
    // Output
    //   rand    Random number uniformly distributed in [0,1)

    const double a = 16807.0;
    const double m = 2147483647.0;
    double temp = a * seed;
    seed = (int)(fmod(temp,m));
    rx = double(seed)/m;
    //return(seed/m);
  }
  /*--- end LCG_rand() ---*/  

  //Random number generator; Normal (Gaussian) dist.
  //Box-Muller transformation
  //  void randn(double &rx) {
  void randn(double &rx) {
    // Input
    //   seed    Integer seed  (DO NOT USE A SEED OF ZERO)
    // Output
    //	 randn   Random number, Gaussian distributed
    //static int seedn = 1;
    static int use_last = 0;
    static double rn2;
    double x1,x2,w;
    double lrx;
    // if(use_last){
    //   use_last = 0;
    //   rx = rn2;
    // }else{
      do{
	LCG_rand(lrx,seedt);
	x1 = 2.0*lrx - 1.0;
	LCG_rand(lrx,seedt);
	x2 = 2.0*lrx - 1.0;
	w = x1*x1 + x2*x2;
      } while ( w>=1.0 );
    
      w = sqrt( (-2.0*log(w) ) / w );
      // rn2 = x2*w;
      // use_last=1;
      rx = x2*w;
    // }
  }
  /*--- end randn() ---*/  

}
