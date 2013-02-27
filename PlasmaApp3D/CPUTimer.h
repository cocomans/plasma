/*-------------------------------------------------------------------------*/
/**
	@file		CPUTimer.h
	@author	J. Payne
	@date		12/21/2012
	@brief	Nanosecond timer for the CPU based on clock_gettime()

*/
/*--------------------------------------------------------------------------*/

#ifndef CPU_TIMER_H
#define CPU_TIMER_H
#include <iostream>
#include <sys/time.h>
#include <math.h>
#ifndef __MACH__
#include <sys/times.h>
#include <sys/vtimes.h>
#endif
#include <stdint.h>

/*-------------------------------------------------------------------------*/
/**
	@brief Wrapper for the RDTSC instruction
*/
/*--------------------------------------------------------------------------*/
static inline uint64_t RDTSC(void) {

	  unsigned int hi, lo;

	  __asm__ volatile("rdtsc" : "=a" (lo), "=d" (hi));

	  return ((uint64_t)hi << 32) | lo;

}

#ifdef __MACH__
#include <mach/mach_time.h>
#define ORWL_NANO (+1.0E-9)
#define ORWL_GIGA UINT64_C(1000000000)

static double orwl_timebase = 0.0;
static uint64_t orwl_timestart = 0;

struct timespec orwl_gettime(void) {
  // be more careful in a multithreaded environement
  if (!orwl_timestart) {
    mach_timebase_info_data_t tb = { 0 };
    mach_timebase_info(&tb);
    orwl_timebase = tb.numer;
    orwl_timebase /= tb.denom;
    orwl_timestart = mach_absolute_time();
  }
  struct timespec t;
  double diff = (mach_absolute_time() - orwl_timestart) * orwl_timebase;
  t.tv_sec = diff * ORWL_NANO;
  t.tv_nsec = diff - (t.tv_sec * ORWL_GIGA);
  return t;
}
#endif

const int NANO_SECONDS_IN_SEC = 1000000000;

/*-------------------------------------------------------------------------*/
/**
	@class CPUTimer
	@author	J. Payne
	@date		12/21/2012
	@brief	High resolution timer class

	This class is used for nanosecond timing and determining cpu frequency.
	When constructed the timer calibrates the number of cpu clock cycles per
	 NanoSecond. This value is then used as a calibration for the frequency calculations.

*/
/*--------------------------------------------------------------------------*/
class CPUTimer
{
public:
	/// Start stime
	timespec begin;
	/// End time
	timespec end;

	/// @brief Cummulative time since first call of start()
	timespec cummulative;

	/// @name counters for measuring clock frequency
	//@{
	uint64_t ts0, ts1, ts_cumm;
	//@}

	/// Calibration for measuring clock frequency
	double g_TicksPerNanoSec;

	double bogus;

	/*-------------------------------------------------------------------------*/
	/**
		@brief Constructs a CPUTimer object and calibrates it for measuring clock
		frequency.

	*/
	/*--------------------------------------------------------------------------*/
	CPUTimer(){cummulative.tv_nsec =0;cummulative.tv_sec=0;ts_cumm = 0;
	CalibrateTicks();}

	/*-------------------------------------------------------------------------*/
	/**
		@brief Calibrates the clock frequency measurement utility of the timer
	*/
	/*--------------------------------------------------------------------------*/
	void CalibrateTicks()
	{

	  struct timespec begints, endts;

	  uint64_t begint = 0, endt = 0;
#ifdef __MACH__
	  begints = orwl_gettime();
#else
	  clock_gettime(CLOCK_REALTIME, &begints);
#endif
	  begint = RDTSC();

	  uint64_t i;
	  double temp = 0;
	  for (i = 0; i < NANO_SECONDS_IN_SEC/50; i++) temp+=sqrt((double)i); /* must be CPU intensive */

	  bogus = temp;

	  endt = RDTSC();
#ifdef __MACH__
	  endts = orwl_gettime();
#else
	  clock_gettime(CLOCK_REALTIME, &endts);
#endif
	  struct timespec tmpts = TimeSpecDiff(&endts, &begints);

	  uint64_t nsecElapsed = tmpts.tv_sec * 1000000000LL + tmpts.tv_nsec;

	  g_TicksPerNanoSec = (double)(endt - begint)/(double)nsecElapsed;

	}

	/*-------------------------------------------------------------------------*/
	/**
		@brief Internal function used to measure the difference in time between
		two timespecs, ts = ts1 - ts2
		\param[in] ts1 First timespec
		\param[in] ts2 Second timespec
		\return ts \sa ts1 - ts2
	*/
	/*--------------------------------------------------------------------------*/
	struct timespec TimeSpecDiff(struct timespec *ts1, struct timespec *ts2)
	{

	timespec ts;

	ts.tv_sec = ts1->tv_sec - ts2->tv_sec;

	ts.tv_nsec = ts1->tv_nsec - ts2->tv_nsec;

	if (ts.tv_nsec < 0) {

	ts.tv_sec--;

	ts.tv_nsec += NANO_SECONDS_IN_SEC;

	}

	return ts;

	}


	/*-------------------------------------------------------------------------*/
	/**
		@brief Function used to start the timer.

		Calls clock_gettime() to get the current cpu time for the current process ID.
	*/
	/*--------------------------------------------------------------------------*/
	void start()
	{
//		clock_gettime(CLOCK_REALTIME, &begin);
		ts0 = RDTSC();
	}

	/*-------------------------------------------------------------------------*/
	/**
		@brief Function used to stop the timer

		Calls clock_gettime() to get the end time, adds the elapsed time between
		 begin and end to the cumulative time.
	*/
	/*--------------------------------------------------------------------------*/
	void stop()
	{
		ts1 = RDTSC();
		ts_cumm += ts1-ts0;

//		clock_gettime(CLOCK_REALTIME, &end);
//
//		timespec tempt = tdiff();
//		double temp_sec;
//
//		cummulative.tv_sec += tempt.tv_sec;
//		cummulative.tv_nsec += tempt.tv_nsec;
//
////		ts_cumm += ts1-ts0;
//
//		temp_sec = floor(cummulative.tv_nsec/1.0e9);
//
//		cummulative.tv_nsec -= temp_sec*1.0e9;
//		cummulative.tv_sec += temp_sec;

	}

	/*-------------------------------------------------------------------------*/
	/**
		@brief Returns the cumulative time in milliseconds.
	*/
	/*--------------------------------------------------------------------------*/
	double get_cummulative()
	{
//		double temp1,temp2;
//		temp1 = cummulative.tv_nsec*1.0e-6 + cummulative.tv_sec*1.0e3;
//		temp2 = 1.0e-6*((double)ts_cumm)/g_TicksPerNanoSec;
//
//		printf("*************Timer difference is %e\n",fabs(temp1-temp2)/temp1);

//		return cummulative.tv_nsec*1.0e-6 + cummulative.tv_sec*1.0e3;
		return 1.0e-6*((double)ts_cumm)/g_TicksPerNanoSec;
	}

	/*-------------------------------------------------------------------------*/
	/**
		@brief Returns the difference between begin and end in microseconds.
	*/
	/*--------------------------------------------------------------------------*/
	double diff()
	{
		// Return the time difference in milliseconds
		timespec temp;
		if ((end.tv_nsec-begin.tv_nsec)<0) {
			temp.tv_sec = end.tv_sec-begin.tv_sec-1;
			temp.tv_nsec = NANO_SECONDS_IN_SEC+end.tv_nsec-begin.tv_nsec;
		} else {
			temp.tv_sec = end.tv_sec-begin.tv_sec;
			temp.tv_nsec = end.tv_nsec-begin.tv_nsec;
		}
		return temp.tv_nsec*1.0e-6 + temp.tv_sec*1.0e3;
	}

	/*-------------------------------------------------------------------------*/
	/**
		@brief Returns the difference between begin and end as a timespec
	*/
	/*--------------------------------------------------------------------------*/
	timespec tdiff()
	{
		// Return the time difference in nanosecnods
		timespec temp;
		if ((end.tv_nsec-begin.tv_nsec)<0) {
			temp.tv_sec = end.tv_sec-begin.tv_sec-1;
			temp.tv_nsec = NANO_SECONDS_IN_SEC+end.tv_nsec-begin.tv_nsec;
		} else {
			temp.tv_sec = end.tv_sec-begin.tv_sec;
			temp.tv_nsec = end.tv_nsec-begin.tv_nsec;
		}
		return temp;
	}

	/*-------------------------------------------------------------------------*/
	/**
		@brief Returns the cpu clock frequency averaged over all timed periods
	*/
	/*--------------------------------------------------------------------------*/
	double getFrequency()
	{
		return ((double)ts_cumm)/(1.0e-3*get_cummulative()) * 1.0e-9;
	}



};

struct CpuInfo {
	char vendor_id[50];
	int family;
	char model[50];
	float freq;
	char cache[20];
};



class CPUSpeed
{
public:



};









#endif /* TIMER_H */
