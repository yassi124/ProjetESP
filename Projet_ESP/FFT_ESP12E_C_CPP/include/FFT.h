#ifndef FFT_h /* Prevent loading library twice */
#define FFT_h
#ifdef ARDUINO
	#if ARDUINO >= 100
		#include "Arduino.h"
	#else
		#include "WProgram.h" /* Standard arduino  */
	#endif
#else
	#include <stdlib.h>
	#include <stdio.h>
	#ifdef __AVR__
		#include <avr/io.h>/* For avr microControllers  */
	#endif
	#include <math.h>
#endif


/* Custom constants */
#define FFT_FORWARD 0x01
#define FFT_REVERSE 0x00

/* Windowing type */
#define FFT_WIN_TYP_RECTANGLE 0x00 /* rectangle (Box car) */
#define FFT_WIN_TYP_HAMMING 0x01 /* hamming */
#define FFT_WIN_TYP_HANN 0x02 /* hann */
#define FFT_WIN_TYP_TRIANGLE 0x03 /* triangle (Bartlett) */
#define FFT_WIN_TYP_BLACKMAN 0x04 /* blackmann */
#define FFT_WIN_TYP_FLT_TOP 0x05 /* flat top */
#define FFT_WIN_TYP_WELCH 0x06 /* welch */
/*Mathematial constants*/
#define twoPi 2*PI
#define fourPi 4*PI

class FFT {
public:
	/* Constructor */
	FFT(void);
	/* Destructor */
	~FFT(void);
	/* Functions */

	void FFT_PACK(double* rData,uint16_t num_samples,uint8_t window,double* fPeak,double* APeak);


private:
	/* Variables */
	//uint16_t _samples;
	//double _samplingFrequency;
	//double *_vReal;
	//double *_vImag;
	//uint8_t _power;
	/* Functions */
	void Swap(double *x, double *y);
	void Compute(double *vReal, double *vImag, uint16_t samples, uint8_t power, uint8_t dir);
	void Compute(double *vReal, double *vImag, uint16_t samples, uint8_t dir);
	uint8_t Exponent(uint16_t value);
	void ComplexToMagnitude(double *vReal, double *vImag, uint16_t samples);
	void MajorPeak(double *vD, uint16_t samples, double samplingFrequency,double* fqPeak, double* aPeak);
	void Windowing(double *vData, uint16_t samples, uint8_t windowType, uint8_t dir);
	void ComplexToMagnitude();
};

#endif
