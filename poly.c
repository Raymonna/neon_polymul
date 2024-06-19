#include <stdint.h>
#include <stdio.h>
#include "params.h"

#include "poly.h"

int32_t mu = (1UL << 16) / NTRUP_Q;

static int16_t cmod(int32_t a, int16_t mod){
    int16_t t;
    t = a % mod;
    if(t > (mod >> 1)){
        t -= mod;
    }
    if(t < -(mod >> 1)){
        t += mod;
    }
    return t;
}
int32_t m = (int32_t)NTRUP_Q;
static int16_t bred(int32_t a) {
    int32_t q = ((int32_t)a * mu) >> 16;
    int32_t t = a - q * m;
    if (t > (m >> 1)) {
        t -= m;
    }
    if (t < -(m >> 1)) {
        t += m;
    }
    return (int16_t)t;
}


void poly_Rq_mul_small(int16_t *h, const int16_t *f,const int8_t *g)
{
    int16_t fg[NTRUP_P + NTRUP_P - 1];
    int32_t result, r1, r2, r3, r4, r5, r6, r7, r8;
	int i,j;
	//h[0]
	fg[0] = f[0]*g[0];
	fg[1] = f[1]*(int32_t)g[0] + f[0]*(int32_t)g[1];
	fg[2] = bred(f[2]*(int32_t)g[0] + f[1]*(int32_t)g[1] + f[0]*(int32_t)g[2]);
	fg[3] = bred(f[3]*(int32_t)g[0] + f[2]*(int32_t)g[1] + f[1]*(int32_t)g[2] + f[0]*(int32_t)g[3]);
	fg[4] = bred(f[4]*(int32_t)g[0] + f[3]*(int32_t)g[1] + f[2]*(int32_t)g[2] 
				+ f[1]*(int32_t)g[3] + f[0]*(int32_t)g[4]);
	fg[5] = bred(f[5] * (int32_t)g[0] + f[4] * (int32_t)g[1] + f[3] * (int32_t)g[2] + 
				f[2] * (int32_t)g[3] + f[1] * (int32_t)g[4] + f[0] * (int32_t)g[5]);
	fg[6] = bred(f[6] * (int32_t)g[0] + f[5] * (int32_t)g[1] + f[4] * (int32_t)g[2] + 
				f[3] * (int32_t)g[3] + f[2] * (int32_t)g[4] + f[1] * (int32_t)g[5] + 
				f[0] * (int32_t)g[6]);
	int16_t rem = 0;
	for (i = 7; i < NTRUP_P; i++) {
      result = 0;
      for (j = 0;j <= i;j+=7){
	  	  rem = i - j;
		  if(rem < 6){
			while(rem >= 0){	
				result = result + (int32_t)(f[i-rem]*(int32_t)g[rem]);
				rem--;
			}	
			result = bred(result);
			break;
		  }
	  	  r1 = (int32_t)f[j]*(int32_t)g[i-j];
	  	  r2 = (int32_t)f[j+1]*(int32_t)g[i-j-1];
	  	  r3 = (int32_t)f[j+2]*(int32_t)g[i-j-2];
	  	  r4 = (int32_t)f[j+3]*(int32_t)g[i-j-3];

	  	  r5 = (int32_t)f[j+4]*(int32_t)g[i-j-4];
	  	  r6 = (int32_t)f[j+5]*(int32_t)g[i-j-5];
	  	  r7 = (int32_t)f[j+6]*(int32_t)g[i-j-6];
	  	  //r8 = (int32_t)(f[j+7]*(int32_t)g[i-j-7]);
		  		   
          result = bred(result+r1+r2+r3+r4+r5+r6+r7);
	  }
      fg[i] = result;

    }
	
	/*
	for (i = 0; i < NTRUP_P; i++) {
      result = 0;
      for (j = 0;j <= i;++j) 
          result = bred(result+f[j]*(int32_t)g[i-j]);
      fg[i] = result;
    }
	*/
	/*
	for (i = NTRUP_P;i < 2 * NTRUP_P - 1;++i) {
      result = 0;
      for (j = i - NTRUP_P + 1; j < NTRUP_P; ++j)
          result = bred(result+f[j]*(int32_t)g[i-j]);
      fg[i] = result;
    }
    */
	for (i = NTRUP_P;i < 2 * NTRUP_P - 1;++i) {
      result = 0;
      for (j = i - NTRUP_P + 1; j < NTRUP_P; j+=7) { 
		  if(j > NTRUP_P-7){
			while(j < NTRUP_P){	
				result = result + (int32_t)(f[j]*(int32_t)g[i - j]);
				j++;
			}	
			result = bred(result);
			break;
		  }
	  	  r1 = (int32_t)f[j]*(int32_t)g[i-j];
	  	  r2 = (int32_t)f[j+1]*(int32_t)g[i-j-1];
	  	  r3 = (int32_t)f[j+2]*(int32_t)g[i-j-2];
	  	  r4 = (int32_t)f[j+3]*(int32_t)g[i-j-3];

	  	  r5 = (int32_t)f[j+4]*(int32_t)g[i-j-4];
	  	  r6 = (int32_t)f[j+5]*(int32_t)g[i-j-5];
	  	  r7 = (int32_t)f[j+6]*(int32_t)g[i-j-6];
	  	  //r8 = (int32_t)(f[j+7]*(int32_t)g[i-j-7]);
          result = bred(result+r1+r2+r3+r4+r5+r6+r7);

          //result = bred(result+f[j]*(int32_t)g[i-j]);
      }
      fg[i] = result;
    }


    for (i = 2 * NTRUP_P - 2; i >= NTRUP_P; i--) {
      fg[i - NTRUP_P] = (int16_t)bred(fg[i - NTRUP_P] + fg[i]);
	  fg[i - NTRUP_P + 1] = (int16_t)bred(fg[i - NTRUP_P + 1] + fg[i]);
    }

    for (i = 0; i < NTRUP_P; ++i) 
        h[i] = fg[i];
}
