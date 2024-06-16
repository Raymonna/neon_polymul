
#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include "params.h"

#include "poly.h"

//=============== Begin ===============
typedef int64_t ll;
#define Max 256
const ll Mod = 689522689;//998244353;//689522689;//677795329;   // 模数
const ll SpMul = 45903;//3;//34286;//45903;//34286;  // 原根
ll twiddles[7] = { 1, 1, 1, 689522688, 88277704, 38177824, 363540571, 523078855};

void swap(ll *a, ll *b)
{
	ll tmp = *a;
	*a = *b;
	*b = tmp;
}

ll qpow(ll a, ll k) {
    ll res = 1;
    while (k > 0) {
        if (k & 1)res = res * a % Mod;
        a = a * a % Mod;
        k >>= 1;
    }
    return res;
}

void Change(ll y[], int len) {
    for (int i = 1, j = len / 2; i < len - 1; i++) {
        if (i < j)swap(&y[i], &y[j]);
        int k = len / 2;
        while (j >= k) {
            j -= k;
            k /= 2;
        }
        if (j < k)j += k;
    }
}


void NTT(ll y[], int len, int on) {
    Change(y, len);
    for (int h = 2, tw_i = 0; h <= len; h <<= 1) {
        ll wn = qpow(SpMul, (Mod - 1) / h);
		//ll wn = twiddles[tw_i++];
		if (on == -1)wn = qpow(wn, Mod - 2);
        for (int j = 0; j < len; j += h) {
            ll w = 1LL;
            for (int k = j; k < j + h / 2; k++) {
                ll u = y[k];
                ll t = (w * y[k + h / 2]) % Mod;
                y[k] = (u + t) % Mod;
                y[k + h / 2] = (u - t + Mod) % Mod;
                w = w * wn % Mod;
            }
        }
    }
    if (on == -1) {
        ll t = qpow(len, Mod - 2);
        for (int i = 0; i < len; i++)
            y[i] = (y[i] * t) % Mod;
    }
}

ll data1[6][Max], data2[6][Max];
ll tmp[36][Max];
//===============  End ===============



static int64_t bred(int64_t a, int32_t m, uint64_t mu) {
    int64_t q = ((uint64_t)a * mu) >> 32;
    int64_t t = a - q * m;
    if (t > (m >> 1)) {
        t -= m;
    }
    if (t < -(m >> 1)) {
        t += m;
    }
    return t;
}

static int16_t bred_smallm(int64_t a, int16_t m, uint32_t mu) {
    int64_t q = ((uint64_t)a * mu) >> 16;
    int64_t t = a - q * m;
    if (t > (m >> 1)) {
        t -= m;
    }
    if (t < -(m >> 1)) {
        t += m;
    }
    return t;
}

static int16_t cmod(int32_t a, int16_t m){
    int16_t t;
    t = a % m;
    if(t > (m >> 1)){
        t -= m;
    }
    if(t < -(m >> 1)){
        t += m;
    }
    return t;
}


void poly_Rq_mul_small(int16_t *h, const int16_t *f,const int8_t *g)
{
	
	uint32_t mu = (1UL << 16) / NTRUP_Q;
	uint64_t mu_big = (1UL << 32) / Mod;
	//1. Initialize
	int16_t len = 128;
	
	for(int16_t d = 0; d < 5; d++){
		memset(data1[d], 0, sizeof(data1[0]));
		for(int16_t i = 0; i < len; i++){
			data1[d][i] = f[d*len + i];
		}
	}
	for(int i = 0; i < 10; i++){
		printf("data1[0][%d] = %ld\n", i, data1[0][i]);
	}
	printf("============\n");
	memset(data1[5], 0, sizeof(data1[0]));
	for(int16_t i = 0; i < 121; i++){
		data1[5][i] = f[640+i];
	}
	
	for(int16_t d = 0; d < 5; d++){
		memset(data2[d], 0, sizeof(data2[0]));
		for(int16_t i = 0; i < len; i++){
			data2[d][i] = g[d*len + i];
		}
	}
	memset(data2[5], 0, sizeof(data2[0]));
	for(int16_t i = 0; i < 121; i++){
		data2[5][i] = g[640+i];
	}
	int16_t len256 = 256;
	
	//2. NTT
    for(int16_t i = 0; i < 6; i++){
    	NTT(data1[i], len256, 1);
    }
	for(int16_t i = 0; i < 6; i++){
    	NTT(data2[i], len256, 1);
    }
    for(int i = 0; i < 10; i++){
		printf("data1[0][%d] = %ld\n", i, data1[0][i]);
	}

	//3. element-wise product
	for(int16_t i = 0; i < 6; i++){
		for(int16_t j = 0; j < 6; j++){
			for(int16_t e = 0; e < len256; e++){
				tmp[j + 6*i][e] = (data1[i][e] * data2[j][e]) % Mod;//bred(data1[i][e] * data2[j][e] , Mod, mu_big);
			}
		}
	}
	
	//4. INTT
	for(int16_t i = 0; i < 36; i++){
		NTT(tmp[i], len256, -1);
	}
	int32_t tmp_h[(len<<3) + (len<<2)];
	//5. Construct h from tmp
	for(int16_t i = 0; i < 6; i++){
		for(int16_t j = 0; j < 6; j++){
			for(int16_t e = 0; e < len256; e++){
				tmp_h[((i + j)<<7) + e] = bred_smallm(tmp[(i<<2) + (i<<1) + j][e], NTRUP_Q, mu); //tmp[(i<<2) + (i<<1) + j][e] % NTRUP_Q;//bred(tmp[(i<<2) + (i<<1) + j][e], NTRUP_Q, mu); 
			}
		}	
	}
	
	for(int i = 0; i < 10; i++){
		printf("tmp_h[%d]= %ld\n", i, tmp_h[i]);
	}
	
	for (int i = 2 * NTRUP_P - 2; i >= NTRUP_P; i--) {
      h[i - NTRUP_P] = bred_smallm(tmp_h[i - NTRUP_P] + tmp_h[i], NTRUP_Q, mu);//(tmp_h[i - NTRUP_P] + tmp_h[i]) % NTRUP_Q;//bred_smallm(tmp_h[i - NTRUP_P] + tmp_h[i], NTRUP_Q, mu);
      h[i - NTRUP_P + 1] = (tmp_h[i - NTRUP_P + 1] + tmp_h[i]) % NTRUP_Q;//bred_smallm(tmp_h[i - NTRUP_P + 1] + tmp_h[i], NTRUP_Q, mu);
    }
	
}
























