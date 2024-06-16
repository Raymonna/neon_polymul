
#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include "params.h"

#include "poly.h"

//=============== Begin ===============
typedef int64_t ll;
#define Max 500100
const ll Mod = 677795329;   // 模数
const ll SpMul = 34286;  // 原根
ll twiddles[7] = { 677795328, 65706969, 505652000, 553072407, 570556391, 13494887, 467151492};

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
        //ll wn = qpow(SpMul, (Mod - 1) / h);
		ll wn = twiddles[tw_i++];
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

ll data1_1[Max], data1_2[Max];
ll data1_3[Max], data1_4[Max];
ll data1_5[Max], data1_6[Max];

ll data2_1[Max], data2_2[Max];
ll data2_3[Max], data2_4[Max];
ll data2_5[Max], data2_6[Max];
//===============  End ===============
typedef long long ll;
const int mod=NTRUP_Q,g=3;




inline int qpow(int x,int k)
{
    int ans=1;
    while(k)
    {
        if(k&1)
            ans=(ll)ans*x%mod;
        x=(ll)x*x%mod,k>>=1;
    }
    return ans;
}

inline int module(int x,int y)
{
    x+=y;
    if(x>=mod)
        x-=mod;
    return x;
}

int rev[4*NTRUP_P];
inline void NTT(int*t,int lim,int type)
{
    
    for(int i=0;i<lim;++i)
        if(i<rev[i])
            swap(&t[i],&t[rev[i]]);
    
    for(int i=1;i<lim;i<<=1)
    {
        int gn=qpow(g,(mod-1)/(i<<1));
        if(type==-1)
            gn=qpow(gn,mod-2);
        for(int j=0;j<lim;j+=(i<<1))
        {
            int gi=1;
            for(int k=0;k<i;++k,gi=(ll)gi*gn%mod)
            {
                int x=t[j+k],y=(ll)gi*t[j+i+k]%mod;
                t[j+k]=module(x,y);
                t[j+i+k]=module(x,mod-y);
            }
        }
    }
    if(type==-1)
    {
        int inv=qpow(lim,mod-2);
        for(rg int i=0;i<lim;++i)
            t[i]=(ll)t[i]*inv%mod;
    }
    
}

int X[4*NTRUP_P],Y[4*NTRUP_P];
inline void mul(int*x, int*y, int n, int m)
{
    
    memset(X,0,sizeof(X));
    memset(Y,0,sizeof(Y));
    int lim = 1, L = 0;  //L=0必须写，局部变量默认值很可能不是0
    while(lim <= n + m) {lim <<= 1; L++;}   //lim为大于(n+m)的2的幂，所以最多需要4倍空间
    printf("in mul, lim = %d\n", lim);
    for(int i = 0; i < lim; i++) rev[i] = (rev[i >> 1] >> 1) | ((i & 1) << (L - 1));
    
    for(rg int i=0;i<lim;++i) {X[i]=x[i];Y[i]=y[i];}
    
    NTT(X,lim,1);
    NTT(Y,lim,1);
    
    for(int i=0;i<lim;++i) X[i]=(ll)X[i]*Y[i]%mod;
    NTT(X,lim,-1);
    for(int i=0;i<lim;++i) x[i]=X[i];
}



int a[4*NTRUP_P], b[4*NTRUP_P];


static int16_t bred(int32_t a, int16_t m, uint32_t mu) {
    int32_t q = ((uint64_t)a * mu) >> 16;
    int16_t t = a - q * m;
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
    /*
    int16_t fg[NTRUP_P + NTRUP_P - 1];
    int16_t result;
    int i,j;
	uint32_t mu = (1UL << 16) / NTRUP_Q;
    for (i = 0; i < NTRUP_P; i++) {
      result = 0;
      for (j = 0;j <= i;++j) 
          result = bred(result+f[j]*(int32_t)g[i-j], NTRUP_Q, mu);
      fg[i] = result;
    }
    for (i = NTRUP_P;i < 2 * NTRUP_P - 1;++i) {
      result = 0;
      for (j = i - NTRUP_P + 1; j < NTRUP_P; ++j) 
          result = bred(result+f[j]*(int32_t)g[i-j], NTRUP_Q, mu);
      fg[i] = result;
    }

    for (i = 2 * NTRUP_P - 2; i >= NTRUP_P; i--) {
      fg[i - NTRUP_P] = bred(fg[i - NTRUP_P] + fg[i], NTRUP_Q, mu);
      fg[i - NTRUP_P + 1] =bred(fg[i - NTRUP_P + 1] + fg[i], NTRUP_Q, mu);
    }
	*/
	for(int i = 0; i < NTRUP_P; i++){
		a[i] = f[i];
		b[i] = g[i];
	}
	mul(a, b, NTRUP_P, NTRUP_P);
	for(int i = 0; i < 2*NTRUP_P - 1; i++){
		h[i] = a[i];
	}
	
	uint32_t mu = (1UL << 16) / NTRUP_Q;
	for (int i = 2 * NTRUP_P - 2; i >= NTRUP_P; i--) {
      h[i - NTRUP_P] = bred(h[i - NTRUP_P] + h[i], NTRUP_Q, mu);
      h[i - NTRUP_P + 1] =bred(h[i - NTRUP_P + 1] + h[i], NTRUP_Q, mu);
    }
	
}
























