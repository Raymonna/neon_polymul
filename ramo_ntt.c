#include <string.h>
#include <stdio.h>
#include <stdlib.h>
typedef int64_t ll;
#define Max 500100
const ll Mod = 689522689;//677795329;   // 模数
const ll SpMul = 34286;//45903;  // 原根
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
    /*
    for (int i = 1, j = 0; i < len; i++) {
		for (int k = len >> 1; !((j ^= k) & k); k >>= 1);
		if (i < j) swap(y[i], y[j]);
	}
	*/
    printf("twiddles = { ");
    for (int h = 2, tw_i = 0; h <= len; h <<= 1) {
        ll wn = qpow(SpMul, (Mod - 1) / h);
		//ll wn = twiddles[tw_i++];
		printf("%ld, ", wn);
		if (on == -1)wn = qpow(wn, Mod - 2);
        for (int j = 0; j < len; j += h) {
            ll w = 1LL;
            for (int k = j; k < j + h / 2; k++) {
                ll u = y[k];
                ll t = (w * y[k + h / 2]) % Mod;
                y[k] = (u + t) % Mod;
                y[k + h / 2] = (u - t + Mod) % Mod;
                w = (w * wn) % Mod;
            }
        }
    }
    printf("}\n");
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
char str1[Max], str2[Max];
int main() {
    //while(~scanf("%s%s", str1, str2)){
        int len1=strlen(str1);
        int len2=strlen(str2);
        ll len=1;
        memset(data1_1, 0, sizeof(data1_1));
        memset(data1_2, 0, sizeof(data2_2));
        memset(data1_3, 0, sizeof(data2_2));
        memset(data1_4, 0, sizeof(data2_2));
        memset(data2_1, 0, sizeof(data1_2));
        memset(data2_2, 0, sizeof(data2_2));
        memset(data2_3, 0, sizeof(data2_2));
        memset(data2_4, 0, sizeof(data2_2));
 
		//while(len<len1*2||len<len2*2)len<<=1;
        len = 128;
		for(int i= 0;i < 10; i++){
			data1_1[i] = i+1;//4590/2;
			data1_2[i] = 0;//4590/2;
			data1_3[i] = 1;//4590/2;
			data1_4[i] = 0;//4590/2;
		}


		//for(int a=0; a<len1; a++)data1[a]=str1[len1-a-1]-'0';
        //for(int a=0; a<len2; a++)data2[a]=str2[len2-a-1]-'0';
        NTT(data1_1, len*2, 1);
        NTT(data1_2, len*2, 1);
        NTT(data1_3, len*2, 1);
        NTT(data1_4, len*2, 1);
		//NTT(data1_1, len, 1);
        //NTT(data1_2, len, 1);
		for(int i = 0; i < 10; i++){
			printf("data1_1[%d] = %ld\n", i, data1_1[i]);
		}

		//0~127
		for(int a=0; a<len*2; a++) data2_1[a]=((data1_1[a]%Mod)*(data1_3[a]%Mod))%Mod;
		//128~256
		for(int a=0; a<len*2; a++) data2_2[a]=data1_2[a]*data1_3[a]%Mod;
		
		for(int a=0; a<len*2; a++) data2_3[a]=data1_1[a]*data1_4[a]%Mod;
		for(int a=0; a<len*2; a++) data2_4[a]=data1_2[a]*data1_4[a]%Mod;
        //for(int a=len; a>=0; a--) printf("%d, ", data1[a]);
        NTT(data2_1, len*2, -1);
        NTT(data2_2, len*2, -1);
        NTT(data2_3, len*2, -1);
        NTT(data2_4, len*2, -1);
        
        for(int i = 0; i < 10; i++){
			printf("data2_1[%d] = %ld\n", i, data2_1[i]);
		}
        
		//for(int a=0; a<len; a++){
        //    data1[a+1]+=data1[a]/10;data1[a]%=10;
        //}
        //while(data1[len]==0&&len>0)len--;
        for(int a=10; a>=0; a--) printf("%lld, ", data2_1[a]%4591);
        printf("\n");
		for(int a=10; a>=0; a--) printf("%lld, ", (data2_2[a] + data2_3[a])%4591);
        printf("\n");
		for(int a=10; a>=0; a--) printf("%lld, ", data2_4[a]%4591);
        printf("\n");
    //}
    return 0;
}
