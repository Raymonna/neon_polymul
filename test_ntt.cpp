#include <bits/stdc++.h>

using namespace std;

typedef long long ll;
typedef double db;

const int inf=0x3f3f3f3f;

int getint()
{
    int f=1,g=0;char c=getchar();
    while(c<'0' || c>'9'){if(c=='-')f=-1;c=getchar();}
    while(c>='0' && c<='9')g=(g<<3)+(g<<1)+c-'0',c=getchar();
    return f*g;
}

const int maxn=300005;
const int mod=998244353;

const int G=3;

int a[maxn];
int b[maxn];
int c[maxn];
int n,m;
int rev[maxn];
int N;
int len;
int inv;

int power(ll x,ll y)
{
    ll res=1ll;
    for(;y;y>>=1,x=(x*x)%mod)
    {
        if(y&1)res=(res*x)%mod;
    }
    return res;
}

void init()
{
    while((n+m)>=(1<<len))len++;
    N=(1<<len);
    inv=power(N,mod-2);
    for(int i=0;i<N;i++)
    {
        int pos=0;
        int temp=i;
        for(int j=1;j<=len;j++)
        {
            pos<<=1;pos |= temp&1;temp>>=1;
        }
        rev[i]=pos;
    }
}

void ntt(int *a,int n,int re)
{
    for(int i=0;i<n;i++)
    {
        if(rev[i]>i)
        {
            swap(a[i],a[rev[i]]);
        }
    }
    for(int i=2;i<=n;i<<=1)
    {
        int mid=i>>1;

        int wn=power(G,(mod-1)/i);
        if(re) wn=power(wn,(mod-2));
        for(int j=0;j<n;j+=i)
        {
            int w=1;
            for(int k=0;k<mid;k++)
            {
                int temp1=a[j+k];
                int temp2=(ll)a[j+k+mid]*w%mod;
                a[j+k]=(temp1+temp2);if(a[j+k]>=mod)a[j+k]-=mod;
                a[j+k+mid]=(temp1-temp2);if(a[j+k+mid]<0)a[j+k+mid]+=mod;
                w=(ll)w*wn%mod;
            }
        }
    }
    if(re)
    {
        for(int i=0;i<n;i++)
        {
            a[i]=(ll)a[i]*inv%mod;
        }
    }
}

int main()
{
    n=761;
    m=761;

    for(int i=0;i<=n;i++)
    {
        a[i]=i+1;
    }
    for(int i=0;i<=m;i++)
    {
        b[i]=i+1;     
    }

    init();

    ntt(a,N,0);
    ntt(b,N,0);
    for(int i=0;i<=N;i++)
    {
        c[i]=(ll)a[i]*b[i]%mod;
    }
    ntt(c,N,1);
    for(int i=0;i<=n+m;i++)
    {
        printf("%d%c",c[i]%4591," \n"[i==n+m]%4591);
    }

    return 0;
}
