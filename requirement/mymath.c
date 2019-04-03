#include"mymath.h"
#include"math.h"
#include"stdlib.h"
#include"stdio.h"



void ppp(double* a,double* e,double* s,double* v,int m,int n)
{
	int i,j,p,q;
    double d;
    if (m>=n)
		i=n;
    else
		i=m;
    for (j=1; j<=i-1; j++)
	{ 
		a[(j-1)*n+j-1]=s[j-1];
        a[(j-1)*n+j]=e[j-1];
	}
    a[(i-1)*n+i-1]=s[i-1];
    if (m<n) 
		a[(i-1)*n+i]=e[i-1];
    for (i=1; i<=n-1; i++)
		for (j=i+1; j<=n; j++)
		{ 
			p=(i-1)*n+j-1; 
			q=(j-1)*n+i-1;
			d=v[p];
			v[p]=v[q];
			v[q]=d;
		}
    return;
  }


 void sss(double fg[2],double cs[2])
{
	double r,d;
    if ((fabs(fg[0])+fabs(fg[1]))==0.0)
	{ 
		cs[0]=1.0;
		cs[1]=0.0;
		d=0.0;
	}
    else 
    { 
		d=sqrt(fg[0]*fg[0]+fg[1]*fg[1]);
        if (fabs(fg[0])>fabs(fg[1]))
        { 
			d=fabs(d);
            if (fg[0]<0.0) 
				d=-d;
         }
		if (fabs(fg[1])>=fabs(fg[0]))
        { 
			d=fabs(d);
            if (fg[1]<0.0)
				d=-d;
         }
        cs[0]=fg[0]/d;
		cs[1]=fg[1]/d;
	}
    r=1.0;
    if (fabs(fg[0])>fabs(fg[1])) 
		r=cs[1];
    else
		if (cs[0]!=0.0)
			r=1.0/cs[0];
    fg[0]=d; 
	fg[1]=r;
    return;
  }


 int bmuav(double* a,int m,int n,double* u,double* v,double eps,int ka)
{
	int i,j,k,l,it,ll,kk,ix,iy,mm,nn,iz,m1,ks;
    double d,dd,t,sm,sm1,em1,sk,ek,b,c,shh,fg[2],cs[2];
    double *s,*e,*w;

    s=(double*)malloc(ka*sizeof(double));
    e=(double*)malloc(ka*sizeof(double));
    w=(double*)malloc(ka*sizeof(double));
    it=60;
	k=n;
    if (m-1<n) 
		k=m-1;
    l=m;
    if (n-2<m) 
		l=n-2;
    if (l<0) 
		l=0;
    ll=k;
    if (l>k)
		ll=l;
    if (ll>=1)
	{ 
		for (kk=1; kk<=ll; kk++)
        {
			if (kk<=k)
            { 
				d=0.0;
                for (i=kk; i<=m; i++)
                { 
					ix=(i-1)*n+kk-1;
					d=d+a[ix]*a[ix];
				}
                s[kk-1]=sqrt(d);
                if (s[kk-1]!=0.0)
                { 
					ix=(kk-1)*n+kk-1;
                    if (a[ix]!=0.0)
                    { 
						s[kk-1]=fabs(s[kk-1]);
                        if (a[ix]<0.0) 
							s[kk-1]=-s[kk-1];
					}
                    for (i=kk; i<=m; i++)
                    { 
						iy=(i-1)*n+kk-1;
                        a[iy]=a[iy]/s[kk-1];
                     }
                    a[ix]=1.0+a[ix];
				}
                s[kk-1]=-s[kk-1];
			}
            if (n>=kk+1)
            { 
				for (j=kk+1; j<=n; j++)
                {
					if ((kk<=k)&&(s[kk-1]!=0.0))
                    { 
						d=0.0;
                        for (i=kk; i<=m; i++)
                        {
							ix=(i-1)*n+kk-1;
                            iy=(i-1)*n+j-1;
                            d=d+a[ix]*a[iy];
                        }
                        d=-d/a[(kk-1)*n+kk-1];
                        for (i=kk; i<=m; i++)
                        { 
							ix=(i-1)*n+j-1;
                            iy=(i-1)*n+kk-1;
                            a[ix]=a[ix]+d*a[iy];
                        }
                      }
					e[j-1]=a[(kk-1)*n+j-1];
                  }
              }
            if (kk<=k)
            { 
				for (i=kk; i<=m; i++)
                { 
					ix=(i-1)*m+kk-1; iy=(i-1)*n+kk-1;
                    u[ix]=a[iy];
                 }
             }
            if (kk<=l)
            { 
				d=0.0;
                for (i=kk+1; i<=n; i++)
                  d=d+e[i-1]*e[i-1];
                e[kk-1]=sqrt(d);
                if (e[kk-1]!=0.0)
                { 
					if (e[kk]!=0.0)
                    { 
						e[kk-1]=fabs(e[kk-1]);
                        if (e[kk]<0.0)
							e[kk-1]=-e[kk-1];
                    }
                    for (i=kk+1; i<=n; i++)
						e[i-1]=e[i-1]/e[kk-1];
                    e[kk]=1.0+e[kk];
                 }
                e[kk-1]=-e[kk-1];
                if ((kk+1<=m)&&(e[kk-1]!=0.0))
                { 
					for (i=kk+1; i<=m; i++) 
						w[i-1]=0.0;
                    for (j=kk+1; j<=n; j++)
						for (i=kk+1; i<=m; i++)
                        w[i-1]=w[i-1]+e[j-1]*a[(i-1)*n+j-1];
                    for (j=kk+1; j<=n; j++)
                        for (i=kk+1; i<=m; i++)
                        { 
							ix=(i-1)*n+j-1;
							a[ix]=a[ix]-w[i-1]*e[j-1]/e[kk];
                        }
                 }
                for (i=kk+1; i<=n; i++)
                  v[(i-1)*n+kk-1]=e[i-1];
			}
		}
	}
    mm=n;
    if (m+1<n)
		mm=m+1;
    if (k<n)
		s[k]=a[k*n+k];
    if (m<mm)
		s[mm-1]=0.0;
    if (l+1<mm)
		e[l]=a[l*n+mm-1];
    e[mm-1]=0.0;
    nn=m;
    if (m>n) 
		nn=n;
    if (nn>=k+1)
    { 
		for (j=k+1; j<=nn; j++)
        { 
			for (i=1; i<=m; i++)
				u[(i-1)*m+j-1]=0.0;
				u[(j-1)*m+j-1]=1.0;
         }
    }
	if (k>=1)
    { 
		for (ll=1; ll<=k; ll++)
        { 
			kk=k-ll+1;
			iz=(kk-1)*m+kk-1;
            if (s[kk-1]!=0.0)
            { 
				if (nn>=kk+1)
                  for (j=kk+1; j<=nn; j++)
                  { 
					  d=0.0;
                      for (i=kk; i<=m; i++)
                      {
						  ix=(i-1)*m+kk-1;
                          iy=(i-1)*m+j-1;
                          d=d+u[ix]*u[iy]/u[iz];
                       }
                      d=-d;
                      for (i=kk; i<=m; i++)
                      { 
						  ix=(i-1)*m+j-1;
                          iy=(i-1)*m+kk-1;
                          u[ix]=u[ix]+d*u[iy];
                       }
                  }
                  for (i=kk; i<=m; i++)
                  { 
					  ix=(i-1)*m+kk-1; 
					  u[ix]=-u[ix];
				  }
				  u[iz]=1.0+u[iz];
                  if (kk-1>=1)
					  for (i=1; i<=kk-1; i++)
						  u[(i-1)*m+kk-1]=0.0;
			}
            else
            { 
				for (i=1; i<=m; i++)
					u[(i-1)*m+kk-1]=0.0;
                u[(kk-1)*m+kk-1]=1.0;
             }
		}
	}
    for (ll=1; ll<=n; ll++)
    { 
		kk=n-ll+1;
		iz=kk*n+kk-1;
        if ((kk<=l)&&(e[kk-1]!=0.0))
        { 
			for (j=kk+1; j<=n; j++)
            { 
				d=0.0;
                for (i=kk+1; i<=n; i++)
                { 
					ix=(i-1)*n+kk-1; iy=(i-1)*n+j-1;
                    d=d+v[ix]*v[iy]/v[iz];
                 }
                d=-d;
                for (i=kk+1; i<=n; i++)
                { 
					ix=(i-1)*n+j-1; iy=(i-1)*n+kk-1;
                    v[ix]=v[ix]+d*v[iy];
                }
            }
        }
        for (i=1; i<=n; i++)
			v[(i-1)*n+kk-1]=0.0;
        v[iz-n]=1.0;
	}
    for (i=1; i<=m; i++)
		for (j=1; j<=n; j++)
			a[(i-1)*n+j-1]=0.0;
    m1=mm;
	it=60;
    while (1==1)
    {
		if (mm==0)
		{ 
			ppp(a,e,s,v,m,n);
            free(s); free(e); free(w); return(1);
        }
        if (it==0)
        { 
			ppp(a,e,s,v,m,n);
            free(s); free(e); free(w); return(-1);
        }
        kk=mm-1;
	while ((kk!=0)&&(fabs(e[kk-1])!=0.0))
    { 
		d=fabs(s[kk-1])+fabs(s[kk]);
        dd=fabs(e[kk-1]);
        if (dd>eps*d)
			kk=kk-1;
        else 
			e[kk-1]=0.0;
     }
     if (kk==mm-1)
     {
		 kk=kk+1;
         if (s[kk-1]<0.0)
         { 
			 s[kk-1]=-s[kk-1];
             for (i=1; i<=n; i++)
             {
				 ix=(i-1)*n+kk-1; v[ix]=-v[ix];}
             }
		 while ((kk!=m1)&&(s[kk-1]<s[kk]))
         { 
			 d=s[kk-1];
			 s[kk-1]=s[kk];
			 s[kk]=d;
             if (kk<n)
				 for (i=1; i<=n; i++)
                 { 
					 ix=(i-1)*n+kk-1; iy=(i-1)*n+kk;
                     d=v[ix]; v[ix]=v[iy]; v[iy]=d;
                  }
              if (kk<m)
				  for (i=1; i<=m; i++)
                  { 
					  ix=(i-1)*m+kk-1; iy=(i-1)*m+kk;
                      d=u[ix]; u[ix]=u[iy]; u[iy]=d;
                  }
				  kk=kk+1;
		 }
         it=60;
          mm=mm-1;
         }
        else
        { 
			ks=mm;
            while ((ks>kk)&&(fabs(s[ks-1])!=0.0))
            { 
				d=0.0;
                if (ks!=mm)
					d=d+fabs(e[ks-1]);
                if (ks!=kk+1)
					d=d+fabs(e[ks-2]);
                dd=fabs(s[ks-1]);
                if (dd>eps*d)
					ks=ks-1;
                else
					s[ks-1]=0.0;
            }
            if (ks==kk)
            { 
				kk=kk+1;
                d=fabs(s[mm-1]);
                t=fabs(s[mm-2]);
                if (t>d)
					d=t;
                t=fabs(e[mm-2]);
                if (t>d)
					d=t;
                t=fabs(s[kk-1]);
                if (t>d)
					d=t;
                t=fabs(e[kk-1]);
                if (t>d)
					d=t;
                sm=s[mm-1]/d; 
				sm1=s[mm-2]/d;
                em1=e[mm-2]/d;
                sk=s[kk-1]/d; 
				ek=e[kk-1]/d;
                b=((sm1+sm)*(sm1-sm)+em1*em1)/2.0;
                c=sm*em1;
				c=c*c;
				shh=0.0;
                if ((b!=0.0)||(c!=0.0))
                { 
					shh=sqrt(b*b+c);
                    if (b<0.0) 
						shh=-shh;
                    shh=c/(b+shh);
                }
                fg[0]=(sk+sm)*(sk-sm)-shh;
                fg[1]=sk*ek;
                for (i=kk; i<=mm-1; i++)
                { 
					sss(fg,cs);
                    if (i!=kk)
						e[i-2]=fg[0];
                    fg[0]=cs[0]*s[i-1]+cs[1]*e[i-1];
                    e[i-1]=cs[0]*e[i-1]-cs[1]*s[i-1];
                    fg[1]=cs[1]*s[i];
                    s[i]=cs[0]*s[i];
                    if ((cs[0]!=1.0)||(cs[1]!=0.0))
						for (j=1; j<=n; j++)
                        { 
							ix=(j-1)*n+i-1;
							iy=(j-1)*n+i;
							d=cs[0]*v[ix]+cs[1]*v[iy];
							 v[iy]=-cs[1]*v[ix]+cs[0]*v[iy];
							 v[ix]=d;
                        }
						sss(fg,cs);
						s[i-1]=fg[0];
						fg[0]=cs[0]*e[i-1]+cs[1]*s[i];
						s[i]=-cs[1]*e[i-1]+cs[0]*s[i];
						fg[1]=cs[1]*e[i];
						e[i]=cs[0]*e[i];
						if (i<m)
							if ((cs[0]!=1.0)||(cs[1]!=0.0))
								for (j=1; j<=m; j++)
								{ 
									ix=(j-1)*m+i-1;
									iy=(j-1)*m+i;
									d=cs[0]*u[ix]+cs[1]*u[iy];
									u[iy]=-cs[1]*u[ix]+cs[0]*u[iy];
									u[ix]=d;
								}
				}
                e[mm-2]=fg[0];
                it=it-1;
			}
            else
            { 
				if (ks==mm)
                { 
					kk=kk+1;
                    fg[1]=e[mm-2];
					e[mm-2]=0.0;
                    for (ll=kk; ll<=mm-1; ll++)
                    { 
						i=mm+kk-ll-1;
                        fg[0]=s[i-1];
                        sss(fg,cs);
                        s[i-1]=fg[0];
                        if (i!=kk)
                        { 
							fg[1]=-cs[1]*e[i-2];
                            e[i-2]=cs[0]*e[i-2];
                        }
                        if ((cs[0]!=1.0)||(cs[1]!=0.0))
							for (j=1; j<=n; j++)
                            { 
								ix=(j-1)*n+i-1;
								iy=(j-1)*n+mm-1;
								d=cs[0]*v[ix]+cs[1]*v[iy];
								v[iy]=-cs[1]*v[ix]+cs[0]*v[iy];
								v[ix]=d;
                            }
					}
				}
                else
                {
					kk=ks+1;
                    fg[1]=e[kk-2];
                    e[kk-2]=0.0;
                    for (i=kk; i<=mm; i++)
                    { 
						fg[0]=s[i-1];
                        sss(fg,cs);
                        s[i-1]=fg[0];
                        fg[1]=-cs[1]*e[i-1];
                        e[i-1]=cs[0]*e[i-1];
                        if ((cs[0]!=1.0)||(cs[1]!=0.0))
							for (j=1; j<=m; j++)
                            { 
								ix=(j-1)*m+i-1;
								iy=(j-1)*m+kk-2;
								d=cs[0]*u[ix]+cs[1]*u[iy];
								u[iy]=-cs[1]*u[ix]+cs[0]*u[iy];
								u[ix]=d;
                            }
					}
				}
			}
		}
	}
	return(1);
}


 

int maxIndex(const double *x,int width)
{
	double max,temp;
	int i,index;
	max = fabs(x[0]);
	index = 0;
	for (i=1;i<width;i++)
	{
		temp = fabs(x[i]);
		if ( temp > max )
		{
			max = temp;
			index = i;
		}
	}
	return index;
}


double Norm(double *x,int width)
{
	double temp=0;
	int i;
	for(i=0;i<width;i++)
	{
		temp = temp +  x[i] * x[i];
	}
	temp = sqrt(temp);
	return  temp;
}


void copyVec(double *desc,int dInc,const double *src,int sInc,int size)
{
	int i;
	for(i = 0;i < size; i++)
	{
		memcpy(desc+i * dInc,src +i * sInc,sizeof(double));
	}
}


void subVec(const double *op1,const double *op2,double *result,int size)
{
	int i;
	for(i = 0;i < size; i++)
	{
		result[i] = op1[i];
		result[i]=result[i] - op2[i];
	}
}


void Submat(double *x, double *x_sub, int m,int n)
{
	int i,j;
	for(i=0;i<m;i++)
		for(j=0;j<n;j++)
		{
			x_sub[j*m+i] = x[i*n+j];
		}
}

void MatMultiVec(double *A, double* x, double* y, int row, int col)
{
	int i,j;
	for(i=0;i<row;i++)
	{
		y[i] = 0;
		for(j=0;j<col;j++)
		{
			y[i] =  y[i] + (A[i*col+j] * x[j]) ;
		}
	}
}


void MatColCopyMat(double *A,double *B,int *index,int num,int m,int n)
{
	int i,j;
	for (i=0;i<m;i++)
		for (j=0;j<num;j++)
		{
			B[num*i+j] = A[n*i+index[j]];
		}
}


void Sort(double *a,int *index,int n)
{
	int i,j,x,*flag;
	double max,temp,*a_new;


	flag = (int *)malloc(n * sizeof(int));
	a_new = (double *)malloc(n * sizeof(double));
	
	for (i=0;i<n;i++)
	{
		a_new[i] = fabs(a[i]);
		flag[i] = 1;
	}

	max = a_new[0];
	x = 0;
	for (i=0;i<n;i++)
	{
		for (j=0;j<n;j++)
		{
			if ( flag[j]==0 )
				continue;
			temp = a_new[j];
			if ( temp >= max )
			{
				max = temp;
				x = j;
			}
		}
		flag[x] = 0; 
		index[i] = x;
		max = 0;
	}
	free(a_new);
}


void BubbleSort(int *a,int num)
{
	int i,j,temp;
	for (i=0;i<num-1;i++)
	{
		for (j=num-1;j>i;j--)
		{
			if (a[j-1]>a[j])
			{
				temp = a[j-1];
				a[j-1] = a[j];
				a[j] = temp;
			}
		}
	}
}


int Check(int *a,int num)
{
	int i,j,result;
	result = 0;		//�ظ�����
	for (i=0;i<num-1;i++)
	{
		if (a[i]==a[i+1])
		{
			for (j=i+1;j<num-1-result;j++)
				a[j] = a[j+1];
			a[num-1-result] = -1-result;
			result++;
		}
	}
	return result;
}