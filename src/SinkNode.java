import java.io.FileInputStream;
import java.util.Arrays;
import java.util.Scanner;

import static java.lang.Math.abs;
import static java.lang.Math.sqrt;


public class SinkNode {
    private double[] measure;
    private double[] recSignal;
    private double[][] A;     //n行m列
    private int m;
    private int n;
    private int k;
    private final int alfa=2;
    private double[] theta;
    private double RRE;

    public static void main(String[] args) {
        SinkNode sinkNode = new SinkNode(5);
        System.out.print("恢复数据为：");
        double[][] psi = DCT_1D.getPsi(sinkNode.n);
        double[][] psiT = MatrixUtil.transpose(psi);
//        double[] recover = MatrixUtil.matMultiVec(psiT, sinkNode.recSignal, sinkNode.n, sinkNode.n);
        double[] recover = MatrixUtil.multiply(psiT, sinkNode.recSignal);
        SourceNode.testMatrix(recover);

//        double averRec = average(recover);

        double[] data= SourceNode.readData();
        double sum=0;

//        for (int i = 0; i < sinkNode.n; i++) {
//            sum += Math.abs(data[i] - recover[i])/Math.abs(data[i]);
//            System.out.print(Math.abs(data[i] - recover[i])+"  ");
//        }

        double[] diff = new double[sinkNode.n];
        for (int i = 0; i < sinkNode.n; i++) {
            diff[i] = data[i] - recover[i];
        }
        double rRE = norm(diff, sinkNode.n) / norm(data, sinkNode.n);

        System.out.print("相对恢复误差：");
        System.out.println(rRE);
    }

    private static double average(double[] recover) {

        return 0;
    }

//    public SinkNode() {
//        SourceNode sourceNode = new SourceNode();
//
//
//    }

    public SinkNode(double s) {
        SourceNode sourceNode = new SourceNode(s);

        m = sourceNode.getM();
        n = sourceNode.getN();
        k = sourceNode.getK();
        measure = sourceNode.getResult();
//        theta=sourceNode.getData();
//        double[][]  A = readA(600, 600);
//        measurementMatrix = sourceNode.readMeasurementMatrix(n, m);
//        measurementMatrix = new double[n][m];
//        double[][] A = sourceNode.getA();
//        subMat(A, measurementMatrix, m, n);
//        soudrceNode.testMatrix(measurementMatrix);
        A=sourceNode.getA();
        double[][] subA=new double[n][m];
        subMat(A,subA,m,n);

        int supportNum = k + alfa * k;
        int[] support = new int[m];
        int[] hIndex = new int[n];
        int[] measurementMatrixIndex = new int[n];
        int[] measurementMatrixIndex2 = new int[k];
        double[] rError = Arrays.copyOf(measure, m); //copyVex
        double[] h = new double[n];
        double[] temp1 = new double[m];
        double[] temp2 = new double[m];
        recSignal = new double[n];      //rec=0
        double eps=0.000001;    //求伪逆的误差上限
        int ka=n>m?n+1:m+1;
        int num=0;
        int j=0,i,result=0;
        while (j++<40) {
            matMultiVec(subA, rError,h,n,m);
//            sourceNode.testMatrix(h);

//            Arrays.sort(h);
            sort(h,hIndex,n);
//            for (int i : hIndex) {
//                System.out.print(h[i]+"  ");
//            }
            if (num++ == 0) {
                for (i = 0; i < alfa * k; i++) {
                    support[i] = hIndex[i];
                }
                supportNum = alfa * k;
                support=bubbleSort(support,alfa*k);
            }else {
                for (i = 0; i < k; i++) {
                    support[i] = measurementMatrixIndex2[i];
                }
                for (i = 0; i < alfa * k; i++) {
                    support[k + i] = hIndex[i];     //s>3 todo
                }
                supportNum = k + alfa * k;
                support = bubbleSort(support, (alfa + 1) * k);
            }
            result = check(support, supportNum);       //重复个数
            supportNum = supportNum - result;
//            double[][] temp3 = new double[supportNum][m];
            double[][] temp3 = new double[m][supportNum];
            double[][] temp4 = new double[supportNum][m];
            double[][] temp5 = new double[m][m];
            double[][] temp6 = new double[supportNum][supportNum];

            matColCopyMat(A, temp3, support, supportNum,m,n);
            bginv(temp3, m, supportNum, temp4, eps, temp5, temp6, ka);
            matMultiVec(temp4, measure,temp1,supportNum,m);
//            sourceNode.testMatrix(temp1);
//            sourceNode.testMatrix(temp4);
//            sourceNode.testMatrix(temp5);
//            sourceNode.testMatrix(temp6);

            if (num == 1) {
                for (i = 0; i < supportNum; i++) {
                    recSignal[support[i]] = temp1[i];
                }
            } else {
                for (i = 0; i < supportNum; i++) {
                    recSignal[support[i]] = temp1[i];
                }
            }

//            sourceNode.testMatrix(recSignal);

            sort(recSignal, measurementMatrixIndex, n);
            for (i = 0; i < n - k; i++) {
                recSignal[measurementMatrixIndex[i + k]] = 0;
            }
            for (i = 0; i < k; i++) {
                measurementMatrixIndex2[i] = measurementMatrixIndex[i];
            }
            matMultiVec(A, recSignal, temp2, m, n);
            subVec(measure, temp2, rError, m);
        }
//        sourceNode.testMatrix(temp1);
//        sourceNode.testMatrix(temp2);
//        System.out.print("相对恢复误差为：");
//        sourceNode.testMatrix(rError);
        System.out.print("重构后的信号为：");
        sourceNode.testMatrix(recSignal);
        RRE = getRRE();
    }

    public double getRRE() {
        double[][] psi = DCT_1D.getPsi(n);
        double[][] psiT = MatrixUtil.transpose(psi);
//        double[] recover = MatrixUtil.matMultiVec(psiT, sinkNode.recSignal, sinkNode.n, sinkNode.n);
        double[] recover = MatrixUtil.multiply(psiT, recSignal);
        System.out.print("恢复数据为：");
        SourceNode.testMatrix(recover);

//        double averRec = average(recover);

        double[] data= SourceNode.readData();

//        for (int i = 0; i < sinkNode.n; i++) {
//            sum += Math.abs(data[i] - recover[i])/Math.abs(data[i]);
//            System.out.print(Math.abs(data[i] - recover[i])+"  ");
//        }

        double[] diff = new double[n];
        for (int i = 0; i < n; i++) {
            diff[i] = data[i] - recover[i];
        }
        double rRE = norm(diff, n) / norm(data,n);

        return rRE;
    }

    private void subVec(double[] op1,double[] op2,double[] result,int size)
        {
            int i;
            for(i = 0;i < size; i++)
            {
                result[i] = op1[i];
                result[i]=result[i] - op2[i];
            }
        }

    private int bginv(double[][] a,int m,int n,double[][] aa,double eps,double[][] u,double[][] v,int ka) {
        int i,j,k,l,t,p,q,f;
        i=bmuav(a,m,n,u,v,eps,ka);
        if (i<0)
            return(-1);
        j=n;
        if (m<n)
            j=m;
        j=j-1;
        k=0;
//        while ((k<=j)&&(a[k*n+k]!=0.0))
        while ((k<=j)&&(a[k][k]!=0.0))
            k=k+1;
        k=k-1;
        for (i=0; i<=n-1; i++)
            for (j=0; j<=m-1; j++)
            {
//                t=i*m+j; aa[t]=0.0;
                aa[i][j]=0.0;
                for (l=0; l<=k; l++)
                {
//                    f=l*n+i;
//                    p=j*m+l;
//                    q=l*n+l;
//                    aa[t]=aa[t]+v[f]*u[p]/a[q];
                    aa[i][j] = aa[i][j] + v[l][i] * u[j][l] / a[l][l];
                }
            }
        return(1);
    }

    int bmuav(double[][] a,int m,int n,double[][] u,double[][] v,double eps,int ka)
    {
//        int i,j,k,l,it,ll,kk,ix,iy,mm,nn,iz,m1,ks;
//        double d,dd,t,sm,sm1,em1,sk,ek,b,c,shh,fg[2],cs[2];
//        double *s,*e,*w;
//
//        s=(double*)malloc(ka*sizeof(double));
//        e=(double*)malloc(ka*sizeof(double));
//        w=(double*)malloc(ka*sizeof(double));

        int i,j,k,l,it,ll,kk,ix,iy,mm,nn,iz,m1,ks;
        double d,dd,t,sm,sm1,em1,sk,ek,b,c,shh;
        double[] fg = new double[2];
        double[] cs = new double[2];
        double[] s = new double[ka];
        double[] e = new double[ka];
        double[] w = new double[ka];


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
//                        d=d+a[ix]*a[ix];
                        d=d+a[i-1][kk-1]*a[i-1][kk-1];
                    }
                    s[kk-1]=sqrt(d);
                    if (s[kk-1]!=0.0)
                    {
                        ix=(kk-1)*n+kk-1;
//                        if (a[ix]!=0.0)
                        if (a[kk-1][kk-1]!=0.0)
                        {
                            s[kk-1]=Math.abs(s[kk-1]);
                            if (a[kk-1][kk-1]<0.0)
                                s[kk-1]=-s[kk-1];
                        }
                        for (i=kk; i<=m; i++)
                        {
                            iy=(i-1)*n+kk-1;
                            a[i-1][kk-1]=a[i-1][kk-1]/s[kk-1];
                        }
                        a[kk-1][kk-1]=1.0+a[kk-1][kk-1];
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
//                                d=d+a[ix]*a[iy];
                                d=d+a[i-1][kk-1]*a[i-1][j-1];
                            }
                            d=-d/a[(kk-1)][kk-1];
                            for (i=kk; i<=m; i++)
                            {
                                ix=(i-1)*n+j-1;
                                iy=(i-1)*n+kk-1;
                                a[i-1][j-1]=a[i-1][j-1]+d*a[i-1][kk-1];
                            }
                        }
                        e[j-1]=a[(kk-1)][j-1];
                    }
                }
                if (kk<=k)
                {
                    for (i=kk; i<=m; i++)
                    {
                        ix=(i-1)*m+kk-1; iy=(i-1)*n+kk-1;
                        u[i-1][kk-1]=a[i-1][kk-1];
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
                            e[kk-1]=Math.abs(e[kk-1]);
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
                                w[i-1]=w[i-1]+e[j-1]*a[(i-1)][j-1];
                        for (j=kk+1; j<=n; j++)
                            for (i=kk+1; i<=m; i++)
                            {
                                ix=(i-1)*n+j-1;
                                a[i-1][j-1]=a[i-1][j-1]-w[i-1]*e[j-1]/e[kk];
                            }
                    }
                    for (i=kk+1; i<=n; i++)
                        v[(i-1)][kk-1]=e[i-1];
                }
            }
        }
        mm=n;
        if (m+1<n)
            mm=m+1;
        if (k<n)
            s[k]=a[k][k];
        if (m<mm)
            s[mm-1]=0.0;
        if (l+1<mm)
            e[l]=a[l][mm-1];
        e[mm-1]=0.0;
        nn=m;
        if (m>n)
            nn=n;
        if (nn>=k+1)
        {
            for (j=k+1; j<=nn; j++)
            {
                for (i=1; i<=m; i++)
                    u[(i-1)][j-1]=0.0;
                u[(j-1)][j-1]=1.0;
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
                                d=d+u[i-1][kk-1]*u[i-1][j-1]/u[kk-1][kk-1];
                            }
                            d=-d;
                            for (i=kk; i<=m; i++)
                            {
                                ix=(i-1)*m+j-1;
                                iy=(i-1)*m+kk-1;
                                u[i-1][j-1]=u[i-1][j-1]+d*u[i-1][kk-1];
                            }
                        }
                    for (i=kk; i<=m; i++)
                    {
                        ix=(i-1)*m+kk-1;
                        u[i-1][kk-1]=-u[i-1][kk-1];
                    }
                    u[kk-1][kk-1]=1.0+u[kk-1][kk-1];
                    if (kk-1>=1)
                        for (i=1; i<=kk-1; i++)
                            u[(i-1)][kk-1]=0.0;
                }
                else
                {
                    for (i=1; i<=m; i++)
                        u[(i-1)][kk-1]=0.0;
                    u[(kk-1)][kk-1]=1.0;
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
                        d=d+v[i-1][kk-1]*v[i-1][j-1]/v[kk][kk-1];
                    }
                    d=-d;
                    for (i=kk+1; i<=n; i++)
                    {
                        ix=(i-1)*n+j-1; iy=(i-1)*n+kk-1;
                        v[i-1][j-1]=v[i-1][j-1]+d*v[i-1][kk-1];
                    }
                }
            }
            for (i=1; i<=n; i++)
                v[(i-1)][kk-1]=0.0;
//    todo        v[iz-n]=1.0;
            v[kk-1][kk-1]=1.0;
        }
        for (i=1; i<=m; i++)
            for (j=1; j<=n; j++)
                a[(i-1)][j-1]=0.0;
        m1=mm;
        it=60;
        while (1==1)
        {
            if (mm==0)
            {
                ppp(a,e,s,v,m,n);
                return(1);
            }
            if (it==0)
            {
                ppp(a,e,s,v,m,n);
                return(-1);
            }
            kk=mm-1;
            while ((kk!=0)&&(Math.abs(e[kk-1])!=0.0))
            {
                d=Math.abs(s[kk-1])+Math.abs(s[kk]);
                dd=Math.abs(e[kk-1]);
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
                        ix=(i-1)*n+kk-1;
                        v[i-1][kk-1]=-v[i-1][kk-1];
                    }
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
                            d=v[i-1][kk-1];
                            v[i-1][kk-1]=v[i-1][kk];
                            v[i-1][kk]=d;
                        }
                    if (kk<m)
                        for (i=1; i<=m; i++)
                        {
                            ix=(i-1)*m+kk-1; iy=(i-1)*m+kk;
                            d=u[i-1][kk-1];
                            u[i-1][kk-1]=u[i-1][kk];
                            u[i-1][kk]=d;
                        }
                    kk=kk+1;
                }
                it=60;
                mm=mm-1;
            }
            else
            {
                ks=mm;
                while ((ks>kk)&&(Math.abs(s[ks-1])!=0.0))
                {
                    d=0.0;
                    if (ks!=mm)
                        d=d+Math.abs(e[ks-1]);
                    if (ks!=kk+1)
                        d=d+Math.abs(e[ks-2]);
                    dd=Math.abs(s[ks-1]);
                    if (dd>eps*d)
                        ks=ks-1;
                    else
                        s[ks-1]=0.0;
                }
                if (ks==kk)
                {
                    kk=kk+1;
                    d=Math.abs(s[mm-1]);
                    t=Math.abs(s[mm-2]);
                    if (t>d)
                        d=t;
                    t=Math.abs(e[mm-2]);
                    if (t>d)
                        d=t;
                    t=Math.abs(s[kk-1]);
                    if (t>d)
                        d=t;
                    t=Math.abs(e[kk-1]);
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
                                d=cs[0]*v[j-1][i-1]+cs[1]*v[j-1][i];
                                v[j-1][i]=-cs[1]*v[j-1][i-1]+cs[0]*v[j-1][i];
                                v[j-1][i-1]=d;
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
                                    d=cs[0]*u[j-1][i-1]+cs[1]*u[j-1][i];
                                    u[j-1][i]=-cs[1]*u[j-1][i-1]+cs[0]*u[j-1][i];
                                    u[j-1][i-1]=d;
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
                                    d=cs[0]*v[j-1][i-1]+cs[1]*v[j-1][mm-1];
                                    v[j-1][mm-1]=-cs[1]*v[j-1][i-1]+cs[0]*v[j-1][mm-1];
                                    v[j-1][i-1]=d;
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
                                    d=cs[0]*u[j-1][i-1]+cs[1]*u[j-1][kk-2];
                                    u[j-1][kk-2]=-cs[1]*u[j-1][i-1]+cs[0]*u[j-1][kk-2];
                                    u[j-1][i-1]=d;
                                }
                        }
                    }
                }
            }
        }
//        return(1);
    }

    private void sss(double[] fg, double[] cs){
        double r, d;
        if ((abs(fg[0]) + abs(fg[1])) == 0.0) {
            cs[0] = 1.0;
            cs[1] = 0.0;
            d = 0.0;
        } else {
            d = sqrt(fg[0] * fg[0] + fg[1] * fg[1]);
            if (abs(fg[0]) > abs(fg[1])) {
                d = abs(d);
                if (fg[0] < 0.0)
                    d = -d;
            }
            if (abs(fg[1]) >= abs(fg[0])) {
                d = abs(d);
                if (fg[1] < 0.0)
                    d = -d;
            }
            cs[0] = fg[0] / d;
            cs[1] = fg[1] / d;
        }
        r = 1.0;
        if (abs(fg[0]) > abs(fg[1]))
            r = cs[1];
        else if (cs[0] != 0.0)
            r = 1.0 / cs[0];
        fg[0] = d;
        fg[1] = r;
        return;
    }


    private void ppp(double[][] a,double[] e,double[] s,double[][] v,int m,int n){
        int i,j,p,q;
        double d;
        if (m>=n)
            i=n;
        else
            i=m;
        for (j=1; j<=i-1; j++)
        {
            a[(j-1)][j-1]=s[j-1];
            a[(j-1)][j]=e[j-1];
        }
        a[(i-1)][i-1]=s[i-1];
        if (m<n)
            a[(i-1)][i]=e[i-1];
        for (i=1; i<=n-1; i++)
            for (j=i+1; j<=n; j++)
            {
                p=(i-1)*n+j-1;
                q=(j-1)*n+i-1;
                d=v[i-1][j-1];
                v[i-1][j-1]=v[j-1][i-1];
                v[j-1][i-1]=d;
            }
        return;
    }



    private void subMat(double[][] a, double[][] measurementMatrix, int m, int n) {
        int i,j;
        for(i=0;i<m;i++)
            for(j=0;j<n;j++)
            {
                measurementMatrix[j][i] = a[i][j];
            }
    }

    public void matColCopyMat(double[][] measurementMatrix, double[][] temp3, int[] support, int supportNum, int m, int n) {
        int i,j;
        for (i=0;i<m;i++)
            for (j=0;j<supportNum;j++)
            {
                temp3[i][j] = measurementMatrix[i][support[j]];
            }
    }

    private int check(int[] support, int supportNum) {
        int i,j,result;
        result = 0;
        for (i=0;i<supportNum-1;i++)
        {
            if (support[i]== support[i+1])
            {
                for (j=i+1;j<supportNum-1-result;j++)
                    support[j] = support[j+1];
                support[supportNum-1-result] = -1-result;
                result++;
            }
        }
        return result;
    }

    public double getrRE() {
        double[][] psi = DCT_1D.getPsi(n);
        double[][] psiT = MatrixUtil.transpose(psi);
//        double[] recover = MatrixUtil.matMultiVec(psiT, sinkNode.recSignal, sinkNode.n, sinkNode.n);
        double[] recover = MatrixUtil.multiply(psiT, recSignal);
//        double rRE = process(RRE);
//        double averRec = average(recover);

        double[] data= SourceNode.readData();

//        for (int i = 0; i < sinkNode.n; i++) {
//            sum += Math.abs(data[i] - recover[i])/Math.abs(data[i]);
//            System.out.print(Math.abs(data[i] - recover[i])+"  ");
//        }
        return RRE;
    }

    private int[] bubbleSort(int[] support, int num) {
        int i,j,temp;
        for (i=0;i<num-1;i++)
        {
            for (j=num-1;j>i;j--)
            {
                if (support[j-1]>support[j])
                {
                    temp = support[j-1];
                    support[j-1] = support[j];
                    support[j] = temp;
                }
            }
        }
        return support;
    }

    private void sort(double[] a,int[] index,int n)
    {
        int i,j,x;
        int[] flag = new int[n];
        double max,temp;
        double[] a_new = new double[n];

        for (i=0;i<n;i++)
        {
            a_new[i] = Math.abs(a[i]);
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
    }



    private void matMultiVec(double[][] A, double[] x, double[] y, int row, int col)
    {
        int i,j;
        for(i=0;i<row;i++)
        {
            y[i] = 0;
            for(j=0;j<col;j++)
            {
                y[i] =  y[i] + (A[i][j] * x[j]) ;
            }
        }
    }

    public static double[][] readA(int m, int length) {
        double[][] matrix = new double[m][length];
        try (Scanner sc=new Scanner(new FileInputStream("requirement/phi.txt"))) {
            for(int i=0;i<m;i++) {
                for(int j=0;j<length;j++) {
                    if(sc.hasNext())
                        matrix[i][j] = Double.parseDouble(sc.next());
                }
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
        return matrix;
    }

    static double norm(double[] x, int width)
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
}
