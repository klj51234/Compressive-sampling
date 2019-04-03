import java.io.FileInputStream;
import java.util.Arrays;
import java.util.Scanner;



public class SourceNode {
    private int m,n;      //row of the measurement matrix
    private double[] data;
    private double[][] measurementMatrix;
    private double[] theta;
    private double[] result;
    private double[][] A;
    private double[][] psi;
    private double[][] phi;
    private double[][] psiT;
    private int k;

    private final int DATA_SIZE=1000;

    public static void main(String[] args) {
        SourceNode sn = new SourceNode(5);
//        testMatrix(sn.theta);
    }
    public SourceNode(double s) {

        data=readData();
        System.out.print("原始数据为：");
        testMatrix(data);
        n = data.length;

//        theta=DCT_1D.process(data);
//        testMatrix(theta);

        psi = new double[n][n];
        psi = DCT_1D.getPsi(n);
        psiT = MatrixUtil.transpose(psi);

        theta = MatrixUtil.multiply(psi, data);
        System.out.print("离散余弦变换后：");
        testMatrix(theta);
//        k = getK(theta);
        k = getKN(theta, 40);

//        m=(int)(60*s);
        m = (int) (k * (s + 1));
//        System.out.println(m);

//        data = MatrixUtil.multiply(psiT, theta);

        /**   to test x=psiT*theta         */
//        double[] dataR = MatrixUtil.multiply(psiT, theta);
//        int err=0;
//        testMatrix(dataR);
//        for (int i = 0; i < n; i++) {
//            dataR[i] = data[i] - dataR[i];
//            if (Math.abs(dataR[i])>1) {
//                err++;
//            }
//        }
//        System.out.println(err);
//        testMatrix(dataR);


//        A = matrixMulti(psi, psiT);    证明psi为正交矩阵

        measurementMatrix=new double[m][n];
        double[][]  phi = SinkNode.readA(600,600 );
        subMat(phi,measurementMatrix,m,n);
        A = matrixMulti(measurementMatrix, psiT);
//        result = MatrixUtil.multiply(measurementMatrix, data);
        result = MatrixUtil.multiply(A, theta);
        System.out.print("传输数据y：");
        testMatrix(result);

//        double[][] AT = transpose(A);
//        testMatrix(matrixMulti(AT, A));
//        result = multiply(measurementMatrix, data);

//        this.A = DCT.process(measurementMatrix);
//        result = multiply(A, theta);

//        testMatrix(theta);
    }

    public double[] getData() {
        return data;
    }

    public double[][] matrixMulti(double[][] m1, double[][] m2) {
        int row=m1.length;
        int col = m2[0].length;
        int k = m1[0].length;
        double sum=0.0;
        double[][] result = new double[row][col];
        for (int i = 0; i < row; i++) {
            for (int j = 0; j < col; j++) {
                for (int m = 0; m < k; m++) {
                    sum += m1[i][m] * m2[m][j];
                }
                result[i][j] = sum;
                sum=0.0;
            }
        }
        return result;
    }

    public double[][] getA() {
        return A;
    }

    private int getK(double[] theta) {
        int k=0;
        for (int i = 0; i < theta.length; i++) {
            if (Math.abs(theta[i]) > 1){
                k++;
            }else {
                theta[i] = 0;
            }
        }
        return k;
    }
    private int getKN(double[] theta,int N) {
        int k=0;
//        double sum = 0.0;
        int n = theta.length;
        double[] temp = new double[n];
        for (int i = 0; i < n; i++) {
            temp[i] = Math.abs(theta[i]);
        }
        Arrays.sort(temp);
        double crit = temp[n-N];
        for (int i = 0; i < n; i++) {
            if (Math.abs(theta[i]) >= crit) {
                k++;
            }else {
//                sum += Math.abs(theta[i]);
                theta[i] = 0;
            }
        }
//        sum = sum /(n-k);
        return k;
    }

    public int getK() {
        return k;
    }

    private void subMat(double[][] a, double[][] measurementMatrix, int m, int n) {
        int i,j;
        for(i=0;i<m;i++)
            for(j=0;j<n;j++)
            {
                measurementMatrix[i][j] = a[i][j];
            }
    }


//    public double[][] readMeasurementMatrix(int m, int length) {
//        double[][] matrix = new double[m][length];
//        try (Scanner sc=new Scanner(new FileInputStream("requirement/phi.txt"))) {
//            for(int i=0;i<m;i++) {
//                for(int j=0;j<length;j++) {
//                    if(sc.hasNext())
//                        matrix[i][j] = Double.parseDouble(sc.next());
//                }
//            }
//        } catch (Exception e) {
//            e.printStackTrace();
//        }
//        return matrix;
//    }

    public static double[] readData() {
        double[] data=new double[1000];
        int i=0;
        try (Scanner sc=new Scanner(new FileInputStream("requirement/Tem_data.txt"))) {
            while (sc.hasNext()){
                data[i++] =Double.parseDouble(sc.next());
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
//        System.out.println(String.valueOf(i));
        data = Arrays.copyOf(data, i);
        return data;
    }

    public static void testMatrix(double[][] test) {
        for (double[] d : test) {
            for (double dd : d) {
                System.out.print(dd+"  ");
            }
            System.out.println();
        }
    }

    public static void testMatrix(double[] test){
        for (double d : test) {
            System.out.print(d+"  ");
        }
        System.out.println();
    }

    public static void testMatrix(int[] test){
        for (int d : test) {
            System.out.print(d+"  ");
        }
        System.out.println();
    }

    public double[] getResult() {
        return result;
    }

    public int getM() {
        return m;
    }

    public int getN() {
        return data.length;
    }

//    public double[] getTheta() {
//        return theta;
//    }
}
