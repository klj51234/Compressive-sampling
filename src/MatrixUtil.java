public class MatrixUtil {
    public static double[] matMultiVec(double[][] A, double[] x, int row, int col)
    {
        double[] y=new double[row];
        int i,j;
        for(i=0;i<row;i++)
        {
            y[i] = 0;
            for(j=0;j<col;j++)
            {
                y[i] =  y[i] + (A[i][j] * x[j]) ;
            }
        }
        return y;
    }

    public static double[][] transpose(double[][] psi) {
        int row = psi[0].length;
        int col = psi.length;
        double[][] trans = new double[row][col];
        for (int i = 0; i < row; i++) {
            for (int j = 0; j < col; j++) {
                trans[i][j] = psi[j][i];
            }
        }
        return trans;
    }

    public static double[] multiply(double[][] measurementMatrix, double[] theta) {
        int length = measurementMatrix.length;
        double[] result = new double[length];
        int n = theta.length;
        double sum;
        for (int i=0;i<length;i++) {
            sum=0.0;
            for(int j=0;j<n;j++) {
                sum += measurementMatrix[i][j] * theta[j];
            }
            result[i]=sum;
        }
        return result;
    }

}
