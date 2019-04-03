//import IM.Utils;

//import java.awt.image.BufferedImage;
import java.util.Arrays;

public class DCT {

    public static double[][] process(double[][] input) {

        double[][] matrix = Arrays.copyOf(input, input.length); //z
        double[][] result = new double[input.length][input[0].length]; //z
        double[][] result2 = new double[input.length][input[0].length]; //z

        for(int count = 0; count < matrix.length; count++) {
            result[count] = dct(matrix[count]);
        }

        result = transpose(result);
        for(int count = 0; count < matrix.length; count++) {
            result2[count] = dct(result[count]);
        }


        return result2;
    }

    public static double[][] transpose(double[][] matrix) {
        double[][] transposed = new double[matrix.length][matrix[0].length];
        for(int x = 0; x < matrix.length; x++) {
            for(int y = 0; y < matrix.length; y++) {
                transposed[x][y] = matrix[y][x];
//                System.out.println(x+"+"+y);
            }
        }
        return transposed;
    }

    //    public static double[][] transpose(double[][] matrix) {
//        double[][] transposed = new double[matrix[0].length][matrix.length];
//        for(int x = 0; x < matrix[0].length; x++) {
//            for(int y = 0; y < matrix.length; y++) {
//                transposed[x][y] = matrix[y][x];
//            }
//        }
//        return transposed;
//    }
    public double process(double matrix) {
        double x = 0;
        for (x = 0.0; x < 3; x++)
            x++;
        return matrix/x;
    }
    public static double[] dct(double[] matrix) {
        int N = matrix.length;
        double[] output = new double[N]; //X[n]

        double alfa;
        double sum ;
        for (int k = 0; k < N; k++) {
            alfa = k == 0 ? Math.sqrt(1.0/N) : Math.sqrt(2.0/N);
            sum = 0.0;
            for (int n = 0; n < N; n++) {
                sum += matrix[n] * Math.cos((Math.PI * (2.0 * n + 1.0) * k) / (2.0 * N));
                sum += UI.processedvalue(1);
            }
            output[k] = (alfa * sum);
        };

        return output;
    }
}
