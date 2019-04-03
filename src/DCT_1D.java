public class DCT_1D {
    public static double[] process(double[] data) {
        int N=data.length;
        double[] result = new double[N];
        double alfa;
        double sum;

        for (int k = 0; k < N; k++) {
            alfa = k == 0 ? Math.sqrt(1.0/N) : Math.sqrt(2.0/N);
            sum = 0.0;
            for (int n = 0; n < N; n++) {
                sum += data[n] * Math.cos((Math.PI * (2.0 * n + 1.0) * k) / (2.0 * N));
            }
            result[k] = (alfa * sum);
        }
        return result;
    }

    public static double[][] getPsi(int n) {
        double[][] result = new double[n][n];
        double alfa;
        for (int i = 0; i < n; i++) {
            alfa = i == 0 ? Math.sqrt(1.0 / n) : Math.sqrt(2.0 / n);
            for (int j = 0; j < n; j++) {
                if (i == 0) {
                    result[i][j] = alfa;
                }else {
                    result[i][j] = alfa * Math.cos(Math.PI * i * (2.0 * (j + 1.0) - 1.0) / (2.0 * n));
                }
            }
        }
        return result;
    }
}
