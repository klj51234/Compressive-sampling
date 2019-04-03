import org.jfree.chart.ChartFactory;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.category.DefaultCategoryDataset;

import javax.swing.*;

public class UI extends MyFrame {
    public static void main(String[] args) {
        JFrame jFrame = new JFrame();
        jFrame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE);
        jFrame.setSize(800, 800);
        jFrame.setVisible(true);

        DefaultCategoryDataset dataset = new DefaultCategoryDataset();

        for (double i = 2; i <= 5; i += 0.5) {
            SinkNode sinkNode = new SinkNode(i);
            double value = sinkNode.getrRE();
            System.out.println("相对恢复误差为："+value);
            System.out.println();
            dataset.setValue((Number)value,"a",i);
        }
        JFreeChart chart = ChartFactory.createLineChart("Relative recovery error", "m/s", "RRE", dataset, PlotOrientation.VERTICAL, true, false, false);
//        try {
//            OutputStream os = new FileOutputStream("requirement/test.jpg");
//            ChartUtilities.writeChartAsJPEG(os, chart, 800, 800);
//            os.close();
//        } catch (Exception e) {
//            e.printStackTrace();
//        }

//        ImageIcon icon=new ImageIcon(chart)

        JLabel label = new JLabel();
        label.setIcon(new ImageIcon(chart.createBufferedImage(800,800)));
        label.setBounds(0,0,800,800);
        jFrame.add(label);
        jFrame.repaint();
    }

    public static double processedvalue(double value) {
        double result=0;
        return result;
    }
}
