import org.jfree.chart.ChartFactory;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.category.DefaultCategoryDataset;

import javax.swing.*;

public class MyFrame extends JFrame {
    public JFrame jFrame;
    public JLabel jLabel;
    public double index;
    public void init() {
        jFrame = new JFrame();
        jFrame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE);
        jFrame.setSize(800, 800);
        jFrame.setVisible(true);

        DefaultCategoryDataset dataset = new DefaultCategoryDataset();

        for (double i = 2; i <= 5; i += 0.5) {
            SinkNode sinkNode = new SinkNode(i);
            double value = sinkNode.getRRE();
            System.out.println(value+"+"+i);
            dataset.setValue((Number)value,"a",i);
        }
        JFreeChart chart = ChartFactory.createLineChart("Relative recovery error", "m/s", "10E-2", dataset, PlotOrientation.VERTICAL, true, false, false);
//        try {
//            OutputStream os = new FileOutputStream("requirement/test.jpg");
//            ChartUtilities.writeChartAsJPEG(os, chart, 800, 800);
//            os.close();
//        } catch (Exception e) {
//            e.printStackTrace();
//        }
//
//        ImageIcon icon=new ImageIcon(chart)
//
        jLabel = new JLabel();
        jLabel.setIcon(new ImageIcon(chart.createBufferedImage(800,800)));
        jLabel.setBounds(0,0,800,800);
        jFrame.add(jLabel);
//        jFrame.repaint();
    }

    public static double processedValue(double value) {
        double result=0;
        return result;
    }
}
