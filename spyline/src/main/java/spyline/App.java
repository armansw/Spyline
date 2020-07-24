package spyline;

import java.util.ArrayList;
import java.util.Arrays;
import java.io.FileWriter;
import java.io.IOException;
import java.io.File;
import java.util.Collections;



public class App {
    public static void main( String[] args )
    {
        ArrayList<Double> x = new ArrayList<>(Arrays.asList(1.0, 3.0, 7.0, 15.0, 31.0, 63.0, 127.0, 255.0, 383.0, 639.0, 895.0, 1151.0, 1663.0, 2175.0, 3199.0, 4223.0, 6271.0, 10367.0, 14463.0, 22655.0, 39039.0));
        ArrayList<Double> y = new ArrayList<>(Arrays.asList(900.0, 779.0, 692.0, 617.0, 551.0, 476.0, 387.0, 283.0, 255.0, 221.0, 200.0, 185.0, 164.0, 150.0, 81.0, 65.0, 46.0, 26.0, 16.0, 7.0, 1.0));
        ArrayList<Double> w = new ArrayList<>(Collections.nCopies(x.size(), 1.0));
        Points points = new Points(x, y, w);
        ArrayList<Double> xCurve = new ArrayList<>(Collections.nCopies(1000, 0.0));
        for(int i = 0; i < xCurve.size(); ++i){
            double d = (points.getItem(x.size() - 1).getX() - points.getItem(0).getX()) / (xCurve.size() - 1);
            xCurve.set(i, points.getItem(0).getX() + d*i);
        }
        ArrayList<Double> yCurve = new ArrayList<>(Collections.nCopies(xCurve.size(), 0.0));
        ArrayList<Double> yCurveOpt = new ArrayList<>(Collections.nCopies(xCurve.size(), 0.0));
        ArrayList<Double> knots = new ArrayList<>(Arrays.asList(points.getItem(0).getX(), 4.0, 8.0, 12.0, 16.0, points.getItem(x.size() - 1).getX()));
        ArrayList<Double> coefficients = new ArrayList<>(Collections.nCopies(knots.size() - 2, 0.0));
        knots.add(knots.get(knots.size() - 2) + 1);
        Collections.sort(knots);
        Spline s = new Spline(coefficients, knots, 3);
        CurveFitter c = new CurveFitter(s);
        CurveFitter.initiateGrid(s, points);
        double q = 1e-9;
        CurveFitter.approximate(s, points, q);


        int index = 0;
        for (Double point : xCurve) {
            yCurve.set(index, s.getValue(point));
            index += 1;
        }

        index = 0;
        knots = s.getKnots();
        ArrayList<Double> knotsYUni = new ArrayList<>(Collections.nCopies(knots.size(), 0.0));
        for (Double point : knots) {
            knotsYUni.set(index, s.getValue(point));
            index += 1;
        }

        c.approximateWithOptimalGrid(s, points, q, 1e-3, 1e-3);
        
        index = 0;
        for (Double point : xCurve) {
            yCurveOpt.set(index, s.getValue(point));
            index += 1;
        }

        index = 0;
        ArrayList<Double> knots2 = s.getKnots();
        ArrayList<Double> knotsY = new ArrayList<>(Collections.nCopies(knots2.size(), 0.0));
        for (Double point : knots2) {
            knotsY.set(index, s.getValue(point));
            index += 1;
        }

        File myObj = new File("result.txt");
        try {
            myObj.createNewFile();
        } catch (IOException e) {
            e.printStackTrace();
        }
        try {
            FileWriter myWriter = new FileWriter("result.txt");

            myWriter.write("x_curve:  ");
            for (Double val: xCurve){
                myWriter.write(val + " ");
            }
            myWriter.write(System.lineSeparator() + System.lineSeparator() + System.lineSeparator());


            myWriter.write("y_curve:  ");
            for (Double val: yCurve){
                myWriter.write(val + " ");
            }
            myWriter.write(System.lineSeparator() + System.lineSeparator() + System.lineSeparator());

            
            myWriter.write("y_curve_opt:  ");
            for (Double val: yCurveOpt){
                myWriter.write(val + " ");
            }
            myWriter.write(System.lineSeparator() + System.lineSeparator() + System.lineSeparator());


            myWriter.write("knots:  ");
            for (Double val: knots){
                myWriter.write(val + " ");
            }
            myWriter.write(System.lineSeparator() + System.lineSeparator() + System.lineSeparator());


            myWriter.write("knots_y_uni:  ");
            for (Double val: knotsYUni){
                myWriter.write(val + " ");
            }
            myWriter.write(System.lineSeparator() + System.lineSeparator() + System.lineSeparator());

            
            myWriter.write("knots2:  ");
            for (Double val: knots2){
                myWriter.write(val + " ");
            }
            myWriter.write(System.lineSeparator() + System.lineSeparator() + System.lineSeparator());


            myWriter.write("knots_y:  ");
            for (Double val: knotsY){
                myWriter.write(val + " ");
            }
            myWriter.write(System.lineSeparator() + System.lineSeparator() + System.lineSeparator());
            myWriter.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}