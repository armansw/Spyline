package spyline;

import java.util.ArrayList;

public class Points {
    Points(ArrayList<Double> x, ArrayList<Double> y, ArrayList<Double> w){
        int minLen = Math.min(Math.min(x.size(), y.size()), w.size());
        this.x = new ArrayList<>(x.subList(0, minLen));
        this.y = new ArrayList<>(y.subList(0, minLen));
        this.w = new ArrayList<>(w.subList(0, minLen));
    }

    public int len(){
        return this.x.size();
    }

    public Point getItem(int index){
        return new Point(this.x.get(index), this.y.get(index), this.w.get(index));
    }

    private ArrayList<Double> x;
    private ArrayList<Double> y;
    private ArrayList<Double> w;
}