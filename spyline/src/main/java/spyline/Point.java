package spyline;


public class Point {
    Point(double x, double y, double w){
        this.x = x;
        this.y = y;
        this.w = w;
    }

    public double getW(){
        return this.w;
    }

    public double getX(){
        return this.x;
    }

    public double getY(){
        return this.y;
    }

    public void setW(double w){
        this.w = w;
    }

    public void setX(double x){
        this.x = x;
    }

    public void setY(double y){
        this.y = y;
    }

    private double x;
    private double y;
    private double w;
}