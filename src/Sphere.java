public class Sphere implements Function {

    private double x0;
    private double y0;
    private double z0;
    private double R;

    Sphere(double R, double x0, double y0, double z0){
	this.x0 = x0;
	this.y0 = y0;
	this.z0 = z0;
	this.R = R;
    }

    public double evaluate (double x, double y, double z){
	return R - Math.sqrt((x-x0)*(x-x0) + (y-y0)*(y-y0) + (z-z0)*(z-z0));
    }
}
