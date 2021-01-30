import java.util.*;

import java.io.FileWriter;
import java.io.IOException;

public class DualContouring{
    public int gX, gY, gZ;
    public double minX, maxX, minY, maxY, minZ, maxZ;
    public double value[][][];
    public QuadMesh mesh;
    public TriMesh t_mesh;
    public Function function;
  
    public DualContouring(){
	value = null;
	mesh = null;
	gX = gY = gZ = 0;
    }
    
    public void generateMesh(){
	sample();
	generateQuadMesh();
	createTriangles();
    }
  
    public void setGridSize(int x, int y, int z){
	gX = x + 2;
	gY = y + 2;
	gZ = z + 2;
    }
  
    public void setBounds(double minX, double minY,
			  double minZ, double maxX,
			  double maxY, double maxZ){
	this.minX = minX;
	this.minY = minY;
	this.minZ = minZ;
	this.maxX = maxX;
	this.maxY = maxY;
	this.maxZ = maxZ;
    }
  
    public void setFunction(Function f){
	this.function = f;
    }

    public void sample(){
	double stepX, stepY, stepZ;
	double x, y, z;
     
	stepX = (maxX - minX)/(gX-1);
	stepY = (maxY - minY)/(gY-1);
	stepZ = (maxZ - minZ)/(gZ-1);
     
	value = new double[gZ][gY][gX];
    
	int i,j,k;
	for(i=0; i<gZ; i++){
	    z = minZ + i*stepZ;
	    for(j=0; j<gY; j++){
		y = minY + j*stepY;
		for(k=0; k<gX; k++){
		    x = minX + k*stepX;
		    value[i][j][k] = function.evaluate(x, y, z);
		}
	    }
	}
    }
  
    public void generateQuadMesh(){
	int i, j, k;
	int cube[][][] = new int[gZ-1][gY-1][gX-1];
	int vertexN = 0;
	for(i=0; i<gZ-1; i++){
	    for(j=0; j<gY-1; j++){
		for(k=0; k<gX-1; k++){
		    boolean flag = (value[i][j][k] > 0);
		    if(flag != (value[i][j][k+1] > 0) ||
		       flag != (value[i][j+1][k] > 0) ||
		       flag != (value[i][j+1][k+1] > 0) ||
		       flag != (value[i+1][j][k] > 0) ||
		       flag != (value[i+1][j][k+1] > 0) ||
		       flag != (value[i+1][j+1][k] > 0) ||
		       flag != (value[i+1][j+1][k+1] > 0)){
			cube[i][j][k] = vertexN++;
		    }
		    else{
			cube[i][j][k] = -1;
		    }
		}
	    }
	}

	mesh = new QuadMesh();
	mesh.setVertexN(vertexN);

	double iX = (maxX-minX)/(gX-1);
	double iY = (maxY-minY)/(gY-1);
	double iZ = (maxZ-minZ)/(gZ-1);
	//Cube in x direction
	for(i=1; i<gZ-1; i++){
	    double z = minZ + i*iZ;
	    for(j=1; j<gY-1; j++){
		double y = minY + j*iY;
		for(k=1; k<gX-2; k++){
		    double v1 = value[i][j][k];
		    double v2 = value[i][j][k+1];
		    if(v1*v2 < 0){
			double x1 = minX + k*iX;
			double x2 = minX + (k+1)*iX;

			double p[] = new double[3];
			root(p, x1, y, z, v1, x2, y, z, v2);
			double n[] = new double[3];
			grad(n, p[0], p[1], p[2]);
			double len = Math.sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);            
			if((double)len == 0)
			    continue;
			n[0] = (double)(n[0]/len);
			n[1] = (double)(n[1]/len);
			n[2] = (double)(n[2]/len);

			int i1 = cube[i-1][j-1][k];
			int i2 = cube[i-1][j][k];
			int i3 = cube[i][j-1][k];
			int i4 = cube[i][j][k];
            
			if(i1 < 0 || i2 < 0 || i3 < 0 || i4 < 0){
			    System.out.println("Error: found a strange edge?");
			    continue;
			}

			if(v1 > v2)
			    mesh.addFace(i1, i2, i4, i3, n);
			else
			    mesh.addFace(i2, i1, i3, i4, n);

			mesh.degree[i1]++;
			mesh.degree[i2]++;
			mesh.degree[i3]++;
			mesh.degree[i4]++;
            
			mesh.vertex[i1][0] += p[0];
			mesh.vertex[i1][1] += p[1];
			mesh.vertex[i1][2] += p[2];
            
			mesh.vertex[i2][0] += p[0];
			mesh.vertex[i2][1] += p[1];
			mesh.vertex[i2][2] += p[2];
            
			mesh.vertex[i3][0] += p[0];
			mesh.vertex[i3][1] += p[1];
			mesh.vertex[i3][2] += p[2];
            
			mesh.vertex[i4][0] += p[0];
			mesh.vertex[i4][1] += p[1];
			mesh.vertex[i4][2] += p[2];
		    }
		}
	    }
	}

	//Cube in y direction
	for(i=1; i<gX-1; i++){
	    double x = minX + i*iX;
	    for(j=1; j<gZ-1; j++){
		double z = minZ + j*iZ;
		for(k=1; k<gY-2; k++){
		    double v1 = value[j][k][i];
		    double v2 = value[j][k+1][i];
		    if(v1*v2 < 0){
			double y1 = minY + k*iY;
			double y2 = minY + (k+1)*iY;
            
			double p[] = new double[3];
			root(p, x, y1, z, v1, x, y2, z, v2);
			double n[] = new double[3];
			grad(n, p[0], p[1], p[2]);
			double len = Math.sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
			if((double)len == 0)
			    continue;
			n[0] = (double)(n[0]/len);
			n[1] = (double)(n[1]/len);
			n[2] = (double)(n[2]/len);
            
			int i1 = cube[j-1][k][i-1];
			int i2 = cube[j][k][i-1];
			int i3 = cube[j-1][k][i];
			int i4 = cube[j][k][i];

			if(i1 < 0 || i2 < 0 || i3 < 0 || i4 < 0){
			    System.out.println( "Error: found a strange edge?");
			    continue;
			}

			if(v1 > v2)
			    mesh.addFace(i1, i2, i4, i3, n);
			else
			    mesh.addFace(i2, i1, i3, i4, n);
            
			mesh.degree[i1]++;
			mesh.degree[i2]++;
			mesh.degree[i3]++;
			mesh.degree[i4]++;
            
			mesh.vertex[i1][0] += p[0];
			mesh.vertex[i1][1] += p[1];
			mesh.vertex[i1][2] += p[2];
            
			mesh.vertex[i2][0] += p[0];
			mesh.vertex[i2][1] += p[1];
			mesh.vertex[i2][2] += p[2];
            
			mesh.vertex[i3][0] += p[0];
			mesh.vertex[i3][1] += p[1];
			mesh.vertex[i3][2] += p[2];
            
			mesh.vertex[i4][0] += p[0];
			mesh.vertex[i4][1] += p[1];
			mesh.vertex[i4][2] += p[2];
		    }
		}
	    }
	}
    
	//Cube in z direction
	for(i=1; i<gY-1; i++){
	    double y = minY + i*iY;
	    for(j=1; j<gX-1; j++){
		double x = minX + j*iX;
		for(k=1; k<gZ-2; k++){
		    double v1 = value[k][i][j];
		    double v2 = value[k+1][i][j];
		    if(v1*v2 < 0){
			double z1 = minZ + k*iZ;
			double z2 = minZ + (k+1)*iZ;
            
			double p[] = new double[3];
			root(p, x, y, z1, v1, x, y, z2, v2);
			double n[] = new double[3];
			grad(n, p[0], p[1], p[2]);
			double len = Math.sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
			if((double)len == 0)
			    continue;
			n[0] = (double)(n[0]/len);
			n[1] = (double)(n[1]/len);
			n[2] = (double)(n[2]/len);
            

			int i1 = cube[k][i-1][j-1];
			int i2 = cube[k][i-1][j];
			int i3 = cube[k][i][j-1];
			int i4 = cube[k][i][j];
            

			if(i1 < 0 || i2 < 0 || i3 < 0 || i4 < 0){
			    System.out.println("found a strange vertex");
			    continue;
			}
         
			if(v1 > v2)
			    mesh.addFace(i1, i2, i4, i3, n);
			else
			    mesh.addFace(i2, i1, i3, i4, n);
            
			mesh.degree[i1]++;
			mesh.degree[i2]++;
			mesh.degree[i3]++;
			mesh.degree[i4]++;
            
			mesh.vertex[i1][0] += p[0];
			mesh.vertex[i1][1] += p[1];
			mesh.vertex[i1][2] += p[2];
            
			mesh.vertex[i2][0] += p[0];
			mesh.vertex[i2][1] += p[1];
			mesh.vertex[i2][2] += p[2];
            
			mesh.vertex[i3][0] += p[0];
			mesh.vertex[i3][1] += p[1];
			mesh.vertex[i3][2] += p[2];
            
			mesh.vertex[i4][0] += p[0];
			mesh.vertex[i4][1] += p[1];
			mesh.vertex[i4][2] += p[2];
		    }
		}
	    }
	}

	for(i=0; i<vertexN; i++){
	    if(mesh.degree[i] == 0)
		continue;
	    double ideg = 1.0f/mesh.degree[i];
	    //centroid of edge points
	    double px = mesh.vertex[i][0] *= ideg;
	    double py = mesh.vertex[i][1] *= ideg;
	    double pz = mesh.vertex[i][2] *= ideg;
	}
    }

    public void createTriangles(){
	int i;
	int fN = mesh.face_N;
	int new_face[][] = new int[2*fN][3];
	for(i=0; i<fN; i++){
	    int i1 = mesh.face[i][0];
	    int i2 = mesh.face[i][1];
	    int i3 = mesh.face[i][2];
	    int i4 = mesh.face[i][3];
      
	    double a1, a2, min1, min2;
	    //cut along (i1,i3)
	    a1 = angle(mesh.vertex[i1], mesh.vertex[i2], mesh.vertex[i3]);
	    a2 = angle(mesh.vertex[i3], mesh.vertex[i4], mesh.vertex[i1]);
	    if(a1 < a2)
		min1 = a1;
	    else
		min1 = a2;
	    //cut along (i2,i3)
	    a1 = angle(mesh.vertex[i1], mesh.vertex[i2], mesh.vertex[i4]);
	    a2 = angle(mesh.vertex[i3], mesh.vertex[i4], mesh.vertex[i2]);
	    if(a1 < a2)
		min2 = a1;
	    else
		min2 = a2;
      
	    if(min1 > min2){
		new_face[2*i][0] = i1;
		new_face[2*i][1] = i2;
		new_face[2*i][2] = i3;
        
		new_face[2*i+1][0] = i3;
		new_face[2*i+1][1] = i4;
		new_face[2*i+1][2] = i1;
	    }
	    else{
		new_face[2*i][0] = i1;
		new_face[2*i][1] = i2;
		new_face[2*i][2] = i4;
        
		new_face[2*i+1][0] = i3;
		new_face[2*i+1][1] = i4;
		new_face[2*i+1][2] = i2;
	    }
	}
    
	//Compute vertex normals
	int vertexN = mesh.vertex_N;
    
	//Construct output mesh
	fN *= 2;

	// creates a mesh of triangles for output: ply2, povray, stl
	t_mesh = new TriMesh();
	t_mesh.setVertexN(mesh.vertex_N);
	t_mesh.setFaceN(fN);
	for(i=0;i<mesh.vertex_N;i++){
	    t_mesh.vertex[i][0]=mesh.vertex[i][0];
	    t_mesh.vertex[i][1]=mesh.vertex[i][1];
	    t_mesh.vertex[i][2]=mesh.vertex[i][2];
	}
	for(i=0;i<fN;i++){
	    t_mesh.face[i][0]=new_face[i][0];
	    t_mesh.face[i][1]=new_face[i][1];
	    t_mesh.face[i][2]=new_face[i][2];
	}
    }
  
    double angle(double p0[], double p1[], double p2[]){
	double v1[] = new double[] {p1[0]-p0[0], p1[1]-p0[1], p1[2]-p0[2]};
	double v2[] = new double[] {p2[0]-p0[0], p2[1]-p0[1], p2[2]-p0[2]};
	double v3[] = new double[] {p2[0]-p1[0], p2[1]-p1[1], p2[2]-p1[2]};
	double dot1 = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
	double dot2 = -(v1[0]*v3[0] + v1[1]*v3[1] + v1[2]*v3[2]);
	double dot3 = v3[0]*v2[0] + v3[1]*v2[1] + v3[2]*v2[2];
	double dot = dot1;
	if(dot < dot2)
	    dot = dot2;
	if(dot < dot3)
	    dot = dot3;
	return (double)Math.acos((double)dot);
    }
  
    double grad(double g[], double x, double y, double z){
	double f = function.evaluate(x, y, z);
	double diagonal = (double)Math.sqrt((maxX - minX)*(maxX - minX) + 
					    (maxY - minY)*(maxY - minY) +
					    (maxZ - minZ)*(maxZ - minZ));
	double h = 0.1f * diagonal;
	g[0] = (function.evaluate(x+h, y, z) - f)/h;
	g[1] = (function.evaluate(x, y+h, z) - f)/h;
	g[2] = (function.evaluate(x, y, z+h) - f)/h;
	return f;
    }
  
    double root(double p[],
		double x1, double y1, double z1, double v1,
		double x2, double y2, double z2, double v2){
	double x=0, y=0, z=0, v=0;
	int i;

	int IMAX = 10;
	double eps;
	double diagonal = (double)Math.sqrt((maxX - minX)*(maxX - minX) + 
					    (maxY - minY)*(maxY - minY) + 
					    (maxZ - minZ)*(maxZ - minZ));
	eps = 0.0001f * diagonal;

	for(i=0; i<IMAX; i++){
	    double w1 = Math.abs(v2);
	    double w2 = Math.abs(v1);
	    double w = w1 + w2;
	    if(w == 0){
		x = x1;
		y = y1;
		z = z1;
		v = 0;
		break;
	    }
	    w1 /= w;
	    w2 /= w;
      
	    x = w1*x1 + w2*x2;
	    y = w1*y1 + w2*y2;
	    z = w1*z1 + w2*z2;
      
	    v = function.evaluate(x, y, z);
      
	    if(Math.abs(v) < eps)
		break;
	    if(v1*v > 0){
		x1 = x;
		y1 = y;
		z1 = z;
		v1 = v;
	    }
	    else{
		x2 = x;
		y2 = y;
		z2 = z;
		v2 = v;
	    }
	}
	//Change into Bisection
	if(i == IMAX){
	    int j=0;
	    for(j=0; j<2*IMAX; j++){
		x = 0.5f*(x1 + x2);
		y = 0.5f*(y1 + y2);
		z = 0.5f*(z1 + z2);
        
		v = function.evaluate(x, y, z);
        
		if(Math.abs(v) < eps)
		    break;
		if(v1*v > 0){
		    x1 = x;
		    y1 = y;
		    z1 = z;
		    v1 = v;
		}
		else{
		    x2 = x;
		    y2 = y;
		    z2 = z;
		    v2 = v;
		}
	    }
	}
	p[0] = x;
	p[1] = y;
	p[2] = z;
    
	return v;
    }
}
