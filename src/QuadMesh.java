
public class QuadMesh  
{
    public double vertex[][]; //[][3]
    public int face[][]; //[][4]
    public int vertex_N;
    public int face_N;
    public int degree[];
    public int max_face_N;
    public double normal_v[][]; //[][3]
    public double normal_f[][]; //[][3]
  
    public QuadMesh(){
	vertex_N = 0;
	face_N = 0;
	vertex = null;
	face = new int[1][4];
	normal_f = new double[1][3];
	degree = null;
	max_face_N = 1;
    }
    
    public void setVertexN(int N){
	vertex_N = N;
	vertex = new double[N][3];
	normal_v = new double[N][3];
	degree = new int[N];
    }
  
    public void setFaceN(int N){
	face_N = N;
	face = new int[N][4];
	normal_f = new double[N][3];
	max_face_N = N;
    }
  
    public void setVertex(int i, double x, double y, double z){
	vertex[i][0] = x;
	vertex[i][1] = y;
	vertex[i][2] = z;
    }
  
    public void setFace(int i, int i0, int i1, int i2, int i3){
	checkFaceN();
	face[i][0] = i0;
	face[i][1] = i1;
	face[i][2] = i2;
	face[i][3] = i3;
    }
  
    public void addFace(int i0, int i1, int i2, int i3, double n[]){
	face_N++;
	checkFaceN();
	int i = face_N - 1;
	face[i][0] = i0;
	face[i][1] = i1;
	face[i][2] = i2;
	face[i][3] = i3;
	normal_f[i][0] = n[0];
	normal_f[i][1] = n[1];
	normal_f[i][2] = n[2];
    }
  
    public void checkFaceN(){
	if(face_N > max_face_N){
	    max_face_N *= 2;
	    int new_face[][] = new int[max_face_N][4];
	    double new_normal[][] = new double[max_face_N][3];
	
	    for(int i=0; i<face_N-1; i++){
		new_face[i][0] = face[i][0];
		new_face[i][1] = face[i][1];
		new_face[i][2] = face[i][2];
		new_face[i][3] = face[i][3];
	    
		new_normal[i][0] = normal_f[i][0];
		new_normal[i][1] = normal_f[i][1];
		new_normal[i][2] = normal_f[i][2];
	    }
	    face = new_face;
	    normal_f = new_normal;
	}
    }
}
