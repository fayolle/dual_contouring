public class TriMesh{
    public double vertex[][];
    public int face[][];
    public int vertex_N;
    public int face_N;
    public double normal_v[][];
    public double normal_f[][];

    public TriMesh(){
	vertex_N=0;
	face_N=0;
	vertex=null;
	face=null;
	normal_f=null;
	normal_v=null;
    }

    public void setVertexN(int N){
	vertex_N=N;
	vertex=new double[N][3];
	normal_v=new double[N][3];
    }

    public void setFaceN(int N){
	face_N=N;
	face=new int[N][3];
	normal_f=new double[N][3];
    }
}
