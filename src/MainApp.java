import java.io.FileWriter;
import java.io.IOException;

public class MainApp{
    public MainApp(){
    }

    public static void main (String[] args){

	// check number arguments
	if(args.length==0) {
		System.out.println("Specify a filename where the mesh will be written in the OFF file format.");
		return;
	}

	// create a sphere
	// Sphere sphere = new Sphere (1,0,0,0);
	
	// create the test object
	Model model = new Model();

	// create a polygonizer based on dual contouring 
	DualContouring dc = new DualContouring();
	dc.setFunction(model);
	dc.setGridSize(128,128,128);
	dc.setBounds(-2,-2,-2,12,12,12);
	dc.generateMesh();
	// get the mesh of triangles
	TriMesh triangles = dc.t_mesh;

	// prints result out:
	try{
	    FileWriter fw = new FileWriter(args[0]);
	    fw.write("OFF\n");
	    fw.write(triangles.vertex_N + " " + triangles.face_N + " 0\n");
	    for(int i=0; i<triangles.vertex_N; i++){
		fw.write(triangles.vertex[i][0] + " ");
		fw.write(triangles.vertex[i][1] + " ");
		fw.write(triangles.vertex[i][2] + "\n");
	    }
	    for(int i=0; i<triangles.face_N; i++){
		fw.write("3 ");
		fw.write(triangles.face[i][0] + " ");
                fw.write(triangles.face[i][1] + " ");
                fw.write(triangles.face[i][2] + "\n");
	    }
	    fw.close();
	}
	catch (IOException e){
	    System.out.println(e);
	}
    }
}
