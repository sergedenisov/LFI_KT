#include <gmsh.h>
#include <vector>

int main(int argc, char **argv) {
    gmsh::initialize(argc, argv);
    gmsh::model::add("toroidal_chamber");


    double R = 10.0, r_ext = 3.0, r_int = 2.0;
    double delta = r_ext - r_int;
    double mesh_size = delta/4;


    int outter = gmsh::model::occ::addTorus(0, 0, 0, R, r_ext);
    int inner = gmsh::model::occ::addTorus(0, 0, 0, R, r_int);
    gmsh::model::occ::synchronize();


    std::vector<std::pair<int, int>> object = {{3, outter}};
    std::vector<std::pair<int, int>> tool = {{3, inner}};
    std::vector<std::pair<int, int>> diff_result;
    std::vector<std::vector<std::pair<int, int>>> map;
    
    gmsh::model::occ::cut(
        object,
        tool,
        diff_result,
        map,
        -1,
        true,
        true
    );
    
    gmsh::model::occ::synchronize();


    gmsh::option::setNumber("Mesh.CharacteristicLengthMin", mesh_size);
    gmsh::option::setNumber("Mesh.CharacteristicLengthMax", mesh_size);
    gmsh::option::setNumber("Mesh.Algorithm3D", 10);


    gmsh::model::mesh::generate(3);
    gmsh::write("toroidal_chamber.msh");
    gmsh::finalize();
    return 0;
}
