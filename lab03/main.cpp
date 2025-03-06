#include <dolfin.h>
#include "Poisson.h"

using namespace dolfin;

// Source term (right-hand side)
class Source : public Expression
{
public:
    void eval(Array<double>& values, const Array<double>& x) const
    {
        double dx = x[0] - 0.5;
        double dy = x[1] - 0.5;
        values[0] = 100 * exp(-(dx*dx + dy*dy) / 0.02); // Increased amplitude
    }
};

// Boundary condition
class DirichletBoundary : public SubDomain
{
    bool inside(const Array<double>& x, bool on_boundary) const
    {
        return on_boundary;
    }
};

int main()
{
    // Create L-shaped mesh manually
    auto mesh = std::make_shared<Mesh>();
    MeshEditor editor;
    
    // Initialize the mesh editor
    editor.open(*mesh, CellType::Type::triangle, 2, 2); // 2D triangular mesh
    editor.init_vertices(6); // Specify the number of vertices
    editor.init_cells(4);    // Specify the number of cells

    // Add vertices
    editor.add_vertex(0, 0.0, 0.0);
    editor.add_vertex(1, 1.0, 0.0);
    editor.add_vertex(2, 1.0, 0.5);
    editor.add_vertex(3, 0.5, 0.5);
    editor.add_vertex(4, 0.5, 1.0);
    editor.add_vertex(5, 0.0, 1.0);

    // Add cells (triangles)
    editor.add_cell(0, 0, 1, 3); // Triangle 1
    editor.add_cell(1, 1, 2, 3); // Triangle 2
    editor.add_cell(2, 3, 4, 5); // Triangle 3
    editor.add_cell(3, 0, 3, 5); // Triangle 4

    // Close the mesh editor
    editor.close();

    // Refine the mesh
    auto refined_mesh = std::make_shared<Mesh>(refine(*mesh)); // Dereference mesh

    // Save refined mesh to file for inspection
    File mesh_file("L_shaped_mesh_refined.pvd");
    mesh_file << *refined_mesh;

    // Print mesh information
    std::cout << "Number of vertices: " << refined_mesh->num_vertices() << std::endl;
    std::cout << "Number of cells: " << refined_mesh->num_cells() << std::endl;

    // Create function space
    auto V = std::make_shared<Poisson::FunctionSpace>(refined_mesh);

    // Define boundary condition
    auto u0 = std::make_shared<Constant>(0.0);
    auto boundary = std::make_shared<DirichletBoundary>();
    DirichletBC bc(V, u0, boundary);

    // Print boundary information
    std::cout << "Applying boundary condition to all boundaries." << std::endl;

    // Define variational forms
    Poisson::BilinearForm a(V, V);
    Poisson::LinearForm L(V);
    auto f = std::make_shared<Source>();
    L.f = f;

    // Compute solution
    Function u(V);
    solve(a == L, u, bc);

    // Print solution information
    std::cout << "Solution vector size: " << u.vector()->size() << std::endl;
    std::cout << "Solution min value: " << u.vector()->min() << std::endl;
    std::cout << "Solution max value: " << u.vector()->max() << std::endl;

    // Save solution to file
    File file("poisson_L_shape_refined.pvd");
    file << u;

    return 0;
}
