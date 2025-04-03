#include <dolfin.h>
#include "Laplace.h"  // Generated from Laplace.ufl

using namespace dolfin;

// Define electrode positions
class Electrode1 : public SubDomain {
    bool inside(const Array<double>& x, bool on_boundary) const override {
        return (on_boundary && (x[0] >= 0.0 && x[0] <= 0.2) && (x[1] < DOLFIN_EPS));
    }
};

class Electrode2 : public SubDomain {
    bool inside(const Array<double>& x, bool on_boundary) const override {
        return (on_boundary && (x[0] >= 0.2 && x[0] <= 0.4) && (x[1] < DOLFIN_EPS));
    }
};

class Electrode3 : public SubDomain {
    bool inside(const Array<double>& x, bool on_boundary) const override {
        return (on_boundary && (x[0] >= 0.4 && x[0] <= 0.6) && (x[1] < DOLFIN_EPS));
    }
};

class Electrode4 : public SubDomain {
    bool inside(const Array<double>& x, bool on_boundary) const override {
        return (on_boundary && (x[0] >= 0.6 && x[0] <= 0.8) && (x[1] < DOLFIN_EPS));
    }
};

class Electrode5 : public SubDomain {
    bool inside(const Array<double>& x, bool on_boundary) const override {
        return (on_boundary && (x[0] >= 0.8 && x[0] <= 1.0) && (x[1] < DOLFIN_EPS));
    }
};

// Expression to compute electric field magnitude |E| = |∇φ|
class ElectricFieldMagnitude : public Expression {
public:
    ElectricFieldMagnitude(std::shared_ptr<Function> u) 
        : Expression(), u(u) {}
    
    void eval(Array<double>& values, const Array<double>& x) const override {
        // Finite difference approximation with boundary checks
        const double eps = 1e-6;
        Array<double> x_plus(2), x_minus(2), phi(1), phi_plus(1), phi_minus(1);
        
        // Get current point value
        u->eval(phi, x);
        
        // x-derivative
        double Ex = 0.0;
        if (x[0] <= eps) {  // Left boundary
            x_plus[0] = x[0] + eps; x_plus[1] = x[1];
            u->eval(phi_plus, x_plus);
            Ex = (phi_plus[0] - phi[0])/eps;
        }
        else if (x[0] >= 1.0 - eps) {  // Right boundary
            x_minus[0] = x[0] - eps; x_minus[1] = x[1];
            u->eval(phi_minus, x_minus);
            Ex = (phi[0] - phi_minus[0])/eps;
        }
        else {  // Interior point
            x_plus[0] = x[0] + eps; x_plus[1] = x[1];
            x_minus[0] = x[0] - eps; x_minus[1] = x[1];
            u->eval(phi_plus, x_plus);
            u->eval(phi_minus, x_minus);
            Ex = (phi_plus[0] - phi_minus[0])/(2*eps);
        }
        
        // y-derivative
        double Ey = 0.0;
        if (x[1] <= eps) {  // Bottom boundary
            x_plus[0] = x[0]; x_plus[1] = x[1] + eps;
            u->eval(phi_plus, x_plus);
            Ey = (phi_plus[0] - phi[0])/eps;
        }
        else if (x[1] >= 1.0 - eps) {  // Top boundary
            x_minus[0] = x[0]; x_minus[1] = x[1] - eps;
            u->eval(phi_minus, x_minus);
            Ey = (phi[0] - phi_minus[0])/eps;
        }
        else {  // Interior point
            x_plus[0] = x[0]; x_plus[1] = x[1] + eps;
            x_minus[0] = x[0]; x_minus[1] = x[1] - eps;
            u->eval(phi_plus, x_plus);
            u->eval(phi_minus, x_minus);
            Ey = (phi_plus[0] - phi_minus[0])/(2*eps);
        }
        if(sqrt(Ex*Ex + Ey*Ey)>7){
        	values[0] = 7;
        }else{
        	values[0] = sqrt(Ex*Ex + Ey*Ey);  // |E| = |∇φ|
        }
    }
    
private:
    std::shared_ptr<Function> u;
};

int main() {
    // Create mesh (2D unit square)
    auto mesh = std::make_shared<UnitSquareMesh>(64, 64);

    // Define function space
    auto V = std::make_shared<Laplace::FunctionSpace>(mesh);

    // Set up electrodes with different voltages
    auto voltage_RC = std::make_shared<Constant>(100.0);    
    auto voltage_DC = std::make_shared<Constant>(0.0);

    // Create electrode subdomains
    auto electrode_1 = std::make_shared<Electrode1>();
    auto electrode_2 = std::make_shared<Electrode2>();
    auto electrode_3 = std::make_shared<Electrode3>();
    auto electrode_4 = std::make_shared<Electrode4>();
    auto electrode_5 = std::make_shared<Electrode5>();

    // Create boundary conditions
    std::vector<const DirichletBC*> bcs;
    DirichletBC bc_left(V, voltage_DC, electrode_1);
    DirichletBC bc_right(V, voltage_RC, electrode_2);
    DirichletBC bc_top(V, voltage_DC, electrode_3);
    DirichletBC bc_bl(V, voltage_RC, electrode_4);
    DirichletBC bc_br(V, voltage_DC, electrode_5);
    
    bcs.push_back(&bc_left);
    bcs.push_back(&bc_right);
    bcs.push_back(&bc_top);
    bcs.push_back(&bc_bl);
    bcs.push_back(&bc_br);

    // Define variational forms (Laplace equation ∇²φ = 0)
    Laplace::BilinearForm a(V, V);
    Laplace::LinearForm L(V);
    L.f = std::make_shared<Constant>(0.0);  // Zero source term

    // Solve for electric potential
    Function u(V);
    solve(a == L, u, bcs);

    // Save potential solution
    File potential_file("potential.pvd");
    potential_file << u;

    // Create function space for scalar field (using same space as potential)
    auto V_scalar = V;

    // Make a non-const copy of the solution
    auto u_nonconst = std::make_shared<Function>(V);
    *u_nonconst = u;
    u_nonconst->set_allow_extrapolation(true);

    // Compute electric field magnitude
    auto E_mag_expr = std::make_shared<ElectricFieldMagnitude>(u_nonconst);
    Function E_magnitude(V_scalar);
    E_magnitude.interpolate(*E_mag_expr);

    // Save electric field magnitude
    File emag_file("electric_field_magnitude.pvd");
    emag_file << E_magnitude;

    return 0;
}
