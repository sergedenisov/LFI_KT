#include <iostream>
#include <cmath>
#include <vector>

#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkTetra.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkSmartPointer.h>

#include <gmsh.h>

using namespace std;

// Класс расчётной точки
class CalcNode
{
friend class CalcMesh;

protected:
    double x, y, z;
    double smth;  // Скалярное поле (спиральная волна)
    double vx, vy, vz;  // Скорости для вращения

public:
    CalcNode() : x(0.0), y(0.0), z(0.0), smth(0.0), vx(0.0), vy(0.0), vz(0.0) {}

    CalcNode(double x, double y, double z, double smth, double vx, double vy, double vz) 
        : x(x), y(y), z(z), smth(smth), vx(vx), vy(vy), vz(vz) {}

    // Перемещение точки согласно скорости
    void move(double tau) {
        x += vx * tau;
        y += vy * tau;
        z += vz * tau;
    }

    // Обновление скалярного поля (спиральная волна)
    void updateSmth(double t, double x0, double y0) {
        double dx = x - x0;
        double dy = y - y0;
        double r = sqrt(dx*dx + dy*dy);
        double angle = atan2(dy, dx);
        smth = sin(3*angle + 5*r - 10*t); 
    }
};

// Класс элемента сетки
class Element
{
friend class CalcMesh;

protected:
    unsigned long nodesIds[4];
};

// Класс расчётной сетки
class CalcMesh
{
protected:
    vector<CalcNode> nodes;
    vector<Element> elements;
    double current_time = 0.0;  // Текущее время для анимации
    const double omega = 2.0;   // Угловая скорость вращения (рад/с)

public:
    CalcMesh(const vector<double>& nodesCoords, const vector<size_t>& tetrsPoints) {
        nodes.resize(nodesCoords.size() / 3);
        
        // Инициализация точек с вращением вокруг оси Z
        for(size_t i = 0; i < nodes.size(); i++) {
            double x = nodesCoords[i*3];
            double y = nodesCoords[i*3 + 1];
            double z = nodesCoords[i*3 + 2];
            
            // Скорости для вращения
            double vx = -omega * y;
            double vy = omega * x;
            double vz = 0.0;

            nodes[i] = CalcNode(x, y, z, 0.0, vx, vy, vz);
        }

        // Инициализация элементов
        elements.resize(tetrsPoints.size() / 4);
        for(size_t i = 0; i < elements.size(); i++) {
            elements[i].nodesIds[0] = tetrsPoints[i*4] - 1;
            elements[i].nodesIds[1] = tetrsPoints[i*4 + 1] - 1;
            elements[i].nodesIds[2] = tetrsPoints[i*4 + 2] - 1;
            elements[i].nodesIds[3] = tetrsPoints[i*4 + 3] - 1;
        }
    }

    // Шаг по времени: движение точек + обновление поля
    void doTimeStep(double tau) {
        current_time += tau;
        double center_x = 0.0, center_y = 0.0;  // Центр спиральной волны

        for(size_t i = 0; i < nodes.size(); i++) {
            nodes[i].move(tau);               // Двигаем точку
            nodes[i].updateSmth(current_time, center_x, center_y);  // Обновляем поле
        }
    }

    // Запись снапшота в VTK
    void snapshot(size_t snap_number) {
        vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
        vtkSmartPointer<vtkPoints> dumpPoints = vtkSmartPointer<vtkPoints>::New();

        // Скалярное поле (спиральная волна)
        auto smth = vtkSmartPointer<vtkDoubleArray>::New();
        smth->SetName("smth");

        // Векторное поле (скорости)
        auto vel = vtkSmartPointer<vtkDoubleArray>::New();
        vel->SetName("velocity");
        vel->SetNumberOfComponents(3);

        for(size_t i = 0; i < nodes.size(); i++) {
            dumpPoints->InsertNextPoint(nodes[i].x, nodes[i].y, nodes[i].z);
            vel->InsertNextTuple3(nodes[i].vx, nodes[i].vy, nodes[i].vz);
            smth->InsertNextValue(nodes[i].smth);
        }

        unstructuredGrid->SetPoints(dumpPoints);
        unstructuredGrid->GetPointData()->AddArray(vel);
        unstructuredGrid->GetPointData()->AddArray(smth);

        // Запись тетраэдров
        for(size_t i = 0; i < elements.size(); i++) {
            auto tetra = vtkSmartPointer<vtkTetra>::New();
            tetra->GetPointIds()->SetId(0, elements[i].nodesIds[0]);
            tetra->GetPointIds()->SetId(1, elements[i].nodesIds[1]);
            tetra->GetPointIds()->SetId(2, elements[i].nodesIds[2]);
            tetra->GetPointIds()->SetId(3, elements[i].nodesIds[3]);
            unstructuredGrid->InsertNextCell(tetra->GetCellType(), tetra->GetPointIds());
        }

        string fileName = "cone-step-" + to_string(snap_number) + ".vtu";
        auto writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
        writer->SetFileName(fileName.c_str());
        writer->SetInputData(unstructuredGrid);
        writer->Write();
    }
};

int main()
{
    const double tau = 0.05;  // Шаг времени
    const unsigned int GMSH_TETR_CODE = 4;

    // Инициализация Gmsh и загрузка STL
    gmsh::initialize();
    gmsh::model::add("cone");

    try {
        gmsh::merge("../cone.stl");
    } catch(...) {
        cerr << "Could not load STL mesh!" << endl;
        gmsh::finalize();
        return -1;
    }

    // Построение сетки
    gmsh::model::mesh::classifySurfaces(40 * M_PI / 180., true, false, 180 * M_PI / 180.);
    gmsh::model::mesh::createGeometry();
    
    vector<pair<int, int>> surfaces;
    gmsh::model::getEntities(surfaces, 2);
    vector<int> surfaceLoops;
    for(auto surf : surfaces) surfaceLoops.push_back(surf.second);
    int volume = gmsh::model::geo::addVolume({gmsh::model::geo::addSurfaceLoop(surfaceLoops)});
    gmsh::model::geo::synchronize();

    gmsh::model::mesh::field::add("MathEval", 1);
    gmsh::model::mesh::field::setString(1, "F", "4");
    gmsh::model::mesh::field::setAsBackgroundMesh(1);

    gmsh::model::mesh::generate(3);

    // Извлечение данных сетки
    vector<double> nodesCoord;
    vector<size_t> nodeTags;
    vector<double> parametricCoord;
    gmsh::model::mesh::getNodes(nodeTags, nodesCoord, parametricCoord);

    vector<int> elementTypes;
    vector<vector<size_t>> elementTags, elementNodeTags;
    gmsh::model::mesh::getElements(elementTypes, elementTags, elementNodeTags);

    vector<size_t>* tetrsNodesTags = nullptr;
    for(size_t i = 0; i < elementTypes.size(); i++) {
        if(elementTypes[i] == GMSH_TETR_CODE) {
            tetrsNodesTags = &elementNodeTags[i];
            break;
        }
    }

    if(!tetrsNodesTags) {
        cerr << "No tetrahedrons found!" << endl;
        gmsh::finalize();
        return -2;
    }

    // Создание расчётной сетки
    CalcMesh mesh(nodesCoord, *tetrsNodesTags);
    gmsh::finalize();

    // Генерация 50 кадров анимации
    for(int step = 0; step < 50; step++) {
        mesh.doTimeStep(tau);
        mesh.snapshot(step);
    }

    return 0;
}
