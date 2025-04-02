//
// Created by ivan on 05.03.2025.
//

#ifndef FEMPRJCT_FEMSOLVER2D_H
#define FEMPRJCT_FEMSOLVER2D_H

#include "NetGeometry.h"
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
#include <fstream>
#include <memory>

class FEMSolver2D {
private:
    World& geometryWorld;
    Eigen::SparseMatrix<double> globalStiffnessMatrix;
    Eigen::VectorXd globalLoadVector;
    Eigen::VectorXd solution;

    // Структура для хранения параметров задачи
    struct ProblemParameters {
        std::function<double(double, double)> lambda;
        std::function<double(double, double)> sourceTerm;
        std::function<double(double, double)> dirichletBC;
        std::function<double(double, double)> neumannBC;
    } params;

public:
    explicit FEMSolver2D(World& world) : geometryWorld(world) {}

    // Установка параметров задачи
    void setParameters(
            std::function<double(double, double)> lambda,
            std::function<double(double, double)> source,
            std::function<double(double, double)> dirichlet,
            std::function<double(double, double)> neumann = [](double, double){ return 0.0; })
    {
        params = {lambda, source, dirichlet, neumann};
    }

    // Основная функция сборки и решения системы
    void solve() {
        assembleSystem();
        applyBoundaryConditions();
        solveLinearSystem();
    }

    // Экспорт решения в VTU формат
    void exportToVTU(const std::string& filename) const {
        const auto& nodes = geometryWorld.getNodePool().nodes;
        const auto& elements = geometryWorld.getElementPool().elements;

        std::ofstream file(filename);
        file << "<?xml version=\"1.0\"?>\n";
        file << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
        file << "  <UnstructuredGrid>\n";
        file << "    <Piece NumberOfPoints=\"" << nodes.size()
             << "\" NumberOfCells=\"" << elements.size() << "\">\n";

        writePointData(file, nodes);
        writeCellData(file);
        writePoints(file, nodes);
        writeCells(file, elements);

        file << "    </Piece>\n";
        file << "  </UnstructuredGrid>\n";
        file << "</VTKFile>\n";
    }

private:
    void assembleSystem() {
        const auto& elements = geometryWorld.getElementPool().elements;
        const auto& nodes = geometryWorld.getNodePool().nodes;
        const size_t numNodes = nodes.size();

        globalStiffnessMatrix.resize(numNodes, numNodes);
        globalLoadVector.resize(numNodes);
        globalLoadVector.setZero();

        std::vector<Eigen::Triplet<double>> triplets;

        for (const Element& elem : elements) {
            if (elem.dim != 3) throw std::runtime_error("Non-triangular elements not supported");

            Eigen::Matrix3d localStiffness = computeLocalStiffness(elem, nodes);
            Eigen::Vector3d localLoad = computeLocalLoad(elem, nodes);

            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    triplets.emplace_back(elem.nodeIndexes[i],
                                          elem.nodeIndexes[j],
                                          localStiffness(i, j));
                }
                globalLoadVector(elem.nodeIndexes[i]) += localLoad(i);
            }
        }

        globalStiffnessMatrix.setFromTriplets(triplets.begin(), triplets.end());
    }

    Eigen::Matrix3d computeLocalStiffness(const Element& elem,
                                          const std::vector<Node>& nodes) const {
        const Node& n1 = nodes[elem.nodeIndexes[0]];
        const Node& n2 = nodes[elem.nodeIndexes[1]];
        const Node& n3 = nodes[elem.nodeIndexes[2]];

        Eigen::Matrix3d B;
        B << n2.y - n3.y, n3.y - n1.y, n1.y - n2.y,
                n3.x - n2.x, n1.x - n3.x, n2.x - n1.x;

        double area = 0.5 * ((n2.x - n1.x)*(n3.y - n1.y) - (n3.x - n1.x)*(n2.y - n1.y));
        B /= (2.0 * area);

        Eigen::Matrix3d Ke = B.transpose() * B * area;

        // Усреднение коэффициента λ по узлам
        double lambda = (params.lambda(n1.x, n1.y) +
                         params.lambda(n2.x, n2.y) +
                         params.lambda(n3.x, n3.y)) / 3.0;

        return lambda * Ke;
    }

    Eigen::Vector3d computeLocalLoad(const Element& elem,
                                     const std::vector<Node>& nodes) const {
        const Node& n1 = nodes[elem.nodeIndexes[0]];
        const Node& n2 = nodes[elem.nodeIndexes[1]];
        const Node& n3 = nodes[elem.nodeIndexes[2]];

        double area = 0.5 * ((n2.x - n1.x)*(n3.y - n1.y) - (n3.x - n1.x)*(n2.y - n1.y));
        Eigen::Vector3d Fe;

        // Среднее значение источника по элементу
        double q_avg = (params.sourceTerm(n1.x, n1.y) +
                        params.sourceTerm(n2.x, n2.y) +
                        params.sourceTerm(n3.x, n3.y)) / 3.0;

        Fe << q_avg * area / 3.0,
                q_avg * area / 3.0,
                q_avg * area / 3.0;

        return Fe;
    }

    void applyBoundaryConditions() {
        const double penalty = 1e30;
        auto applyBC = [&](const std::vector<int>& nodes) {
            for (int nodeIdx : nodes) {
                for (Eigen::SparseMatrix<double>::InnerIterator it(globalStiffnessMatrix, nodeIdx); it; ++it) {
                    if (it.col() == nodeIdx) it.valueRef() = penalty;
                    else it.valueRef() = 0.0;
                }
                globalLoadVector[nodeIdx] = penalty * params.dirichletBC(
                        geometryWorld.getNodePool().nodes[nodeIdx].x,
                        geometryWorld.getNodePool().nodes[nodeIdx].y
                );
            }
        };

        applyBC(geometryWorld.getBoundaryLeftNodes());
        applyBC(geometryWorld.getBoundaryRightNodes());
        applyBC(geometryWorld.getBoundaryTopNodes());
        applyBC(geometryWorld.getBoundaryBottomNodes());
    }

    void solveLinearSystem() {
        Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver;
        solver.compute(globalStiffnessMatrix);
        solution = solver.solve(globalLoadVector);
    }

    void writePointData(std::ofstream& file, const std::vector<Node>& nodes) const {
        file << "      <PointData Scalars=\"Solution\">\n";
        file << "        <DataArray type=\"Float64\" Name=\"Solution\" format=\"ascii\">\n";
        for (int i = 0; i < solution.size(); ++i) {
            file << solution[i] << " ";
        }
        file << "\n        </DataArray>\n      </PointData>\n";
    }

    void writeCellData(std::ofstream& file) const {
        file << "      <CellData>\n";
        file << "        <DataArray type=\"Int32\" Name=\"Material\" format=\"ascii\">\n";
        for (size_t i = 0; i < geometryWorld.getElementPool().elements.size(); ++i) {
            file << "0 ";
        }
        file << "\n        </DataArray>\n      </CellData>\n";
    }

    void writePoints(std::ofstream& file, const std::vector<Node>& nodes) const {
        file << "      <Points>\n";
        file << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        for (const Node& n : nodes) {
            file << n.x << " " << n.y << " 0.0 ";
        }
        file << "\n        </DataArray>\n      </Points>\n";
    }

    void writeCells(std::ofstream& file, const std::vector<Element>& elements) const {
        file << "      <Cells>\n";
        file << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
        for (const Element& e : elements) {
            file << e.nodeIndexes[0] << " " << e.nodeIndexes[1] << " " << e.nodeIndexes[2] << " ";
        }
        file << "\n        </DataArray>\n";

        file << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
        int offset = 3;
        for (size_t i = 0; i < elements.size(); ++i) {
            file << offset << " ";
            offset += 3;
        }
        file << "\n        </DataArray>\n";

        file << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
        for (size_t i = 0; i < elements.size(); ++i) {
            file << "5 "; // VTK_TRIANGLE
        }
        file << "\n        </DataArray>\n";
        file << "      </Cells>\n";
    }
};


#endif //FEMPRJCT_FEMSOLVER2D_H
