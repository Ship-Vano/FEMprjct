//
// Created by ivan on 05.03.2025.
//

#include "FEMSolver2D.h"

//FEMSolver2D::FEMSolver2D(const World &world): geometryWorld(world){}
//
//void FEMSolver2D::writeVTU(const std::string& filename) {
//    const NodePool& np = geometryWorld.getNodePool();
//    const ElementPool& ep = geometryWorld.getElementPool();
//
//
//    std::ofstream file(filename);
//    if (!file.is_open()) {
//        throw std::runtime_error("Cannot open file: " + filename);
//    }
//
//    file << "<?xml version=\"1.0\"?>\n";
//    file << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
//    file << "  <UnstructuredGrid>\n";
//    file << "    <Piece NumberOfPoints=\"" << np.nodeCount
//         << "\" NumberOfCells=\"" << ep.elCount << "\">\n";
//
//    // Write points (nodes)
//    file << "      <Points>\n";
//    file << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
//    for (const auto& node : np.nodes) {
//            file << "          " << node.x << " " << node.y << " " << node.z << "\n";
//    }
//    file << "        </DataArray>\n";
//    file << "      </Points>\n";
//
//    // Write cells (elements)
//    file << "      <Cells>\n";
//
//    // Connectivity (fix 1-based index issue)
//    file << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
//    for (const auto& element : ep.elements) {
//            for (const auto &nodeIndex : element.nodeIndexes) {
//                file << "          " << (nodeIndex) << " ";  // Convert to 0-based indexing
//            }
//            file << "\n";
//    }
//    file << "        </DataArray>\n";
//
//    // Offsets
//    file << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
//    int offset = 0;
//    for (const auto& element : ep.elements) {
//            offset += element.dim;
//            file << "          " << offset << "\n";
//    }
//    file << "        </DataArray>\n";
//
//    // Types
//    file << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
//    for (const auto& element : ep.elements) {
//            file << "          " << ((element.dim == 3) ? 5 : 9) << "\n";  // 5 for Triangle, 9 for Quadrilateral
//    }
//    file << "        </DataArray>\n";
//    file << "      </Cells>\n";
//
//    // Write solution data (elemUs) - Fix NumberOfComponents to 8
//    file << "      <CellData Scalars=\"elemUs\">\n";
//    file << "        <DataArray type=\"Float64\" NumberOfComponents=\"9\" Name=\"elemUs\" format=\"ascii\">\n";
//    for (const auto& U : elemUs) {
//        for (const auto& value : U) {
//            file << "          " << value << " ";
//        }
//        file << "\n";
//    }
//    file << "        </DataArray>\n";
//    file << "      </CellData>\n";
//
//    file << "    </Piece>\n";
//    file << "  </UnstructuredGrid>\n";
//    file << "</VTKFile>\n";
//
//    file.close();
//}
