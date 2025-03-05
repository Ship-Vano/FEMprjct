//
// Created by ivan on 05.03.2025.
//

#ifndef FEMPRJCT_FEMSOLVER2D_H
#define FEMPRJCT_FEMSOLVER2D_H

#include "NetGeometry.h"

class FEMSolver2D {
public:
    World geometryWorld;

    std::vector<double> elemUs;

    int task_type = 1;
    double finalTime = 0.0;
    int iterationsPerFrame = 10.0;

    //г.у.
    int bcTypeBot = 1;
    int bcTypeTop = 1;
    int bcTypeLeft = 1;
    int bcTypeRight = 1;

    FEMSolver2D(const World& world);
    double runSolver();
    double lambda(const double& x, const double& y){
        return 1.0;
    }
    double q(const double& x, const double& y){
        return 0.0;
    }
    double exactSol(const double& x, const double& y){
        return x*x - y*y + 2.0;
    }

    void writeVTU(const std::string& filename);
};


#endif //FEMPRJCT_FEMSOLVER2D_H
