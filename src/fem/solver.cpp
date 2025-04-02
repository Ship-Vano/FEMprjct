//
// Created by ivan on 05.03.2025.
//
#include "NetGeometry.h"
#include "FEMSolver2D.h"
#include "json/json.h"

int main() {
    std::string configPath = "InputData/solverConfig.json";
    std::ifstream in(configPath, std::ios::in);

    Json::Reader json_reader;
    Json::Value json_root;
    bool read_succeeded = json_reader.parse(in, json_root);
    assert(read_succeeded);

    std::string meshBinFileName = json_root.get("importFileName", "").asString();

    // Создание и загрузка сетки
    World world(meshBinFileName, true);
    std::cout << "Corners: minX = " << world.minX << " , maxX = " << world.maxX << " , minY = " << world.minY << " , maxY = " << world.maxY << std::endl;
    std::cout << "Loaded boundaries: "
              << world.getBoundaryLeftNodes().size() << " left, "
              << world.getBoundaryRightNodes().size() << " right, "
              << world.getBoundaryTopNodes().size() << " top, "
              << world.getBoundaryBottomNodes().size() << " bottom nodes\n";

    // Инициализация решателя
    FEMSolver2D solver(world);

    // Установка параметров задачи (пример)
    solver.setParameters(
            [](double x, double y){ return 1.0; },           // λ
            [](double x, double y){ return 0.0; },           // Источник
            [](double x, double y){ return x*x - y*y + 2; } // Условие Дирихле
    );

    // Решение задачи
    solver.solve();

    // Экспорт результатов
    solver.exportToVTU("solution.vtu");

    return 0;
}


/**
 * orbtvcfgdfg()
 * frifrfrfewfdsf()
 * hello there
 * how are you
 * gfh ghbdtn dfsr
 * fds3 как  еоаывафывццловссчяя
 * куываааавв
 * frfsdzcxff
 * fefm ds orint( fgdf )
 * hello wirkd ()  if toy for (int i =0; i < n; ++ i){
 * for i in range(3):
 *     input()
 *     ff()
 *     cxz()
 *     zxc()
 *     sdfqwe = print_f_v_cdv() + 2
 */