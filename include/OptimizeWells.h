#pragma once
#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>

#include <Eigen/Dense>
#include <nlohmann/json.hpp>

#include <TrajectoryTemplate.h>

double OneWellScore(std::vector<TrajectoryTemplate*>& mainWell,
                    double penalty = 1000);

std::vector<TrajectoryTemplate*> well2CHCH(const std::vector<double>& x,
                                           const Eigen::Vector3d& pinit,
                                           const Eigen::Vector3d& pT1,
                                           const Eigen::Vector3d& pT3);

void OptimizeHorizontal(const Eigen::Vector3d& pinit,
                        const Eigen::Vector3d& pT1, const Eigen::Vector3d& pT3);

void OptimizeHorizontals(const std::vector<Eigen::Vector3d>& pinits,
                         const std::vector<Eigen::Vector3d>& targets1,
                         const std::vector<Eigen::Vector3d>& targets3);

void OptimizeTogether(const std::vector<Eigen::Vector3d>& pinits,
                      const std::vector<Eigen::Vector3d>& targets1,
                      const std::vector<Eigen::Vector3d>& targets3);

void testHorizontal(const Eigen::VectorXd& x,
                    std::vector<TrajectoryTemplate*>& tt);

void writeGG(const std::string& filename, const nlohmann::json& jresults);

void aloneOpt();