/*
 * RpcClusterization.cpp
 *
 *  Created on: Jan 14, 2019
 *      Author: kbunkow
 */

#include "L1Trigger/L1TMuonOverlapPhase1/interface/RpcClusterization.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <cmath>
#include <algorithm>

RpcClusterization::~RpcClusterization() {}

std::vector<RpcCluster> RpcClusterization::getClusters(const RPCDetId& roll, std::vector<RPCDigi>& digis) const {
  std::vector<RpcCluster> allClusters;

  std::sort(digis.begin(), digis.end(), [](const RPCDigi& a, const RPCDigi& b) { return a.strip() < b.strip(); });

  typedef std::pair<unsigned int, unsigned int> Cluster;

  for (auto& digi : digis) {
    if (allClusters.empty()) {
      allClusters.emplace_back(digi.strip(), digi.strip());
      allClusters.back().bx = digi.bx();
      allClusters.back().timing = convertTiming(digi.time());
    } else if (digi.strip() - allClusters.back().lastStrip == 1) {
      allClusters.back().lastStrip = digi.strip();
      //TODO update bx and timing in some smart way
    } else if (digi.strip() - allClusters.back().lastStrip > 1) {
      allClusters.emplace_back(digi.strip(), digi.strip());
      allClusters.back().bx = digi.bx();
      allClusters.back().timing = convertTiming(digi.time());
    }
  }

  std::vector<RpcCluster> filteredClusters;

  std::vector<RpcCluster> emptyResult;

  //modification in the FW from Nov 2021
  if (dropAllClustersIfMoreThanMax)
    if (allClusters.size() > maxClusterCnt)
      return filteredClusters;

  //debug printout only
  if (allClusters.size() > maxClusterCnt) {
    LogTrace("l1tOmtfEventPrint") << __FUNCTION__ << ":" << __LINE__ << " allClusters.size() >= maxClusterCnt "
                                  << std::endl;
    for (auto& cluster : allClusters)
      LogTrace("l1tOmtfEventPrint")
          << __FUNCTION__ << ":" << __LINE__ << " roll " << roll << " cluster: firstStrip " << cluster.firstStrip
          << " lastStrip " << cluster.lastStrip << " halfStrip " << cluster.halfStrip() << std::endl;
  }

  //this is very simple filtering of the clusters,
  //till Nov 2021: unfortunately in the firmware it was more complicated and cannot be easily emulated from digi
  //(in principle would required raws, because in the firmware the clusterizaton is based on the 8-bit strip partitions)
  //the FW from from Nov 2021 solved this problem - option dropAllClustersIfMoreThanMax (above and below)
  //beside betterdata-to-emulator agreement it provides better eff for high pt muons
  for (auto& cluster : allClusters) {
    if (cluster.size() <= maxClusterSize)
      filteredClusters.emplace_back(cluster);
    else if (dropAllClustersIfMoreThanMax)
      return emptyResult;
      //modification in the FW from Nov 2021:
      //if a cluster is bigger then 3 strips - all clusters are dropped

    if (filteredClusters.size() >= maxClusterCnt)
      break;
  }

  return filteredClusters;
}

int RpcClusterization::convertTiming(double timing) const {
  return timing;  //TODO implement
}
