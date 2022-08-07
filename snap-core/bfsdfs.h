#include <vector>
#include <set>
#include <iostream>
#include <queue>
#include <limits>
using namespace std;

namespace TSnap {

/////////////////////////////////////////////////
// BFS and DFS
/// Returns a directed Breadth-First-Search tree rooted at StartNId. ##GetBfsTree1
template <class PGraph> PNGraph GetBfsTree(const PGraph& Graph, const int& StartNId, const bool& FollowOut, const bool& FollowIn);
/// Returns the BFS tree size (number of nodes) and depth (number of levels) by following in-links (parameter FollowIn = true) and/or out-links (parameter FollowOut = true) of node StartNId.
template <class PGraph> int GetSubTreeSz(const PGraph& Graph, const int& StartNId, const bool& FollowOut, const bool& FollowIn, int& TreeSzX, int& TreeDepthX);
/// Finds IDs of all nodes that are at distance Hop from node StartNId. ##GetSubTreeSz
template <class PGraph> int GetNodesAtHop(const PGraph& Graph, const int& StartNId, const int& Hop, TIntV& NIdV, const bool& IsDir=false);
/// Returns the number of nodes at each hop distance from the starting node StartNId. ##GetNodesAtHops
template <class PGraph> int GetNodesAtHops(const PGraph& Graph, const int& StartNId, TIntPrV& HopCntV, const bool& IsDir=false);

/////////////////////////////////////////////////
// Shortest paths
/// Returns the length of the shortest path from node SrcNId to node DstNId. ##GetShortPath1
template <class PGraph> int GetShortPath(const PGraph& Graph, const int& SrcNId, const int& DstNId, const bool& IsDir=false);
/// Returns the length of the shortest path from node SrcNId to all other nodes in the network. ##GetShortPath2
template <class PGraph> int GetShortPath(const PGraph& Graph, const int& SrcNId, TIntH& NIdToDistH, const bool& IsDir=false, const int& MaxDist=TInt::Mx);

/////////////////////////////////////////////////
// Diameter

/// Returns the (approximation of the) Diameter (maximum shortest path length) of a graph (by performing BFS from NTestNodes random starting nodes). ##GetBfsFullDiam
template <class PGraph> int GetBfsFullDiam(const PGraph& Graph, const int& NTestNodes, const bool& IsDir=false);
/// Returns the (approximation of the) Effective Diameter (90-th percentile of the distribution of shortest path lengths) of a graph (by performing BFS from NTestNodes random starting nodes). ##GetBfsEffDiam1
template <class PGraph> double GetBfsEffDiam(const PGraph& Graph, const int& NTestNodes, const bool& IsDir=false);
/// Returns the (approximation of the) Effective Diameter and the Diameter of a graph (by performing BFS from NTestNodes random starting nodes). ##GetBfsEffDiam2
template <class PGraph> double GetBfsEffDiam(const PGraph& Graph, const int& NTestNodes, const bool& IsDir, double& EffDiamX, int& FullDiamX);
/// Returns the (approximation of the) Effective Diameter, the Diameter and the Average Shortest Path length in a graph (by performing BFS from NTestNodes random starting nodes). ##GetBfsEffDiam3
template <class PGraph> double GetBfsEffDiam(const PGraph& Graph, const int& NTestNodes, const bool& IsDir, double& EffDiamX, int& FullDiamX, double& AvgSPLX);
/// Returns the (approximation of the) Effective Diameter, the Diameter and the Average Shortest Path length in a graph (by performing BFS from NTestNodes random starting nodes). ##GetBfsEffDiamAll
template <class PGraph> double GetBfsEffDiamAll(const PGraph& Graph, const int& NTestNodes, const bool& IsDir, double& EffDiamX, int& FullDiamX, double& AvgSPLX);
/// Use the whole graph (all edges) to measure the shortest path lengths but only report the path lengths between nodes in the SubGraphNIdV. ##GetBfsEffDiam4
template <class PGraph> double GetBfsEffDiam(const PGraph& Graph, const int& NTestNodes, const TIntV& SubGraphNIdV, const bool& IsDir, double& EffDiamX, int& FullDiamX);

///////////////////////////////////////////////////
    // Minimum spanning tree
    // Returns the tree and the parent ,child hashtable where hashtable[child] = parent/ pi[child] = parent
    template <class PGraph> std::pair<PUNGraph,TIntIntH> GetMinSpanTree(const PGraph& Graph, const int &StartNId);
    template <class nodewType, class edgewType>
    std::pair<PUNGraph,TIntIntH> GetMinSpanTree(const TPt<TNodeEDatNet<nodewType, edgewType>>& Graph, const int &StartNId);
// TODO: Implement in the future
//template <class PGraph> int GetRangeDist(const PGraph& Graph, const int& SrcNId, const int& DstNId, const bool& IsDir=false);
//template <class PGraph> int GetShortPath(const PGraph& Graph, const int& SrcNId, TIntH& NIdToDistH, const bool& IsDir=false, const int& MaxDist=1000);
//template <class PGraph> int GetShortPath(const PGraph& Graph, const int& SrcNId, const TIntSet& TargetSet, const bool& IsDir, TIntV& PathNIdV);
//template <class PGraph> int GetShortPath(TIntH& NIdPrnH, TCcQueue<int>& NIdQ, const PGraph& Graph, const int& SrcNId, const TIntSet& TargetSet, const bool& IsDir, TIntV& PathNIdV);
//template <class PGraph> int GetMxShortDist(const PGraph& Graph, const int& SrcNId, const bool& IsDir=false);
//template <class PGraph> int GetMxShortDist(const PGraph& Graph, const int& SrcNId, const bool& IsDir, int& MxDistNId);
//template <class PGraph> int GetMxShortDist(const PGraph& Graph, const int& SrcNId, const bool& IsDir, int& MxDistNId, TCcQueue<int>& NIdQ, TCcQueue<int>& DistQ, TIntSet& VisitedH);
//template <class PGraph> int GetMxGreedyDist(const PGraph& Graph, const int& SrcNId, const bool& IsDir=false);
//template <class PGraph> int GetMxGreedyDist(const PGraph& Graph, const int& SrcNId, const bool& IsDir, TCcQueue<int>& NIdQ, TCcQueue<int>& DistQ, TIntSet& VisitedH);
//template <class PGraph> PNGraph GetShortPathsSubGraph(const PGraph& Graph, const TIntV& SubGraphNIdV);
//template <class PGraph> PGraph GetWccPathsSubGraph(const PGraph& Graph, const TIntV& NIdV);
//template <class PGraph> void GetSubTreeSz(const PGraph& Graph, const int& StartNId, const bool& FollowOutEdges, int& TreeSz, int& TreeDepth);

} // namespace TSnap
//typedef TNodeEDatNet<TInt,TFlt> Wgraph;
template<class Wgraph>
class TBreathFSW {
    public:
    PIntFltNEDNet Graph; //int node weight, float edge weight
    //PIntNEDNet Graph; // int node weight, int edge weight
    TSnapQueue<int> Queue;
    TInt StartNId;
    TIntFltH NIdDistH;
    TIntQV VQueues; // typedef TVec<TQQueue<TInt> > TIntQV;
    bool concurrent;
    TIntPrFltH smallest_intra_clust_dist;
    std::vector <std::set<int>> cl_boundary;
    TBreathFSW(const Wgraph& net, const bool& concurrent, const bool& InitBigQ=true):Graph(net),Queue(InitBigQ?Graph->GetNodes():1024), NIdDistH(InitBigQ?Graph->GetNodes():1024),concurrent(true) {}
    void DoBfs_concurrent(const TIntV& StartNodes ,vector <TIntV> &Vbfstree_nodes, const int& MxDist=TInt::Mx){
        const int N_startNodes = StartNodes.Len();
        if(concurrent){
            VQueues.Reserve(N_startNodes);
        }
        NIdDistH.Clr(false);
        for (int i =0; i<N_startNodes; i++){
            auto StartNId = StartNodes[i];
            IAssert(Graph->IsNode(StartNId));
            //  const typename PGraph::TObj::TNodeI StartNodeI = Graph->GetNI(StartNode);
            //  IAssertR(StartNodeI.GetOutDeg() > 0, TStr::Fmt("No neighbors from start node %d.", StartNode));
            NIdDistH.AddDat(StartNId, 0.0);
            VQueues[i].Clr(false);  VQueues[i].Push(StartNId);
            Vbfstree_nodes[i].Add(StartNId);
            cl_boundary.push_back(std::set<int>());
        }
        while(true){
            int empty_counter = 0;
            for (int i = 0; i < N_startNodes ; i++){
                //std::cout << "iter " << i << ": "<<std::endl;
                int v;
                double MaxDist = 0;
                if (! VQueues[i].Empty()) {
                    const int NId = VQueues[i].Top();  VQueues[i].Pop();
                    const double Dist = NIdDistH.GetDat(NId);
                    if (Dist == MxDist) { break; } // max distance limit reached
                    //const TNodeEDatNet<TInt, TInt>::TNodeI NodeI = Graph->GetNI(NId);
                    const TNodeEDatNet<TInt, TFlt>::TNodeI NodeI = Graph->GetNI(NId);
                    for (v = 0; v < NodeI.GetOutDeg(); v++) {  // out-links
                        const int DstNId = NodeI.GetOutNId(v);
                        if (! NIdDistH.IsKey(DstNId)) {
                            //std::cout << NId << " -> "<< DstNId << "\n" << std::endl;
                            NIdDistH.AddDat(DstNId, Dist+1);
                            MaxDist = std::max(MaxDist, Dist+1);
                            VQueues[i].Push(DstNId);
                            Vbfstree_nodes[i].Add(DstNId);
                        }
                        else{ // DstNId is already visited by some vertex (may be from cluster i, or other clusters.
                           
                            // Which cluster does DstNId belong to?
                            for(int j=(i+1)%N_startNodes,k=0; k<N_startNodes; k++,j=(j+1)%N_startNodes){
                                TIntV cl_nodes = Vbfstree_nodes[j];
                                // reverse iterate
                                //std::cout << " reverse iterate  \n";
                                int adj_clusterid = -1;
                                for(int m=0; m<cl_nodes.Len(); m++){
                                    // std::cout << cl_nodes[m]<< " ";
                                    if (cl_nodes[m] == DstNId) {
                                        adj_clusterid = j;
                                        break;
                                    }
                                }
                                if (adj_clusterid==i){} // can not be a cluster boundary
                                if(adj_clusterid != i && adj_clusterid!=-1){
                                    cl_boundary[i].insert(NId);
                                    cl_boundary[j].insert(DstNId);
//                                }
//                                // std::cout << i<< "adj "<<adj_clusterid<<"\n";
//                                if (adj_clusterid!=-1){
                                    TIntPr pair;
                                    if(i<adj_clusterid){
                                        pair.Val1 = i;
                                        pair.Val2 = adj_clusterid;
                                        //std::cout << "pair \n";
                                    }
                                    else{
                                        pair.Val1 = adj_clusterid;
                                        pair.Val2 = i;
                                    }
                                    if(smallest_intra_clust_dist.IsKey(pair)){
                                        double prev_dist = smallest_intra_clust_dist.GetDat(pair);
                                        //std::cout <<"after getdat\n";
                                        //auto v1 = Graph->GetNI(NId);
                                        //auto v2 = Graph->GetNI(DstNId);
                                        //std::cout << "boundary: "<<NId << " "<<DstNId<<"\n";
                                        double weight = Graph->GetEDat(NId,DstNId);
                                        smallest_intra_clust_dist.AddDat(pair,TMath::Mn(prev_dist,weight));
                                    }
                                    else
                                        smallest_intra_clust_dist.AddDat(pair,Graph->GetEDat(NId,DstNId));
                                    
                                }
                            }
                        }
                        
                    }
                }
                else{
                    empty_counter++;
                }
            }
            if (empty_counter == N_startNodes)
                break;
        }
    }
    double DoBfs(const int& StartNode, const bool& FollowOut, const bool& FollowIn, const int& TargetNId, const int& MxDist) {
        StartNId = StartNode;
        IAssert(Graph->IsNode(StartNId));
        //  const typename PGraph::TObj::TNodeI StartNodeI = Graph->GetNI(StartNode);
        //  IAssertR(StartNodeI.GetOutDeg() > 0, TStr::Fmt("No neighbors from start node %d.", StartNode));
        NIdDistH.Clr(false);  NIdDistH.AddDat(StartNId, 0.0);
        Queue.Clr(false);  Queue.Push(StartNId);
        int v;
        double MaxDist = 0.0;
        
        while (! Queue.Empty()) {
            const int NId = Queue.Top();  Queue.Pop();
            const double Dist = NIdDistH.GetDat(NId);
            if (Dist == MxDist) { break; } // max distance limit reached
            const TNodeEDatNet<TInt, TFlt>::TNodeI NodeI = Graph->GetNI(NId);
            if (FollowOut) { // out-links
                for (v = 0; v < NodeI.GetOutDeg(); v++) {  // out-links
                    const int DstNId = NodeI.GetOutNId(v);
                    if (! NIdDistH.IsKey(DstNId)) {
                        NIdDistH.AddDat(DstNId, Dist+Graph->GetEDat(NId,DstNId));
                        MaxDist = TMath::Mx(MaxDist, Dist+Graph->GetEDat(NId,DstNId));
                        if (DstNId == TargetNId) { return MaxDist; }
                        Queue.Push(DstNId);
                        //                      std::cout << NId << ":" << DstNId << ":"<< MaxDist << std::endl;
                    }
                }
            }
            /*
            if (FollowIn) { // in-links
                for (v = 0; v < NodeI.GetInDeg(); v++) {
                    const int DstNId = NodeI.GetInNId(v);
                    if (! NIdDistH.IsKey(DstNId)) {
                        NIdDistH.AddDat(DstNId, Dist+1);
                        MaxDist = TMath::Mx(MaxDist, Dist+1);
                        if (DstNId == TargetNId) { return MaxDist; }
                        Queue.Push(DstNId);
                    }
                }
            }*/
        }
        
        return MaxDist;
    }
    double GetHops(const int& SrcNId, const int& DstNId){
        TFlt Dist;
        if (SrcNId!=StartNId) { return -1; }
        if (! NIdDistH.IsKeyGetDat(DstNId, Dist)) { return -1; }
        return Dist.Val;
    }
};


//#//////////////////////////////////////////////
/// Breath-First-Search class.
/// The class is meant for executing many BFSs over a fixed graph. This means that the class can keep the hash tables and queues initialized between different calls of the DoBfs() function.
template<class PGraph>
class TBreathFS {
public:
  PGraph Graph;
  TSnapQueue<int> Queue;
  TInt StartNId;
  TIntH NIdDistH;
  TIntQV VQueues; // typedef TVec<TQQueue<TInt> > TIntQV;
  bool concurrent;
  TIntPrFltH smallest_intra_clust_dist;
  std::vector <TIntV> cl_boundary;
public:
    
    TBreathFS(const PGraph& GraphPt, const bool& concurrent, const bool& InitBigQ=true) :
    Graph(GraphPt), Queue(InitBigQ?Graph->GetNodes():1024), NIdDistH(InitBigQ?Graph->GetNodes():1024),concurrent(true) {
        // We keep only one NiDDistH hashtable. As a vertex is allowed to be the shortest to only one of the StartNodes.
        // But we need keep a queues for each start node.
        //std::cout << "oka"<<"\n";
    }
    
  TBreathFS(const PGraph& GraphPt, const bool& InitBigQ=true) :
    Graph(GraphPt), Queue(InitBigQ?Graph->GetNodes():1024), NIdDistH(InitBigQ?Graph->GetNodes():1024) { }
  /// Sets the graph to be used by the BFS to GraphPt and resets the data structures.
  void SetGraph(const PGraph& GraphPt);
  /// Concurrently perform BFS from a bunch of nodes, stops a path exploration when conflict happens, stores bfs tree nodes in i'th position of Vbfstree_nodes
  void DoBfs_concurrent(const TIntV& StartNodes, vector <TIntV> &Vbfstree_nodes, const int& MxDist=TInt::Mx);
  /// Performs BFS from node id StartNode for at maps MxDist steps by only following in-links (parameter FollowIn = true) and/or out-links (parameter FollowOut = true).
  int DoBfs(const int& StartNode, const bool& FollowOut, const bool& FollowIn, const int& TargetNId=-1, const int& MxDist=TInt::Mx);
  /// Same functionality as DoBfs with better performance.
  int DoBfsHybrid(const int& StartNode, const bool& FollowOut, const bool& FollowIn, const int& TargetNId=-1, const int& MxDist=TInt::Mx);
  /// Returns the number of nodes visited/reached by the BFS.
  int GetNVisited() const { return NIdDistH.Len(); }
  /// Returns the IDs of the nodes visited/reached by the BFS.
  void GetVisitedNIdV(TIntV& NIdV) const { NIdDistH.GetKeyV(NIdV); }
  /// Returns the shortst path distance between SrcNId and DistNId.
  /// Note you have to first call DoBFs(). SrcNId must be equal to StartNode, otherwise return value is -1.
  int GetHops(const int& SrcNId, const int& DstNId) const;
  /// Returns a random shortest path from SrcNId to DstNId.
  /// Note you have to first call DoBFs(). SrcNId must be equal to StartNode, otherwise return value is -1.
  int GetRndPath(const int& SrcNId, const int& DstNId, TIntV& PathNIdV) const;
    void DoBfs_concurrent(const TNodeEDatNet<TInt,TFlt>& net,vector <TIntV> &Vbfstree_nodes, const int& MxDist=TInt::Mx);

/* Private variables and functions for DoBfsHybrid */
private:
  int Stage; // 0, 2: top down, 1: bottom up
  static const unsigned int alpha = 100;
  static const unsigned int beta = 20;
  /* Private functions */
  bool TopDownStep(TIntV &NIdDistV, TIntV *Frontier, TIntV *NextFrontier, int& MaxDist, const int& TargetNId, const bool& FollowOut, const bool& FollowIn);
  bool BottomUpStep(TIntV &NIdDistV, TIntV *Frontier, TIntV *NextFrontier, int& MaxDist, const int& TargetNId, const bool& FollowOut, const bool& FollowIn);
};


///// Start - of - my - code
template <class PGraph>
void TBreathFS<PGraph>::DoBfs_concurrent(const TIntV& StartNodes, std::vector <TIntV> &Vbfstree_nodes, const int& MxDist){
    const int N_startNodes = StartNodes.Len();
    if(concurrent){
        VQueues.Reserve(N_startNodes);
    }
    NIdDistH.Clr(false);
    
    for (int i =0; i<N_startNodes; i++){
        auto StartNId = StartNodes[i];
        IAssert(Graph->IsNode(StartNId));
        //  const typename PGraph::TObj::TNodeI StartNodeI = Graph->GetNI(StartNode);
        //  IAssertR(StartNodeI.GetOutDeg() > 0, TStr::Fmt("No neighbors from start node %d.", StartNode));
        NIdDistH.AddDat(StartNId, 0);
        VQueues[i].Clr(false);  VQueues[i].Push(StartNId);
        Vbfstree_nodes[i].Add(StartNId);
        cl_boundary.push_back(TIntV());
    }
    while(true){
        int empty_counter = 0;
        for (int i = 0; i < N_startNodes ; i++){
//            std::cout << "iter " << i << ": "<<std::endl;
            int v, MaxDist = 0;
            if (! VQueues[i].Empty()) {
                const int NId = VQueues[i].Top();  VQueues[i].Pop();
                const int Dist = NIdDistH.GetDat(NId);
                if (Dist == MxDist) { break; } // max distance limit reached
                const typename PGraph::TObj::TNodeI NodeI = Graph->GetNI(NId);
                for (v = 0; v < NodeI.GetOutDeg(); v++) {  // out-links
                    const int DstNId = NodeI.GetOutNId(v);
                    if (! NIdDistH.IsKey(DstNId)) {
//                        std::cout << NId << " -> "<< DstNId << "\n" << std::endl;
                        NIdDistH.AddDat(DstNId, Dist+1);
                        MaxDist = TMath::Mx(MaxDist, Dist+1);
                        VQueues[i].Push(DstNId);
                        Vbfstree_nodes[i].Add(DstNId);
                    }
                    else{ // DstNId is already in some other cluster. meaning, NodeI is a boundary vertex.
                        cl_boundary[i].Add(NId);
                        // Which cluster does DstNId belong to?
                        for(int j=(i+1)%N_startNodes,k=0; k<N_startNodes; k++,j=(j+1)%N_startNodes){
                            TIntV cl_nodes = Vbfstree_nodes[j];
                            // reverse iterate
                            //std::cout << " reverse iterate  \n";
                            int adj_clusterid = -1;
                            for(int m=0; m<cl_nodes.Len() && j!=i ;m++){
                               // std::cout << cl_nodes[m]<< " ";
                                if (cl_nodes[m] == DstNId) break;
                                adj_clusterid = j;
                            }
                           // std::cout << i<< "adj "<<adj_clusterid<<"\n";
                            if (adj_clusterid!=-1){
                                TIntPr pair;
                                if(i<adj_clusterid){
                                    pair.Val1 = i;
                                    pair.Val2 = adj_clusterid;
                                    //std::cout << "pair \n";
                                }
                                else{
                                    pair.Val1 = adj_clusterid;
                                    pair.Val2 = i;
                                }
                                if(smallest_intra_clust_dist.IsKey(pair)){
                                    double prev_dist = smallest_intra_clust_dist.GetDat(pair);
                                    //std::cout <<"after getdat\n";
                                    smallest_intra_clust_dist.AddDat(pair,TMath::Mn(prev_dist,1.0));
                                }
                                else
                                    smallest_intra_clust_dist.AddDat(pair,1.0);
                                
                            }
                        }
                    }
                }
            }
            else{
                empty_counter++;
            }
        }
        if (empty_counter == N_startNodes)
            break;
    }
}
/////// End - of - my - code
template<class PGraph>
void TBreathFS<PGraph>::SetGraph(const PGraph& GraphPt) {
  Graph=GraphPt;
  const int N=GraphPt->GetNodes();
  if (Queue.Reserved() < N) { Queue.Gen(N); }
  if (NIdDistH.GetReservedKeyIds() < N) { NIdDistH.Gen(N); }
}

template<class PGraph>
int TBreathFS<PGraph>::DoBfs(const int& StartNode, const bool& FollowOut, const bool& FollowIn, const int& TargetNId, const int& MxDist) {
  StartNId = StartNode;
  IAssert(Graph->IsNode(StartNId));
//  const typename PGraph::TObj::TNodeI StartNodeI = Graph->GetNI(StartNode);
//  IAssertR(StartNodeI.GetOutDeg() > 0, TStr::Fmt("No neighbors from start node %d.", StartNode));
  NIdDistH.Clr(false);  NIdDistH.AddDat(StartNId, 0);
  Queue.Clr(false);  Queue.Push(StartNId);
  int v, MaxDist = 0;
  while (! Queue.Empty()) {
    const int NId = Queue.Top();  Queue.Pop();
    const int Dist = NIdDistH.GetDat(NId);
    if (Dist == MxDist) { break; } // max distance limit reached
    const typename PGraph::TObj::TNodeI NodeI = Graph->GetNI(NId);
    if (FollowOut) { // out-links
      for (v = 0; v < NodeI.GetOutDeg(); v++) {  // out-links
        const int DstNId = NodeI.GetOutNId(v);
        if (! NIdDistH.IsKey(DstNId)) {
          NIdDistH.AddDat(DstNId, Dist+1);
          MaxDist = TMath::Mx(MaxDist, Dist+1);
          if (DstNId == TargetNId) { return MaxDist; }
          Queue.Push(DstNId);
//                      std::cout << NId << ":" << DstNId << ":"<< MaxDist << std::endl;
        }
      }
    }
    if (FollowIn) { // in-links
      for (v = 0; v < NodeI.GetInDeg(); v++) {
        const int DstNId = NodeI.GetInNId(v);
        if (! NIdDistH.IsKey(DstNId)) {
          NIdDistH.AddDat(DstNId, Dist+1);
          MaxDist = TMath::Mx(MaxDist, Dist+1);
          if (DstNId == TargetNId) { return MaxDist; }
          Queue.Push(DstNId);
        }
      }
    }
  }
  return MaxDist;
}

template<class PGraph>
int TBreathFS<PGraph>::DoBfsHybrid(const int& StartNode, const bool& FollowOut, const bool& FollowIn, const int& TargetNId, const int& MxDist) {
  StartNId = StartNode;
  IAssert(Graph->IsNode(StartNId));
  if (TargetNId == StartNode) return 0;
  const typename PGraph::TObj::TNodeI StartNodeI = Graph->GetNI(StartNode);

  // Initialize vector
  TIntV NIdDistV(Graph->GetMxNId() + 1);
  for (int i = 0; i < NIdDistV.Len(); i++) {
    NIdDistV.SetVal(i, -1);
  }
  TIntV *Frontier = new TIntV(Graph->GetNodes(), 0);
  TIntV *NextFrontier = new TIntV(Graph->GetNodes(), 0);

  NIdDistV.SetVal(StartNId, 0);
  Frontier->Add(StartNId);
  Stage = 0;
  int MaxDist = -1;
  const unsigned int TotalNodes = Graph->GetNodes();
  unsigned int UnvisitedNodes = Graph->GetNodes();
  while (! Frontier->Empty()) {
    MaxDist += 1;
    NextFrontier->Clr(false);
    if (MaxDist == MxDist) { break; } // max distance limit reached

    UnvisitedNodes -= Frontier->Len();
    if (Stage == 0 && UnvisitedNodes / Frontier->Len() < alpha) {
      Stage = 1;
    } else if (Stage == 1 && TotalNodes / Frontier->Len() > beta) {
      Stage = 2;
    }

    // Top down or bottom up depending on stage
    bool targetFound = false;
    if (Stage == 0 || Stage == 2) {
      targetFound = TopDownStep(NIdDistV, Frontier, NextFrontier, MaxDist, TargetNId, FollowOut, FollowIn);
    } else {
      targetFound = BottomUpStep(NIdDistV, Frontier, NextFrontier, MaxDist, TargetNId, FollowOut, FollowIn);
    }
    if (targetFound) {
      MaxDist = NIdDistV[TargetNId];
      break;
    }

    // swap Frontier and NextFrontier
    TIntV *temp = Frontier;
    Frontier = NextFrontier;
    NextFrontier = temp;
  }

  delete Frontier;
  delete NextFrontier;
  // Transform vector to hash table
  NIdDistH.Clr(false);
  for (int NId = 0; NId < NIdDistV.Len(); NId++) {
    if (NIdDistV[NId] != -1) {
      NIdDistH.AddDat(NId, NIdDistV[NId]);
    }
  }
  return MaxDist;
}

template<class PGraph>
bool TBreathFS<PGraph>::TopDownStep(TIntV &NIdDistV, TIntV *Frontier, TIntV *NextFrontier, int& MaxDist, const int& TargetNId, const bool& FollowOut, const bool& FollowIn) {
  for (TIntV::TIter it = Frontier->BegI(); it != Frontier->EndI(); ++it) { // loop over frontier
    const int NId = *it;
    const int Dist = NIdDistV[NId];
    IAssert(Dist == MaxDist); // Must equal to MaxDist
    const typename PGraph::TObj::TNodeI NodeI = Graph->GetNI(NId);
    if (FollowOut) {
      for (int v = 0; v < NodeI.GetOutDeg(); v++) {
        const int NeighborNId = NodeI.GetOutNId(v);
        if (NIdDistV[NeighborNId] == -1) {
          NIdDistV.SetVal(NeighborNId, Dist+1);
          if (NeighborNId == TargetNId) return true;
          NextFrontier->Add(NeighborNId);
        }
      }
    }
    if (FollowIn) {
      for (int v = 0; v < NodeI.GetInDeg(); v++) {
        const int NeighborNId = NodeI.GetInNId(v);
        if (NIdDistV[NeighborNId] == -1) {
          NIdDistV.SetVal(NeighborNId, Dist+1);
          if (NeighborNId == TargetNId) return true;
          NextFrontier->Add(NeighborNId);
        }
      }
    }
  }
  return false;
}

template<class PGraph>
bool TBreathFS<PGraph>::BottomUpStep(TIntV &NIdDistV, TIntV *Frontier, TIntV *NextFrontier, int& MaxDist, const int& TargetNId, const bool& FollowOut, const bool& FollowIn) {
  for (typename PGraph::TObj::TNodeI NodeI = Graph->BegNI(); NodeI < Graph->EndNI(); NodeI++) {
    const int NId = NodeI.GetId();
    if (NIdDistV[NId] == -1) {
      if (FollowOut) {
        for (int v = 0; v < NodeI.GetInDeg(); v++) {
          const int ParentNId = NodeI.GetInNId(v);
          if (NIdDistV[ParentNId] == MaxDist) {
            NIdDistV[NId] = MaxDist + 1;
            if (NId == TargetNId) return true;
            NextFrontier->Add(NId);
            break;
          }
        }
      }
      if (FollowIn && NIdDistV[NId] == -1) {
        for (int v = 0; v < NodeI.GetOutDeg(); v++) {
          const int ParentNId = NodeI.GetOutNId(v);
          if (NIdDistV[ParentNId] == MaxDist) {
            NIdDistV[NId] = MaxDist + 1;
            if (NId == TargetNId) return true;
            NextFrontier->Add(NId);
            break;
          }
        }
      }
    }
  }
  return false;
}

template<class PGraph>
int TBreathFS<PGraph>::GetHops(const int& SrcNId, const int& DstNId) const {
  TInt Dist;
  if (SrcNId!=StartNId) { return -1; }
  if (! NIdDistH.IsKeyGetDat(DstNId, Dist)) { return std::numeric_limits<int>::max(); }
  return Dist.Val;
}

template<class PGraph>
int TBreathFS<PGraph>::GetRndPath(const int& SrcNId, const int& DstNId, TIntV& PathNIdV) const {
  PathNIdV.Clr(false);
  if (SrcNId!=StartNId || ! NIdDistH.IsKey(DstNId)) { return -1; }
  PathNIdV.Add(DstNId);
  TIntV CloserNIdV;
  int CurNId = DstNId;
  TInt CurDist, NextDist;
  while (CurNId != SrcNId) {
    typename PGraph::TObj::TNodeI NI = Graph->GetNI(CurNId);
    IAssert(NIdDistH.IsKeyGetDat(CurNId, CurDist));
    CloserNIdV.Clr(false);
    for (int e = 0; e < NI.GetDeg(); e++) {
      const int Next = NI.GetNbrNId(e);
      if (NIdDistH.IsKeyGetDat(Next, NextDist)) {
        if (NextDist == CurDist-1) { CloserNIdV.Add(Next); }
      }
    }
    IAssert(! CloserNIdV.Empty());
    CurNId = CloserNIdV[TInt::Rnd.GetUniDevInt(CloserNIdV.Len())];
    PathNIdV.Add(CurNId);
  }
  PathNIdV.Reverse();
  return PathNIdV.Len()-1;
}

/////////////////////////////////////////////////
// Implementation
namespace TSnap {
template <class PGraph>
std::pair<PUNGraph,TIntIntH> GetMinSpanTree(const PGraph& Graph, const int &StartNId){
    std::cout << "minimum spanning tree "<<StartNId<<"\n";
    typedef pair<double, int> iPair;
    PUNGraph Tree = TUNGraph::New();
    TIntIntH parent; // who is who's parent.
    int N_ = Graph->GetNodes();
    std::cout << N_ <<"\n";
    std::vector<bool> visited(N_,false);
    std::vector <double> key(N_,std::numeric_limits<double>::max());
    std::priority_queue<iPair,vector <iPair> , greater<iPair>> Q;
//    for (int i = 0; i < N_; i++){
//        key[i] = std::numeric_limits<double>::max();
//        visited[i] = false;
//        
//    }
//    parent.Clr(false);
    Q.push(make_pair(0.0,StartNId));
    key[StartNId] = 0.0;
//    parent[StartNId] = -1;
    parent.AddDat(StartNId,-1);
    while (!Q.empty()) {
        int u = Q.top().second;
//        std::cout << Q.top().first << " "<<Q.top().second << " \n";
        Q.pop();
        visited[u] = true;
        if(!Tree->IsNode(u))
            Tree->AddNode(u);
        auto NI = Graph->GetNI(u);
        for (int e = 0; e < NI.GetOutDeg(); e++) {
            const int nbr = NI.GetOutNId(e);
            int weight = 1.0;
            if (visited[nbr] == false && key[nbr] > weight) {
                key[nbr] = weight;
                Q.push(make_pair(weight, nbr));
                // Add an edge from u -> nbr
//                parent[nbr] = u;
                parent.AddDat(nbr,u);
                Tree->AddNode(nbr);
                Tree->AddEdge(nbr,u);
            }
            
        }
    }
    std::cout <<"returning \n";
    return make_pair(Tree,parent);
    }

template <class nodewType,class edgewType>
std::pair<PUNGraph,TIntIntH> GetMinSpanTree(const TPt<TNodeEDatNet<nodewType, edgewType>>& Graph, const int &StartNId){
        std::cout << "minimum spanning tree "<<StartNId<<"\n";
        typedef pair<double, int> iPair;
        PUNGraph Tree = TUNGraph::New();
        TIntIntH parent; // who is who's parent.
        int N_ = Graph->GetNodes();
        std::cout << N_ <<"\n";
        std::vector<bool> visited(N_,false);
        std::vector <double> key(N_,std::numeric_limits<double>::max());
        std::priority_queue<iPair,vector <iPair> , greater<iPair>> Q;
        Q.push(make_pair(0.0,StartNId));
        key[StartNId] = 0.0;
        parent.AddDat(StartNId,-1);
        while (!Q.empty()) {
            int u = Q.top().second;
            //        std::cout << Q.top().first << " "<<Q.top().second << " \n";
            Q.pop();
            visited[u] = true;
            if(!Tree->IsNode(u))
                Tree->AddNode(u);
            auto NI = Graph->GetNI(u);
            for (int e = 0; e < NI.GetOutDeg(); e++) {
                const int nbr = NI.GetOutNId(e);
                int weight = 1.0;
                if (visited[nbr] == false && key[nbr] > weight) {
                    key[nbr] = weight;
                    Q.push(make_pair(weight, nbr));
                    // Add an edge from u -> nbr
                    //                parent[nbr] = u;
                    parent.AddDat(nbr,u);
                    Tree->AddNode(nbr);
                    Tree->AddEdge(nbr,u);
                }
                
            }
        }
        std::cout <<"returning \n";
        return make_pair(Tree,parent);
    }
    // end -my -code
template <class PGraph>
PNGraph GetBfsTree(const PGraph& Graph, const int& StartNId, const bool& FollowOut, const bool& FollowIn) {
  TBreathFS<PGraph> BFS(Graph);
  BFS.DoBfs(StartNId, FollowOut, FollowIn, -1, TInt::Mx);
  PNGraph Tree = TNGraph::New();
  BFS.NIdDistH.SortByDat();
  for (int i = 0; i < BFS.NIdDistH.Len(); i++) {
    const int NId = BFS.NIdDistH.GetKey(i);
    const int Dist = BFS.NIdDistH[i];
    typename PGraph::TObj::TNodeI NI = Graph->GetNI(NId);
    if (!Tree->IsNode(NId)) {
      Tree->AddNode(NId);
    }
    if (FollowOut) {
      for (int e = 0; e < NI.GetInDeg(); e++) {
        const int Prev = NI.GetInNId(e);
        if (Tree->IsNode(Prev) && BFS.NIdDistH.GetDat(Prev)==Dist-1) {
          Tree->AddEdge(Prev, NId); }
      }
    }
    if (FollowIn) {
      for (int e = 0; e < NI.GetOutDeg(); e++) {
        const int Prev = NI.GetOutNId(e);
        if (Tree->IsNode(Prev) && BFS.NIdDistH.GetDat(Prev)==Dist-1) {
          Tree->AddEdge(Prev, NId); }
      }
    }
  }
  return Tree;
}

template <class PGraph>
int GetSubTreeSz(const PGraph& Graph, const int& StartNId, const bool& FollowOut, const bool& FollowIn, int& TreeSz, int& TreeDepth) {
  TBreathFS<PGraph> BFS(Graph);
  BFS.DoBfs(StartNId, FollowOut, FollowIn, -1, TInt::Mx);
  TreeSz = BFS.NIdDistH.Len();
  TreeDepth = 0;
  for (int i = 0; i < BFS.NIdDistH.Len(); i++) {
    TreeDepth = TMath::Mx(TreeDepth, BFS.NIdDistH[i].Val);
  }
  return TreeSz;
}

template <class PGraph>
int GetNodesAtHop(const PGraph& Graph, const int& StartNId, const int& Hop, TIntV& NIdV, const bool& IsDir) {
  TBreathFS<PGraph> BFS(Graph);
  BFS.DoBfs(StartNId, true, !IsDir, -1, Hop);
  NIdV.Clr(false);
  for (int i = 0; i < BFS.NIdDistH.Len(); i++) {
    if (BFS.NIdDistH[i] == Hop) {
      NIdV.Add(BFS.NIdDistH.GetKey(i)); }
  }
  return NIdV.Len();
}

template <class PGraph>
int GetNodesAtHops(const PGraph& Graph, const int& StartNId, TIntPrV& HopCntV, const bool& IsDir) {
  TBreathFS<PGraph> BFS(Graph);
  BFS.DoBfs(StartNId, true, !IsDir, -1, TInt::Mx);
  TIntH HopCntH;
  for (int i = 0; i < BFS.NIdDistH.Len(); i++) {
    HopCntH.AddDat(BFS.NIdDistH[i]) += 1;
  }
  HopCntH.GetKeyDatPrV(HopCntV);
  HopCntV.Sort();
  return HopCntV.Len();
}

template <class PGraph>
int GetShortPath(const PGraph& Graph, const int& SrcNId, TIntH& NIdToDistH, const bool& IsDir, const int& MaxDist) {
  TBreathFS<PGraph> BFS(Graph);
  BFS.DoBfs(SrcNId, true, ! IsDir, -1, MaxDist);
  NIdToDistH.Clr();
  NIdToDistH.Swap(BFS.NIdDistH);
  return NIdToDistH[NIdToDistH.Len()-1];
}

template <class PGraph>
int GetShortPath(const PGraph& Graph, const int& SrcNId, const int& DstNId, const bool& IsDir) {
  TBreathFS<PGraph> BFS(Graph);
  BFS.DoBfs(SrcNId, true, ! IsDir, DstNId, TInt::Mx);
  return BFS.GetHops(SrcNId, DstNId);
}

template <class PGraph>
int GetBfsFullDiam(const PGraph& Graph, const int& NTestNodes, const bool& IsDir) {
  int FullDiam;
  double EffDiam;
  GetBfsEffDiam(Graph, NTestNodes, IsDir, EffDiam, FullDiam);
  return FullDiam;
}

template <class PGraph>
double GetBfsEffDiam(const PGraph& Graph, const int& NTestNodes, const bool& IsDir) {
  int FullDiam;
  double EffDiam;
  GetBfsEffDiam(Graph, NTestNodes, IsDir, EffDiam, FullDiam);
  return EffDiam;
}

template <class PGraph>
double GetBfsEffDiam(const PGraph& Graph, const int& NTestNodes, const bool& IsDir, double& EffDiam, int& FullDiam) {
  double AvgDiam;
  EffDiam = -1;  FullDiam = -1;
  return GetBfsEffDiam(Graph, NTestNodes, IsDir, EffDiam, FullDiam, AvgDiam);
}

template <class PGraph>
double GetBfsEffDiam(const PGraph& Graph, const int& NTestNodes, const bool& IsDir, double& EffDiam, int& FullDiam, double& AvgSPL) {
  EffDiam = -1;  FullDiam = -1;  AvgSPL = -1;
  TIntFltH DistToCntH;
  TBreathFS<PGraph> BFS(Graph);
  // shotest paths
  TIntV NodeIdV;
  Graph->GetNIdV(NodeIdV);  NodeIdV.Shuffle(TInt::Rnd);
  for (int tries = 0; tries < TMath::Mn(NTestNodes, Graph->GetNodes()); tries++) {
    const int NId = NodeIdV[tries];
    BFS.DoBfs(NId, true, ! IsDir, -1, TInt::Mx);
    for (int i = 0; i < BFS.NIdDistH.Len(); i++) {
      DistToCntH.AddDat(BFS.NIdDistH[i]) += 1; }
  }
  TIntFltKdV DistNbrsPdfV;
  double SumPathL=0, PathCnt=0;
  for (int i = 0; i < DistToCntH.Len(); i++) {
    DistNbrsPdfV.Add(TIntFltKd(DistToCntH.GetKey(i), DistToCntH[i]));
    SumPathL += DistToCntH.GetKey(i) * DistToCntH[i];
    PathCnt += DistToCntH[i];
  }
  DistNbrsPdfV.Sort();
  EffDiam = TSnap::TSnapDetail::CalcEffDiamPdf(DistNbrsPdfV, 0.9); // effective diameter (90-th percentile)
  FullDiam = DistNbrsPdfV.Last().Key;                // approximate full diameter (max shortest path length over the sampled nodes)
  AvgSPL = SumPathL/PathCnt;                        // average shortest path length
  return EffDiam;
}

template <class PGraph>
double GetBfsEffDiamAll(const PGraph& Graph, const int& NTestNodes, const bool& IsDir, double& EffDiam, int& FullDiam, double& AvgSPL) {
  return GetBfsEffDiam(Graph, NTestNodes, IsDir, EffDiam, FullDiam, AvgSPL);
}

template <class PGraph>
double GetBfsEffDiam(const PGraph& Graph, const int& NTestNodes, const TIntV& SubGraphNIdV, const bool& IsDir, double& EffDiam, int& FullDiam) {
  EffDiam = -1;
  FullDiam = -1;

  TIntFltH DistToCntH;
  TBreathFS<PGraph> BFS(Graph);
  // shotest paths
  TIntV NodeIdV(SubGraphNIdV);  NodeIdV.Shuffle(TInt::Rnd);
  TInt Dist;
  for (int tries = 0; tries < TMath::Mn(NTestNodes, SubGraphNIdV.Len()); tries++) {
    const int NId = NodeIdV[tries];
    BFS.DoBfs(NId, true, ! IsDir, -1, TInt::Mx);
    for (int i = 0; i < SubGraphNIdV.Len(); i++) {
      if (BFS.NIdDistH.IsKeyGetDat(SubGraphNIdV[i], Dist)) {
        DistToCntH.AddDat(Dist) += 1;
      }
    }
  }
  TIntFltKdV DistNbrsPdfV;
  for (int i = 0; i < DistToCntH.Len(); i++) {
    DistNbrsPdfV.Add(TIntFltKd(DistToCntH.GetKey(i), DistToCntH[i]));
  }
  DistNbrsPdfV.Sort();
  EffDiam = TSnap::TSnapDetail::CalcEffDiamPdf(DistNbrsPdfV, 0.9);  // effective diameter (90-th percentile)
  FullDiam = DistNbrsPdfV.Last().Key;                 // approximate full diameter (max shortest path length over the sampled nodes)
  return EffDiam;                                     // average shortest path length
}

template <class PGraph>
int GetShortestDistances(const PGraph& Graph, const int& StartNId, const bool& FollowOut, const bool& FollowIn, TIntV& ShortestDists) {
  PSOut StdOut = TStdOut::New();
  int MxNId = Graph->GetMxNId();
  int NonNodeDepth = 2147483647; // INT_MAX
  int InfDepth = 2147483646; // INT_MAX - 1
  ShortestDists.Gen(MxNId);
  for (int NId = 0; NId < MxNId; NId++) {
    if (Graph->IsNode(NId)) { ShortestDists[NId] = InfDepth; }
    else { ShortestDists[NId] = NonNodeDepth; }
  }

  TIntV Vec1(MxNId, 0); // ensure enough capacity
  TIntV Vec2(MxNId, 0); // ensure enough capacity

  ShortestDists[StartNId] = 0;
  TIntV* PCurV = &Vec1;
  PCurV->Add(StartNId);
  TIntV* PNextV = &Vec2;
  int Depth = 0; // current depth
  while (!PCurV->Empty()) {
    Depth++; // increase depth
    for (int i = 0; i < PCurV->Len(); i++) {
      int NId = PCurV->GetVal(i);
      typename PGraph::TObj::TNodeI NI = Graph->GetNI(NId);
      for (int e = 0; e < NI.GetOutDeg(); e++) {
        const int OutNId = NI.GetOutNId(e);
        if (ShortestDists[OutNId].Val == InfDepth) {
          ShortestDists[OutNId] = Depth;
          PNextV->Add(OutNId);
        }
      }
    }
    // swap pointer, no copying
    TIntV* Tmp = PCurV;
    PCurV = PNextV;
    PNextV = Tmp;
    // clear next
    PNextV->Reduce(0); // reduce length, does not initialize new array
  }
  return Depth-1;
}

#ifdef USE_OPENMP
template <class PGraph>
int GetShortestDistancesMP2(const PGraph& Graph, const int& StartNId, const bool& FollowOut, const bool& FollowIn, TIntV& ShortestDists) {
  int MxNId = Graph->GetMxNId();
  int NonNodeDepth = 2147483647; // INT_MAX
  int InfDepth = 2147483646; // INT_MAX - 1
  ShortestDists.Gen(MxNId);
  #pragma omp parallel for schedule(dynamic,10000)
  for (int NId = 0; NId < MxNId; NId++) {
    if (Graph->IsNode(NId)) { ShortestDists[NId] = InfDepth; }
    else { ShortestDists[NId] = NonNodeDepth; }
  }

  TIntV Vec1(MxNId, 0); // ensure enough capacity
  TIntV Vec2(MxNId, 0); // ensure enough capacity

  ShortestDists[StartNId] = 0;
  TIntV* PCurV = &Vec1;
  PCurV->Add(StartNId);
  TIntV* PNextV = &Vec2;
  int Depth = 0; // current depth

  while (!PCurV->Empty()) {
    Depth++; // increase depth
    #pragma omp parallel for schedule(dynamic,10000)
    for (int i = 0; i < PCurV->Len(); i++) {
      int NId = PCurV->GetVal(i);
      typename PGraph::TObj::TNodeI NI = Graph->GetNI(NId);
      for (int e = 0; e < NI.GetOutDeg(); e++) {
        const int OutNId = NI.GetOutNId(e);
        if (__sync_bool_compare_and_swap(&(ShortestDists[OutNId].Val), InfDepth, Depth)) {
          PNextV->AddMP(OutNId);
        }
      }
    }
//      #pragma omp parallel for schedule(dynamic,10000)
//      for (int NId = 0; NId < MxNId; NId++) {
//        if (ShortestDists[NId] == InfDepth) {
//          typename PGraph::TObj::TNodeI NI = Graph->GetNI(NId);
//          for (int e = 0; e < NI.GetInDeg(); e++) {
//            const int InNId = NI.GetInNId(e);
//            if (ShortestDists[InNId] < Depth) {
//              ShortestDists[NId] = Depth;
//              PNextV->AddMP(NId);
//              break;
//            }
//          }
//        }
//      }
    // swap pointer, no copying
    TIntV* Tmp = PCurV;
    PCurV = PNextV;
    PNextV = Tmp;
    // clear next
    PNextV->Reduce(0); // reduce length, does not initialize new array
  }
  return Depth-1;
}
#endif // USE_OPENMP

} // namespace TSnap
