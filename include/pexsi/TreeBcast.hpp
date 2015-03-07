#ifndef _PEXSI_TREE_HPP_
#define _PEXSI_TREE_HPP_

#include "pexsi/environment.hpp"

#include <vector>

namespace PEXSI{

class TreeBcast{
  protected:
    Int myRoot_;
    MPI_Comm comm_;
    vector<Int> myDests_;
    Int myRank_;
  
    virtual void buildTree(Int * ranks, Int rank_cnt)=0;
  public:
    TreeBcast(const MPI_Comm & pComm, Int * ranks, Int rank_cnt){
      comm_ = pComm;
      MPI_Comm_rank(comm_,&myRank_);
      myRoot_ = -1; 
    }
    

    Int * GetDests(){ return &myDests_[0];}
    Int GetDest(Int i){ return myDests_[i];}
    Int GetDestCount(){ return myDests_.size();}
    Int GetRoot(){ return myRoot_;}
    void ForwardMessage( char * data, size_t size, int tag, MPI_Request * requests ){
                  for( Int idxRecv = 0; idxRecv < myDests_.size(); ++idxRecv ){
                    Int iProcRow = myDests_[idxRecv];
                    // Use Isend to send to multiple targets
                    MPI_Isend( data, size, MPI_BYTE, 
                        iProcRow, tag,comm_, &requests[2*iProcRow+1] );
                  } // for (iProcRow)
    }

};


class BTreeBcast: public TreeBcast{
  protected:
  virtual void buildTree(Int * ranks, Int rank_cnt){
          Int numLevel = floor(log2(rank_cnt));
          Int numRoots = 0;
          for(Int level=0;level<numLevel;++level){
            numRoots = std::min( rank_cnt, numRoots + (Int)pow(2,level));
            Int numNextRoots = std::min(rank_cnt,numRoots + (Int)pow(2,(level+1)));
            Int numReceivers = numNextRoots - numRoots;
            for(Int ip = 0; ip<numRoots;++ip){
              Int p = ranks[ip];
              for(Int ir = ip; ir<numReceivers;ir+=numRoots){
                Int r = ranks[numRoots+ir];
                if(r==myRank_){
                  myRoot_ = p;
                }

                if(p==myRank_){
                  myDests_.push_back(r);
                }
              }
            }
          }
  }
  public:
    BTreeBcast(const MPI_Comm & pComm, Int * ranks, Int rank_cnt):TreeBcast(pComm,ranks,rank_cnt){
      //build the binary tree;
      buildTree(ranks,rank_cnt);
    }


};






}

#endif
