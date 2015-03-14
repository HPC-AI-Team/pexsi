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
    Int msgSize_;
  
    virtual void buildTree(Int * ranks, Int rank_cnt)=0;
  public:
    TreeBcast(const MPI_Comm & pComm, Int * ranks, Int rank_cnt,Int msgSize){
      comm_ = pComm;
      MPI_Comm_rank(comm_,&myRank_);
      myRoot_ = -1; 
      msgSize_ = msgSize;
    }
    

    Int * GetDests(){ return &myDests_[0];}
    Int GetDest(Int i){ return myDests_[i];}
    Int GetDestCount(){ return myDests_.size();}
    Int GetRoot(){ return myRoot_;}
    Int GetMsgSize(){ return msgSize_;}

    void ForwardMessage( char * data, size_t size, int tag, MPI_Request * requests ){
                  for( Int idxRecv = 0; idxRecv < myDests_.size(); ++idxRecv ){
                    Int iProc = myDests_[idxRecv];
                    // Use Isend to send to multiple targets
                    MPI_Isend( data, size, MPI_BYTE, 
                        iProc, tag,comm_, &requests[2*iProc+1] );
                  } // for (iProc)
    }

    void ForwardMessage2( char * data, size_t size, int tag, MPI_Request * requests ){
                  for( Int idxRecv = 0; idxRecv < myDests_.size(); ++idxRecv ){
                    Int iProc = myDests_[idxRecv];
                    // Use Isend to send to multiple targets
                    MPI_Isend( data, size, MPI_BYTE, 
                        iProc, tag,comm_, &requests[iProc] );
                  } // for (iProc)
    }

    template< typename T>
    void SumAndForwardMessage2( T * mydata, T * data, Int size, int tag, MPI_Request * requests ){
                  //add thing to my data
                  blas::Axpy(size, ONE<T>(), data, 1, mydata, 1 );

                  for( Int idxRecv = 0; idxRecv < myDests_.size(); ++idxRecv ){
                    Int iProc = myDests_[idxRecv];
                    // Use Isend to send to multiple targets
                    MPI_Isend( (char*)mydata, size*sizeof(T), MPI_BYTE, 
                        iProc, tag,comm_, &requests[iProc] );
                  } // for (iProc)
    }
};


class BTreeBcast: public TreeBcast{
  protected:
////  virtual void buildTree(Int * ranks, Int rank_cnt){
////          Int numLevel = floor(log2(rank_cnt));
////          Int numRoots = 0;
////          for(Int level=0;level<numLevel;++level){
////            numRoots = std::min( rank_cnt, numRoots + (Int)pow(2,level));
////            Int numNextRoots = std::min(rank_cnt,numRoots + (Int)pow(2,(level+1)));
////            Int numReceivers = numNextRoots - numRoots;
////            for(Int ip = 0; ip<numRoots;++ip){
////              Int p = ranks[ip];
////              for(Int ir = ip; ir<numReceivers;ir+=numRoots){
////                Int r = ranks[numRoots+ir];
////                if(r==myRank_){
////                  myRoot_ = p;
////                }
////
////                if(p==myRank_){
////                  myDests_.push_back(r);
////                }
////              }
////            }
////          }
////  }
////
    virtual void buildTree(Int * ranks, Int rank_cnt){

      Int idxStart = 0;
      Int idxEnd = rank_cnt;



      Int prevRoot = ranks[0];
      while(idxStart<idxEnd){
        Int curRoot = ranks[idxStart];
        Int listSize = idxEnd - idxStart;

        if(listSize == 1){
          if(curRoot == myRank_){
            myRoot_ = prevRoot;
            break;
          }
        }
        else{
          Int halfList = floor(ceil(double(listSize) / 2.0));
          Int idxStartL = idxStart+1;
          Int idxStartH = idxStart+halfList;

          if(curRoot == myRank_){
            if ((idxEnd - idxStartH) > 0 && (idxStartH - idxStartL)>0){
              Int childL = ranks[idxStartL];
              Int childR = ranks[idxStartH];

              myDests_.push_back(childL);
              myDests_.push_back(childR);
            }
            else if ((idxEnd - idxStartH) > 0){
              Int childR = ranks[idxStartH];
              myDests_.push_back(childR);
            }
            else{
              Int childL = ranks[idxStartL];
              myDests_.push_back(childL);
            }
            myRoot_ = prevRoot;
            break;
          } 

          if( myRank_ < ranks[idxStartH]){
            idxStart = idxStartL;
            idxEnd = idxStartH;
          }
          else{
            idxStart = idxStartH;
          }
          prevRoot = curRoot;
        }

      }

#if ( _DEBUGlevel_ >= 1 )
      statusOFS<<"My root is "<<myRoot_<<std::endl;
      statusOFS<<"My dests are ";
      for(int i =0;i<myDests_.size();++i){statusOFS<<myDests_[i]<<" ";}
      statusOFS<<std::endl;
#endif
    }



  public:
    BTreeBcast(const MPI_Comm & pComm, Int * ranks, Int rank_cnt, Int msgSize):TreeBcast(pComm,ranks,rank_cnt,msgSize){
      //build the binary tree;
      buildTree(ranks,rank_cnt);
    }


};






}

#endif
