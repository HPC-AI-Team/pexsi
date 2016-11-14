#ifndef _PEXSI_REDUCE_TREE_V2_HPP_
#define _PEXSI_REDUCE_TREE_V2_HPP_

#include "pexsi/environment.hpp"
#include "pexsi/timer.h"
#include "pexsi/TreeBcast_v2.hpp"

#include <vector>
#include <map>
#include <algorithm>
#include <string>
//#include <random>



namespace PEXSI{



  template< typename T>
    class TreeReduce_v2: public TreeBcast_v2<T>{
      protected:
        bool isAllocated_;
        Int numRecvPosted_;

      public:
        static TreeReduce_v2<T> * Create(const MPI_Comm & pComm, Int * ranks, Int rank_cnt, Int msgSize,double rseed);

        TreeReduce_v2(const MPI_Comm & pComm, Int * ranks, Int rank_cnt, Int msgSize);
        TreeReduce_v2(const TreeReduce_v2 & Tree);
        virtual ~TreeReduce_v2();
        virtual TreeReduce_v2 * clone() const = 0; 
        virtual void Copy(const TreeReduce_v2 & Tree);
        virtual void Reset();


        bool IsAllocated(){return isAllocated_;}
        virtual inline Int GetNumMsgToSend(){return this->myRank_==this->myRoot_?0:1;}
        virtual inline Int GetNumMsgToRecv(){return GetDestCount();}


        virtual void AllocRecvBuffers();
        


        virtual void SetLocalBuffer(T * locBuffer);
        virtual T * GetLocalBuffer();



        //async wait and forward
        virtual bool Progress();
        

        //  void CopyLocalBuffer(T* destBuffer){
        //    std::copy((char*)myData_,(char*)myData_+GetMsgSize(),(char*)destBuffer);
        //  }



      protected:
        virtual void reduce( Int idxRecv, Int idReq);
        virtual void forwardMessage();
        virtual void postRecv();
        virtual bool IsDataReceived();
        virtual bool isMessageForwarded();




    };


template< typename T>
class FTreeReduce_v2: public TreeReduce_v2<T>{
protected:
  virtual void buildTree(Int * ranks, Int rank_cnt){

    Int idxStart = 0;
    Int idxEnd = rank_cnt;

    this->myRoot_ = ranks[0];

    if(this->myRank_==this->myRoot_){
      this->myDests_.insert(this->myDests_.end(),&ranks[1],&ranks[0]+rank_cnt);
    }

#if (defined(REDUCE_VERBOSE))
    statusOFS<<"My root is "<<this->myRoot_<<std::endl;
    statusOFS<<"My dests are ";
    for(int i =0;i<this->myDests_.size();++i){statusOFS<<this->myDests_[i]<<" ";}
    statusOFS<<std::endl;
#endif
  }



public:
  FTreeReduce_v2(const MPI_Comm & pComm, Int * ranks, Int rank_cnt, Int msgSize):TreeReduce_v2<T>(pComm, ranks, rank_cnt, msgSize){
    buildTree(ranks,rank_cnt);
  }

  virtual void postRecv()
  {
    if(this->isAllocated_ && this->GetDestCount()>this->numRecvPosted_){
      MPI_Irecv( (char*)this->recvDataPtrs_[0], this->msgSize_, this->type_, 
          MPI_ANY_SOURCE, this->tag_,this->comm_, &this->recvRequests_[0] );
      this->recvPostedCount_++;
    }
  }

  virtual void AllocRecvBuffers(){
    this->recvDataPtrs_.assign(1,NULL);
    this->recvTempBuffer_.resize(this->msgSize_);

    this->recvDataPtrs_[0] = (T*)&(this->recvTempBuffer_[0]);

    this->recvRequests_.assign(1,MPI_REQUEST_NULL);
    this->recvStatuses_.resize(1);
    this->recvDoneIdx_.resize(1);
    this->sendRequests_.assign(1,MPI_REQUEST_NULL);

    this->isAllocated_ = true;
  }






  virtual bool Progress(){

          bool retVal = false;
          if(done_){
            retVal = true;
          }
          else{

            if(!isAllocated_){
              AllocRecvBuffers();
            }

            if(isAllocated_){
              if(myRank_==myRoot_ && isAllocated_){
                isReady_=true;
              }

              if(isReady_){
                if(IsDataReceived()){


                  //free the unnecessary arrays
                  this->recvTempBuffer_.clear();
                  this->recvRequests_.clear();
                  this->recvStatuses_.clear();
                  this->recvDoneIdx_.clear();

                  if(isMessageForwarded()){
                    retVal = true;
                  }
                }
                else if(this->recvPostedCount_<GetDestCount()){
                  //TODO check this
                  if(this->recvPostedCount_==this->recvCount_){
                    this->postRecv();
                  }
                }
              }
            }
          }

          if(retVal){
            this->done_ = retVal;
            //TODO do some smart cleanup
          }
          return retVal;
  }


  virtual FTreeReduce_v2<T> * clone() const{
    FTreeReduce_v2<T> * out = new FTreeReduce_v2<T>(*this);
    return out;
  }



};



template< typename T>
class BTreeReduce_v2: public TreeReduce_v2<T>{
protected:
  virtual void buildTree(Int * ranks, Int rank_cnt){
    Int idxStart = 0;
    Int idxEnd = rank_cnt;



    Int prevRoot = ranks[0];
    while(idxStart<idxEnd){
      Int curRoot = ranks[idxStart];
      Int listSize = idxEnd - idxStart;

      if(listSize == 1){
        if(curRoot == this->myRank_){
          this->myRoot_ = prevRoot;
          break;
        }
      }
      else{
        Int halfList = floor(ceil(double(listSize) / 2.0));
        Int idxStartL = idxStart+1;
        Int idxStartH = idxStart+halfList;

        if(curRoot == this->myRank_){
          if ((idxEnd - idxStartH) > 0 && (idxStartH - idxStartL)>0){
            Int childL = ranks[idxStartL];
            Int childR = ranks[idxStartH];

            this->myDests_.push_back(childL);
            this->myDests_.push_back(childR);
          }
          else if ((idxEnd - idxStartH) > 0){
            Int childR = ranks[idxStartH];
            this->myDests_.push_back(childR);
          }
          else{
            Int childL = ranks[idxStartL];
            this->myDests_.push_back(childL);
          }
          this->myRoot_ = prevRoot;
          break;
        } 

        if( this->myRank_ < ranks[idxStartH]){
          idxStart = idxStartL;
          idxEnd = idxStartH;
        }
        else{
          idxStart = idxStartH;
        }
        prevRoot = curRoot;
      }

    }

#if (defined(REDUCE_VERBOSE))
    statusOFS<<"My root is "<<this->myRoot_<<std::endl;
    statusOFS<<"My dests are ";
    for(int i =0;i<this->myDests_.size();++i){statusOFS<<this->myDests_[i]<<" ";}
    statusOFS<<std::endl;
#endif
  }
public:
  BTreeReduce_v2(const MPI_Comm & pComm, Int * ranks, Int rank_cnt, Int msgSize):TreeReduce_v2<T>(pComm, ranks, rank_cnt, msgSize){
    buildTree(ranks,rank_cnt);
  }

  virtual BTreeReduce_v2<T> * clone() const{
    BTreeReduce_v2<T> * out = new BTreeReduce_v2<T>(*this);
    return out;
  }
};


template< typename T>
class ModBTreeReduce_v2: public TreeReduce_v2<T>{
protected:
  double rseed_;
  virtual void buildTree(Int * ranks, Int rank_cnt){

    Int idxStart = 0;
    Int idxEnd = rank_cnt;

    //sort the ranks with the modulo like operation
    if(rank_cnt>1){
      //generate a random position in [1 .. rand_cnt]
      //Int new_idx = (int)((rand()+1.0) * (double)rank_cnt / ((double)RAND_MAX+1.0));
      //srand(ranks[0]+rank_cnt);
      //Int new_idx = rseed_%(rank_cnt-1)+1;

      //Int new_idx = (int)((rank_cnt - 0) * ( (double)this->rseed_ / (double)RAND_MAX ) + 0);// (this->rseed_)%(rank_cnt-1)+1;
      //Int new_idx = (Int)rseed_ % (rank_cnt - 1) + 1; 
//      Int new_idx = (int)((rank_cnt - 0) * ( (double)this->rseed_ / (double)RAND_MAX ) + 0);// (this->rseed_)%(rank_cnt-1)+1;
      Int new_idx = (int)(this->rseed_)%(rank_cnt-1)+1;

      Int * new_start = &ranks[new_idx];
      //        for(int i =0;i<rank_cnt;++i){statusOFS<<ranks[i]<<" ";} statusOFS<<std::endl;

      //        Int * new_start = std::lower_bound(&ranks[1],&ranks[0]+rank_cnt,ranks[0]);
      //just swap the two chunks   r[0] | r[1] --- r[new_start-1] | r[new_start] --- r[end]
      // becomes                   r[0] | r[new_start] --- r[end] | r[1] --- r[new_start-1] 
      std::rotate(&ranks[1], new_start, &ranks[0]+rank_cnt);
      //        for(int i =0;i<rank_cnt;++i){statusOFS<<ranks[i]<<" ";} statusOFS<<std::endl;
    }

    Int prevRoot = ranks[0];
    while(idxStart<idxEnd){
      Int curRoot = ranks[idxStart];
      Int listSize = idxEnd - idxStart;

      if(listSize == 1){
        if(curRoot == this->myRank_){
          this->myRoot_ = prevRoot;
          break;
        }
      }
      else{
        Int halfList = floor(ceil(double(listSize) / 2.0));
        Int idxStartL = idxStart+1;
        Int idxStartH = idxStart+halfList;

        if(curRoot == this->myRank_){
          if ((idxEnd - idxStartH) > 0 && (idxStartH - idxStartL)>0){
            Int childL = ranks[idxStartL];
            Int childR = ranks[idxStartH];

            this->myDests_.push_back(childL);
            this->myDests_.push_back(childR);
          }
          else if ((idxEnd - idxStartH) > 0){
            Int childR = ranks[idxStartH];
            this->myDests_.push_back(childR);
          }
          else{
            Int childL = ranks[idxStartL];
            this->myDests_.push_back(childL);
          }
          this->myRoot_ = prevRoot;
          break;
        } 

        //not true anymore ?
        //first half to 
        TIMER_START(FIND_RANK);
        Int * pos = std::find(&ranks[idxStartL], &ranks[idxStartH], this->myRank_);
        TIMER_STOP(FIND_RANK);
        if( pos != &ranks[idxStartH]){
          idxStart = idxStartL;
          idxEnd = idxStartH;
        }
        else{
          idxStart = idxStartH;
        }
        prevRoot = curRoot;
      }

    }

#if (defined(REDUCE_VERBOSE)) 
    statusOFS<<"My root is "<<this->myRoot_<<std::endl;
    statusOFS<<"My dests are ";
    for(int i =0;i<this->myDests_.size();++i){statusOFS<<this->myDests_[i]<<" ";}
    statusOFS<<std::endl;
#endif
  }
public:
  ModBTreeReduce_v2(const MPI_Comm & pComm, Int * ranks, Int rank_cnt, Int msgSize, double rseed):TreeReduce_v2<T>(pComm, ranks, rank_cnt, msgSize){
    this->rseed_ = rseed;
    buildTree(ranks,rank_cnt);
  }

  virtual void Copy(const ModBTreeReduce_v2<T> & Tree){
    ((TreeReduce_v2<T>*)this)->Copy(*((const TreeReduce_v2<T>*)&Tree));
    this->rseed_ = Tree.rseed_;
  }

  virtual ModBTreeReduce_v2<T> * clone() const{
    ModBTreeReduce_v2<T> * out = new ModBTreeReduce_v2<T>(*this);
    return out;
  }

};


template< typename T>
class PalmTreeReduce_v2: public TreeReduce_v2<T>{
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
          if(r==this->myRank_){
            this->myRoot_ = p;
          }

          if(p==this->myRank_){
            this->myDests_.push_back(r);
          }
        }
      }
    }

#if (defined(BCAST_VERBOSE))
    statusOFS<<"My root is "<<this->myRoot_<<std::endl;
    statusOFS<<"My dests are ";
    for(int i =0;i<this->myDests_.size();++i){statusOFS<<this->myDests_[i]<<" ";}
    statusOFS<<std::endl;
#endif
  }



public:
  PalmTreeReduce_v2(const MPI_Comm & pComm, Int * ranks, Int rank_cnt, Int msgSize):TreeReduce_v2<T>(pComm,ranks,rank_cnt,msgSize){
    //build the binary tree;
    buildTree(ranks,rank_cnt);
  }

  virtual PalmTreeReduce_v2<T> * clone() const{
    PalmTreeReduce_v2<T> * out = new PalmTreeReduce_v2<T>(*this);
    return out;
  }

};













}//namespace PEXSI

#include "pexsi/TreeReduce_v2_impl.hpp"
#endif
