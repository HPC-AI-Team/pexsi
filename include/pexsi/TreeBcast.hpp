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
    bool isReady_;
    Int mainRoot_;
    Int tag_;
    Int numRecv_;




  
    virtual void buildTree(Int * ranks, Int rank_cnt)=0;
  public:
    TreeBcast(const MPI_Comm & pComm, Int * ranks, Int rank_cnt,Int msgSize){
      comm_ = pComm;
      MPI_Comm_rank(comm_,&myRank_);
      myRoot_ = -1; 
      msgSize_ = msgSize;

      numRecv_ = 0;
      tag_=-1;
      mainRoot_=ranks[0];
      isReady_ = false;
    }
    

    virtual inline Int GetNumRecvMsg(){return numRecv_;}
    virtual inline Int GetNumMsgToRecv(){return 1;}
    inline void SetDataReady(bool rdy){ isReady_=rdy;}
    inline void SetTag(Int tag){ tag_ = tag;}


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






template< typename T>
class TreeReduce: public TreeBcast{
  protected:

    std::vector<char> myLocalBuffer_;
    T * myData_;
    MPI_Request sendRequest_;

    std::vector<char> myRecvBuffers_;
    std::vector<T *> remoteData_;
    std::vector<MPI_Request> myRequests_;
    std::vector<MPI_Status> myStatuses_;
    std::vector<int> recvIdx_;

  public:
    TreeReduce(const MPI_Comm & pComm, Int * ranks, Int rank_cnt, Int msgSize):TreeBcast(pComm,ranks,rank_cnt,msgSize){
      myData_ = NULL;
      sendRequest_ = MPI_REQUEST_NULL;

    }

    virtual inline Int GetNumMsgToRecv(){return GetDestCount();}

    void AllocRecvBuffers(){
      remoteData_.resize(GetDestCount(),NULL);
      myRecvBuffers_.resize(GetDestCount()*msgSize_);
      for( Int idxRecv = 0; idxRecv < GetDestCount(); ++idxRecv ){
        remoteData_[idxRecv] = (T*)&myRecvBuffers_[idxRecv*msgSize_];
      }

      myRequests_.resize(GetDestCount(),MPI_REQUEST_NULL);
      myStatuses_.resize(GetDestCount());
      recvIdx_.resize(GetDestCount(),-1);

      sendRequest_ = MPI_REQUEST_NULL;
    }

    void CleanupBuffers(){
      {
        std::vector<T *> empty;
        empty.swap(remoteData_);
      }
      {
        std::vector<char> empty;
        empty.swap(myLocalBuffer_);
      }
      {
        std::vector<char> empty;
        empty.swap(myRecvBuffers_);
      }
      {
        std::vector<MPI_Request> empty;
        empty.swap(myRequests_);
      }
      {
        std::vector<MPI_Status> empty;
        empty.swap(myStatuses_);
      }
      {
        std::vector<int> empty;
        empty.swap(recvIdx_);
      }
    }


    void SetLocalBuffer(T * locBuffer){
      //if(myData_!=NULL){
      //  throw std::
      //}
      myData_ = locBuffer;
    }

    inline bool AccumulationDone(){
      return isReady_ && (numRecv_ == GetDestCount());
    }


    inline bool IsDone(){
      bool retVal = AccumulationDone();
      if(retVal && isReady_){
        //test the send request
        int flag = 0;
        MPI_Test(&sendRequest_,&flag,MPI_STATUS_IGNORE);
        retVal = flag==1;
      }
      return retVal;
    }

    //async wait and forward
    bool Progress(){
      bool retVal = AccumulationDone();
      if(isReady_ && !retVal){
        //mpi_test_some on my requests
        int recvCount = -1;
        MPI_Testsome(myRequests_.size(),&myRequests_[0],&recvCount,&recvIdx_[0],&myStatuses_[0]);
        //if something has been received, accumulate and potentially forward it
        for(Int i = 0;i<recvCount;++i ){
          Int idx = recvIdx_[i];
          Int size = 0;
          MPI_Get_count(&myStatuses_[i], MPI_BYTE, &size);
          if(size>0){


            //If myData is 0, allocate to the size of what has been received
            if(myData_==NULL){
              myLocalBuffer_.resize(size);
              myData_ = (T*)&myLocalBuffer_[0];
              std::fill(myData_,myData_+size/sizeof(T),ZERO<T>());
            }

            Reduce(idx);
          }

          MPI_Request_free(&myRequests_[idx]);
        }
      }
      else if (isReady_ && sendRequest_ == MPI_REQUEST_NULL){
        //Forward
        Forward();
        retVal = false;
      }
      else{
        retVal = IsDone();
      }

      return retVal;
    }

    //blocking wait
    void Wait(){
      while(!Progress());
    }

    T * GetLocalBuffer(){
       return myData_;
    }



    void CopyLocalBuffer(T* destBuffer){
       std::copy((char*)myData_,(char*)myData_+GetMsgSize(),destBuffer);
    }


    void PostAllRecv()
    {
      for( Int idxRecv = 0; idxRecv < myDests_.size(); ++idxRecv ){
        Int iProc = myDests_[idxRecv];
        MPI_Irecv( (char*)remoteData_[idxRecv], msgSize_*sizeof(T), MPI_BYTE, 
            iProc, tag_,comm_, &myRequests_[iProc] );
      } // for (iProc)
    }


    protected:
    void Reduce( Int idxRecv ){
      //add thing to my data
      blas::Axpy(msgSize_/sizeof(T), ONE<T>(), remoteData_[idxRecv], 1, myData_, 1 );
    }

    void Forward(){ 
      //forward to my root if I have reseived everything
      Int iProc = myRoot_;
      // Use Isend to send to multiple targets
      if(myData_==NULL){
        MPI_Isend( NULL, 0, MPI_BYTE, 
            iProc, tag_,comm_, &sendRequest_ );
      }
      else{
        MPI_Isend( (char*)myData_, msgSize_, MPI_BYTE, 
            iProc, tag_,comm_, &sendRequest_ );
      }

      //PROFILE_COMM(MYPROC(this->grid_),PNUM(MYROW(this->grid_),PCOL(snode.Index,this->grid_),this->grid_),IDX_TO_TAG(snode.Index,SELINV_TAG_L_REDUCE),redLTree->GetMsgSize());
    }

};



template< typename T>
class BTreeReduce: public TreeReduce<T>,BTreeBcast{
    public:
    BTreeReduce(const MPI_Comm & pComm, Int * ranks, Int rank_cnt, Int msgSize):BTreeBcast(pComm,ranks,rank_cnt,msgSize),TreeReduce(pComm, ranks, rank_cnt, msgSize){
    }
};





}

#endif
