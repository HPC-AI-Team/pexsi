#ifndef _PEXSI_REDUCE_TREE_IMPL_V2_HPP_
#define _PEXSI_REDUCE_TREE_IMPL_V2_HPP_
namespace PEXSI{
  template<typename T>
    TreeReduce_v2<T>::TreeReduce_v2(const MPI_Comm & pComm, Int * ranks, Int rank_cnt, Int msgSize):TreeBcast_v2<T>(pComm,ranks,rank_cnt,msgSize){
      this->sendDataPtrs_.assign(1,NULL);
      this->sendRequests_.assign(1,MPI_REQUEST_NULL);
      isAllocated_=false;
    }




  template<typename T>
    TreeReduce_v2<T>::TreeReduce_v2(const TreeReduce_v2<T> & Tree){
      this->Copy(Tree);
    }
  template<typename T>
    TreeReduce_v2<T>::~TreeReduce_v2(){
      cleanupBuffers();
    }

  template<typename T>
    inline void TreeReduce_v2<T>::Copy(const TreeReduce_v2<T> & Tree){
      ((TreeBcast_v2<T>*)this)->Copy(*(const TreeBcast_v2<T>*)&Tree);

      this->sendDataPtrs_.assign(1,NULL);
      this->sendRequests_.assign(1,MPI_REQUEST_NULL);
      this->isAllocated_= Tree.isAllocated_;

      cleanupBuffers();
    }



  template<typename T>
    inline void TreeReduce_v2<T>::postRecv(){
      if(this->GetDestCount()>this->recvPostedCount_){
        for( Int idxRecv = 0; idxRecv < myDests_.size(); ++idxRecv ){
          Int iProc = myDests_[idxRecv];
          MPI_Irecv( (char*)this->recvDataPtrs_[idxRecv], this->msgSize_, this->type_, 
              iProc, this->tag_,this->comm_, &this->recvRequests_[idxRecv] );
          this->recvPostedCount_++;
        } // for (iProc)
      }
    }




  template<typename T>
    inline void TreeReduce_v2<T>::reduce( Int idxRecv, Int idReq){
      //add thing to my data
      blas::Axpy(msgSize_, ONE<T>(), this->recvDataPtrs_[idxRecv], 1, this->sendDataPtrs_[0], 1 );
    }

  template<typename T>
    inline void TreeReduce_v2<T>::forwardMessage(){ 
      if(isReady_){
        if(this->myRank_!=this->myRoot_){
          //forward to my root if I have reseived everything
          Int iProc = myRoot_;
          // Use Isend to send to multiple targets
          if(this->sendDataPtrs_.size()<1){
            this->sendDataPtrs_.assign(1,NULL);
          }

          int msgsz = this->sendDataPtrs_[0]==NULL?0:this->msgSize_;

          MPI_Isend((char*)this->sendDataPtrs_[0], msgsz, this->type_, 
              iProc, this->tag_,this->comm_, &this->sendRequests_[0] );
          this->sendPostedCount_++;
#ifdef COMM_PROFILE
          PROFILE_COMM(myGRank_,myGRoot_,tag_,msgsz);
#endif

#if ( _DEBUGlevel_ >= 1 ) || defined(REDUCE_VERBOSE)
          statusOFS<<myRank_<<" FWD to "<<iProc<<" on tag "<<this->tag_<<std::endl;
#endif
        }
        this->fwded_ = true;
      }
    }

  template< typename T> 
    inline bool TreeReduce_v2<T>::IsDataReceived(){
      bool retVal = false;
      if(this->isReady_){
        if(this->recvCount_== this->GetDestCount()){
          retVal = true;
        }
        else{
          //mpi_test_some on recvRequests_
          int recvCount = -1;
          int reqCnt = this->recvRequests_.size();//this->recvPostedCount_-this->recvCount_;//GetDestCount();
//          assert(reqCnt <= this->recvRequests_.size());

          MPI_Testsome(reqCnt,&this->recvRequests_[0],&recvCount,&this->recvDoneIdx_[0],&this->recvStatuses_[0]);
          //if something has been received, accumulate and potentially forward it
          for(Int i = 0;i<recvCount;++i ){
            Int idx = this->recvDoneIdx_[i];

            if(idx!=MPI_UNDEFINED){
              Int size = 0;
              MPI_Get_count(&this->recvStatuses_[i], MPI_BYTE, &size);


#if ( _DEBUGlevel_ >= 1 ) || defined(REDUCE_VERBOSE)
              statusOFS<<myRank_<<" RECVD from "<<this->recvStatuses_[i].MPI_SOURCE<<" on tag "<<tag_<<std::endl;
#endif
              if(size>0){
                //resize if needed
                if(this->sendDataPtrs_.size()<1){
                  this->sendDataPtrs_.assign(1,NULL);
                }

                //If sendDataPtrs is 0, allocate to the size of what has been received
                if(this->sendDataPtrs_[0]==NULL){
                  this->sendTempBuffer_.resize(this->msgSize_);
                  this->sendDataPtrs_[0] = (T*)&this->sendTempBuffer_[0];
                  Int nelem = msgSize_;
                  std::fill(this->sendDataPtrs_[0],this->sendDataPtrs_[0]+nelem,ZERO<T>());
                }

                //This is where the handle would be called
                reduce(idx,i);

              }

              this->recvCount_++;
            }
          }

          if(this->recvCount_== this->GetDestCount()){
            retVal = true;
          }
          else{
            retVal = false;
          }
        }
      }
      return retVal;
    }

  template< typename T> 
    inline bool TreeReduce_v2<T>::Progress(){

      bool retVal = false;
      if(done_){
        retVal = true;
      }
      else{
        //Do we need this ?
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
          }
        }
      }

      if(retVal){
        this->done_ = retVal;
        //TODO do some smart cleanup here
      }
      return retVal;
    }

  template< typename T> 
    inline T * TreeReduce_v2<T>::GetLocalBuffer(){ 
    return this->sendDataPtrs_[0];
  }

  template< typename T> 
    inline void TreeReduce_v2<T>::SetLocalBuffer(T * locBuffer){
      if(this->sendDataPtrs_.size()<1){
        this->sendDataPtrs_.assign(1,NULL);
      }

      if(this->sendDataPtrs_[0]!=NULL && this->sendDataPtrs_[0]!=locBuffer){
        blas::Axpy(msgSize_, ONE<T>(), this->sendDataPtrs_[0], 1, locBuffer, 1 );
        this->sendTempBuffer_.clear(); 
      }

      this->sendDataPtrs_[0] = locBuffer;
    }

  template< typename T> 
    inline void TreeReduce_v2<T>::AllocRecvBuffers(){
      this->recvDataPtrs_.assign(GetDestCount(),NULL);
      this->recvTempBuffer_.resize(GetDestCount()*msgSize_);

      for( Int idxRecv = 0; idxRecv < GetDestCount(); ++idxRecv ){
        this->recvDataPtrs_[idxRecv] = (T*)&(this->recvTempBuffer_[idxRecv*msgSize_]);
      }

      this->recvRequests_.assign(GetDestCount(),MPI_REQUEST_NULL);
      this->recvStatuses_.resize(GetDestCount());
      this->recvDoneIdx_.resize(GetDestCount());

      this->sendRequests_.assign(1,MPI_REQUEST_NULL);

      this->isAllocated_ = true;
    }

  template< typename T> 
    inline void TreeReduce_v2<T>::Reset(){
      TreeBcast_v2<T>::Reset();
      this->isAllocated_=false;
    }

  template< typename T> 
    inline bool TreeReduce_v2<T>::isMessageForwarded(){
      bool retVal=false;

      if(!fwded_){
        //If data has been received but not forwarded 
        if(IsDataReceived()){
          forwardMessage();
        }
        retVal = false;
      }
      else{
        //If data has been forwared, check for completion of send requests
        int destCount = this->myRank_==this->myRoot_?0:1;
        int completed = 0;
        if(destCount>0){
          //test the send requests
          int flag = 0;

          this->sendDoneIdx_.resize(destCount);
          MPI_Testsome(destCount,sendRequests_.data(),&completed,sendDoneIdx_.data(),MPI_STATUSES_IGNORE);
        }
        this->sendCount_ += completed;
        retVal = this->sendCount_ == this->sendPostedCount_;

        //MPI_Testall(destCount,sendRequests_.data(),&flag,MPI_STATUSES_IGNORE);
        //retVal = flag==1;
      }
      return retVal;
    }


  template< typename T>
    inline TreeReduce_v2<T> * TreeReduce_v2<T>::Create(const MPI_Comm & pComm, Int * ranks, Int rank_cnt, Int msgSize, double rseed){
      //get communicator size
      Int nprocs = 0;
      MPI_Comm_size(pComm, &nprocs);

#if defined(FTREE)
      return new FTreeReduce_v2<T>(pComm,ranks,rank_cnt,msgSize);
#elif defined(MODBTREE)
      return new ModBTreeReduce_v2<T>(pComm,ranks,rank_cnt,msgSize, rseed);
#elif defined(BTREE)
      return new BTreeReduce<T>(pComm,ranks,rank_cnt,msgSize);
#elif defined(PALMTREE)
      return new PalmTreeReduce_v2<T>(pComm,ranks,rank_cnt,msgSize);
#endif


      if(nprocs<=FTREE_LIMIT){
#if ( _DEBUGlevel_ >= 1 ) || defined(REDUCE_VERBOSE)
        statusOFS<<"FLAT TREE USED"<<endl;
#endif
        return new FTreeReduce_v2<T>(pComm,ranks,rank_cnt,msgSize);
      }
      else{
#if ( _DEBUGlevel_ >= 1 ) || defined(REDUCE_VERBOSE)
        statusOFS<<"BINARY TREE USED"<<endl;
#endif
        return new ModBTreeReduce_v2<T>(pComm,ranks,rank_cnt,msgSize, rseed);
      }
    }



} //namespace PEXSI
#endif
