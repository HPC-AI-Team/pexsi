#ifndef _PEXSI_REDUCE_TREE_IMPL_V2_HPP_
#define _PEXSI_REDUCE_TREE_IMPL_V2_HPP_
namespace PEXSI{
  template<typename T>
    TreeReduce_v2<T>::TreeReduce_v2(const MPI_Comm & pComm, Int * ranks, Int rank_cnt, Int msgSize):TreeBcast_v2<T>(pComm,ranks,rank_cnt,msgSize){
      this->sendDataPtrs_.assign(1,NULL);
      this->sendRequests_.assign(1,MPI_REQUEST_NULL);
      this->isAllocated_=false;
    }




  template<typename T>
    TreeReduce_v2<T>::TreeReduce_v2(const TreeReduce_v2<T> & Tree){
      this->Copy(Tree);
    }

  template<typename T>
    TreeReduce_v2<T>::TreeReduce_v2():TreeBcast_v2<T>(){
    }

  template<typename T>
    TreeReduce_v2<T>::~TreeReduce_v2(){
      this->cleanupBuffers();
    }

  template<typename T>
    inline void TreeReduce_v2<T>::Copy(const TreeReduce_v2<T> & Tree){
      ((TreeBcast_v2<T>*)this)->Copy(*(const TreeBcast_v2<T>*)&Tree);

      this->sendDataPtrs_.assign(1,NULL);
      this->sendRequests_.assign(1,MPI_REQUEST_NULL);
      this->isAllocated_= Tree.isAllocated_;

      this->cleanupBuffers();
    }



  template<typename T>
    inline void TreeReduce_v2<T>::postRecv(){
      if(this->GetDestCount()>this->recvPostedCount_){
        for( Int idxRecv = 0; idxRecv < this->myDests_.size(); ++idxRecv ){
          Int iProc = this->myDests_[idxRecv];
          MPI_Irecv( (char*)this->recvDataPtrs_[idxRecv], this->msgSize_, this->type_, 
              iProc, this->tag_,this->comm_, &this->recvRequests_[idxRecv] );
          this->recvPostedCount_++;
        } // for (iProc)
      }
    }




  template<typename T>
    inline void TreeReduce_v2<T>::reduce( Int idxRecv, Int idReq){
      //add thing to my data
      blas::Axpy(this->msgSize_, ONE<T>(), this->recvDataPtrs_[idxRecv], 1, this->sendDataPtrs_[0], 1 );
    }

  template<typename T>
    inline void TreeReduce_v2<T>::forwardMessage(){ 
      if(this->isReady_){
        if(this->myRank_!=this->myRoot_){
          //forward to my root if I have reseived everything
          Int iProc = this->myRoot_;
          // Use Isend to send to multiple targets
          if(this->sendDataPtrs_.size()<1){
            this->sendDataPtrs_.assign(1,NULL);
          }

          int msgsz = this->sendDataPtrs_[0]==NULL?0:this->msgSize_;

          MPI_Isend((char*)this->sendDataPtrs_[0], msgsz, this->type_, 
              iProc, this->tag_,this->comm_, &this->sendRequests_[0] );
          this->sendPostedCount_++;
#ifdef COMM_PROFILE
          PROFILE_COMM(this->myGRank_,this->myGRoot_,this->tag_,msgsz);
#endif

#if ( _DEBUGlevel_ >= 1 ) || defined(REDUCE_VERBOSE)
          statusOFS<<this->myRank_<<" FWD to "<<iProc<<" on tag "<<this->tag_<<std::endl;
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
        else if(this->recvCount_<this->recvPostedCount_){
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
              statusOFS<<this->myRank_<<" RECVD from "<<this->recvStatuses_[i].MPI_SOURCE<<" on tag "<<this->tag_<<std::endl;
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
                  Int nelem = this->msgSize_;
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
        else if(this->recvPostedCount_<this->GetDestCount()){
          this->postRecv();
          retVal = false;
        }
      }
      return retVal;
    }

  template< typename T> 
    inline bool TreeReduce_v2<T>::Progress(){

      bool retVal = false;
      if(this->done_){
        retVal = true;
      }
      else{
        //Do we need this ?
        if(!this->isAllocated_){
          AllocRecvBuffers();
        }

        if(this->isAllocated_){
          if(this->myRank_==this->myRoot_ && this->isAllocated_){
            this->isReady_=true;
          }

          if(this->isReady_){
            if(this->IsDataReceived()){

              //free the unnecessary arrays
              this->recvTempBuffer_.clear();
              this->recvRequests_.clear();
              this->recvStatuses_.clear();
              this->recvDoneIdx_.clear();

              if(this->isMessageForwarded()){
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


      if(!this->IsRoot()){
        //if not root, we need to allocate a temp buffer anyway
        if(this->sendDataPtrs_[0]==NULL){
          this->sendTempBuffer_.resize(this->msgSize_);
          this->sendDataPtrs_[0] = (T*)&this->sendTempBuffer_[0];
          Int nelem = this->msgSize_;
          std::fill(this->sendDataPtrs_[0],this->sendDataPtrs_[0]+nelem,ZERO<T>());
        }
        blas::Axpy(this->msgSize_, ONE<T>(), locBuffer, 1, this->sendDataPtrs_[0], 1 );
      }
      else{

        if(this->sendDataPtrs_[0]!=NULL && this->sendDataPtrs_[0]!=locBuffer){
          blas::Axpy(this->msgSize_, ONE<T>(), this->sendDataPtrs_[0], 1, locBuffer, 1 );
          this->sendTempBuffer_.clear(); 
        }

        this->sendDataPtrs_[0] = locBuffer;
      }
    }

  template< typename T> 
    inline void TreeReduce_v2<T>::AllocRecvBuffers(){
      this->recvDataPtrs_.assign(this->GetDestCount(),NULL);
      this->recvTempBuffer_.resize(this->GetDestCount()*this->msgSize_);

      for( Int idxRecv = 0; idxRecv < this->GetDestCount(); ++idxRecv ){
        this->recvDataPtrs_[idxRecv] = (T*)&(this->recvTempBuffer_[idxRecv*this->msgSize_]);
      }

      this->recvRequests_.assign(this->GetDestCount(),MPI_REQUEST_NULL);
      this->recvStatuses_.resize(this->GetDestCount());
      this->recvDoneIdx_.resize(this->GetDestCount());

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

      if(!this->fwded_){
        //If data has been received but not forwarded 
        if(this->IsDataReceived()){
          this->forwardMessage();
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
          MPI_Testsome(destCount,this->sendRequests_.data(),&completed,this->sendDoneIdx_.data(),MPI_STATUSES_IGNORE);
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

template< typename T>
  FTreeReduce_v2<T>::FTreeReduce_v2(const MPI_Comm & pComm, Int * ranks, Int rank_cnt, Int msgSize):TreeReduce_v2<T>(pComm, ranks, rank_cnt, msgSize){
    buildTree(ranks,rank_cnt);
  }

template< typename T>
  inline FTreeReduce_v2<T> * FTreeReduce_v2<T>::clone() const{
    FTreeReduce_v2<T> * out = new FTreeReduce_v2<T>(*this);
    return out;
  }



template< typename T>
  inline void FTreeReduce_v2<T>::buildTree(Int * ranks, Int rank_cnt){

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

template< typename T>
  inline void FTreeReduce_v2<T>::postRecv()
  {
    if(this->isAllocated_ && this->GetDestCount()>this->recvPostedCount_){
      MPI_Irecv( (char*)this->recvDataPtrs_[0], this->msgSize_, this->type_, 
          MPI_ANY_SOURCE, this->tag_,this->comm_, &this->recvRequests_[0] );
      this->recvPostedCount_++;
    }
  }

template< typename T>
  inline void FTreeReduce_v2<T>::AllocRecvBuffers(){
    this->recvDataPtrs_.assign(1,NULL);
    this->recvTempBuffer_.resize(this->msgSize_);

    this->recvDataPtrs_[0] = (T*)&(this->recvTempBuffer_[0]);

    this->recvRequests_.assign(1,MPI_REQUEST_NULL);
    this->recvStatuses_.resize(1);
    this->recvDoneIdx_.resize(1);
    this->sendRequests_.assign(1,MPI_REQUEST_NULL);

    this->isAllocated_ = true;
  }


template< typename T>
  inline bool FTreeReduce_v2<T>::Progress(){

          bool retVal = false;
          if(this->done_){
            retVal = true;
          }
          else{

            if(!this->isAllocated_){
              this->AllocRecvBuffers();
            }

            if(this->isAllocated_){
              if(this->myRank_==this->myRoot_ && this->isAllocated_){
                this->isReady_=true;
              }

              if(this->isReady_){
                if(this->IsDataReceived()){


                  //free the unnecessary arrays
                  this->recvTempBuffer_.clear();
                  this->recvRequests_.clear();
                  this->recvStatuses_.clear();
                  this->recvDoneIdx_.clear();

                  if(this->isMessageForwarded()){
                    retVal = true;
                  }
                }
                //else if(this->recvPostedCount_<this->GetDestCount()){
                //  //TODO check this
                //  if(this->recvPostedCount_==this->recvCount_){
                //    this->postRecv();
                //  }
                //}
              }
            }
          }

          if(retVal){
            this->done_ = retVal;
            //TODO do some smart cleanup
          }
          return retVal;
  }


template< typename T>
  BTreeReduce_v2<T>::BTreeReduce_v2(const MPI_Comm & pComm, Int * ranks, Int rank_cnt, Int msgSize):TreeReduce_v2<T>(pComm, ranks, rank_cnt, msgSize){
    buildTree(ranks,rank_cnt);
  }

template< typename T>
  inline BTreeReduce_v2<T> * BTreeReduce_v2<T>::clone() const{
    BTreeReduce_v2<T> * out = new BTreeReduce_v2<T>(*this);
    return out;
  }

template< typename T>
  inline void BTreeReduce_v2<T>::buildTree(Int * ranks, Int rank_cnt){
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


template< typename T>
  ModBTreeReduce_v2<T>::ModBTreeReduce_v2(const MPI_Comm & pComm, Int * ranks, Int rank_cnt, Int msgSize, double rseed):TreeReduce_v2<T>(pComm, ranks, rank_cnt, msgSize){
    this->rseed_ = rseed;
    buildTree(ranks,rank_cnt);
  }

template< typename T>
  inline void ModBTreeReduce_v2<T>::Copy(const ModBTreeReduce_v2<T> & Tree){
    ((TreeReduce_v2<T>*)this)->Copy(*((const TreeReduce_v2<T>*)&Tree));
    this->rseed_ = Tree.rseed_;
  }

template< typename T>
  inline ModBTreeReduce_v2<T> * ModBTreeReduce_v2<T>::clone() const{
    ModBTreeReduce_v2<T> * out = new ModBTreeReduce_v2<T>(*this);
    return out;
  }

template< typename T>
inline void ModBTreeReduce_v2<T>::buildTree(Int * ranks, Int rank_cnt){

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

template< typename T>
  PalmTreeReduce_v2<T>::PalmTreeReduce_v2(const MPI_Comm & pComm, Int * ranks, Int rank_cnt, Int msgSize):TreeReduce_v2<T>(pComm,ranks,rank_cnt,msgSize){
    //build the binary tree;
    buildTree(ranks,rank_cnt);
  }

template< typename T>
  inline PalmTreeReduce_v2<T> * PalmTreeReduce_v2<T>::clone() const{
    PalmTreeReduce_v2<T> * out = new PalmTreeReduce_v2<T>(*this);
    return out;
  }


template< typename T>
  inline void PalmTreeReduce_v2<T>::buildTree(Int * ranks, Int rank_cnt){
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




} //namespace PEXSI
#endif
