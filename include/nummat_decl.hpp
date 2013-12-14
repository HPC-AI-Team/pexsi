/*
	 Copyright (c) 2012 The Regents of the University of California,
	 through Lawrence Berkeley National Laboratory.  

   Authors: Lexing Ying and Lin Lin
	 
   This file is part of PEXSI. All rights reserved.

	 Redistribution and use in source and binary forms, with or without
	 modification, are permitted provided that the following conditions are met:

	 (1) Redistributions of source code must retain the above copyright notice, this
	 list of conditions and the following disclaimer.
	 (2) Redistributions in binary form must reproduce the above copyright notice,
	 this list of conditions and the following disclaimer in the documentation
	 and/or other materials provided with the distribution.
	 (3) Neither the name of the University of California, Lawrence Berkeley
	 National Laboratory, U.S. Dept. of Energy nor the names of its contributors may
	 be used to endorse or promote products derived from this software without
	 specific prior written permission.

	 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
	 ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
	 WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
	 DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
	 ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
	 (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
	 LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
	 ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
	 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
	 SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

	 You are under no obligation whatsoever to provide any bug fixes, patches, or
	 upgrades to the features, functionality or performance of the source code
	 ("Enhancements") to anyone; however, if you choose to make your Enhancements
	 available either publicly, or directly to Lawrence Berkeley National
	 Laboratory, without imposing a separate written license agreement for such
	 Enhancements, then you hereby grant the following license: a non-exclusive,
	 royalty-free perpetual license to install, use, modify, prepare derivative
	 works, incorporate into other computer software, distribute, and sublicense
	 such enhancements or derivative works thereof, in binary and source code form.
*/
/// @file nummat_decl.hpp
/// @brief Numerical matrix.
/// @date 2010-09-27
#ifndef _NUMMAT_DECL_HPP_
#define _NUMMAT_DECL_HPP_

#include "environment.hpp"

// TODO Move the things from decl to impl

//#ifdef _NUMMAT_VECTOR_
//#define bool char
//#endif


namespace  PEXSI{


#ifdef _NUMMAT_VECTOR_
  //	template <typename F>
  //		struct NumMatType
  //		{
  //      using type = F;
  //    };
  //
  //
  //	template <>
  //		struct NumMatType<bool>
  //		{
  //      using type = char;
  //    };
  //
  //template<typename t, typename... p>
  //using fixed_vector = std::vector<typename NumMatType<F>::type, p...>;


  //template <class F>class NumMat;
  //typedef NumMat<bool> NumMat<char>;

#endif


  /// @class NumMat
  ///
  /// @brief Numerical matrix.
  ///
  /// NumMat is a portable encapsulation of a pointer to represent a 2D
  /// matrix, which can either own (owndata == true) or view (owndata ==
  /// false) a piece of data.  
  template <class F>
    class NumMat
    {
    public:


      Int bufsize_; 

      /// @brief Whether it owns the data.
      bool owndata_;

      /// @brief The size of the first dimension.
      Int m_; 

      /// @brief The size of second dimension.
      Int n_;

#ifdef _NUMMAT_VECTOR_
      std::vector<F> * container_;
      /// @brief The pointer for the actual data.
      F* data_;
#else


      /// @brief The pointer for the actual data.
      F* data_;



      inline void allocate(F* data=NULL) {
        if(owndata_) {
#ifdef _NUMMAT_VECTOR_
          container_ = new std::vector<F>(m_*n_);
          if(data!=NULL){std::copy(data,data+m_*n_,container_->begin());}
          data_=&(*container_)[0];
#else
          if(m_>0 && n_>0) { data_ = new F[m_*n_]; if( data_ == NULL ) throw std::runtime_error("Cannot allocate memory."); } else data_=NULL;
          if(data!=NULL){std::copy(data,data+m_*n_,data_);}
#endif
        } else {
          data_ = data;
        }
#ifndef _NUMMAT_VECTOR_
        bufsize_ = m_*n_;
#endif
      }
      inline void deallocate(){
        if(owndata_) {
#ifdef _NUMMAT_VECTOR_
          delete container_;
#else
          if(bufsize_>0) { delete[] data_; data_ = NULL; }
#endif
        }
      }

    public:
      NumMat(Int m=0, Int n=0): m_(m), n_(n), owndata_(true) {
        this->allocate();
      }

      NumMat(Int m, Int n, bool owndata, F* data): m_(m), n_(n), owndata_(owndata) {
        this->allocate(data);
      }

      NumMat(const NumMat& C): m_(C.m_), n_(C.n_), owndata_(C.owndata_) {
        this->allocate(C.data_);
      }
      ~NumMat() {
        this->deallocate();
      }

      NumMat& Copy(const NumMat& C) {
        this->deallocate();
        m_ = C.m_; n_=C.n_; owndata_=C.owndata_;
        this->allocate(C.data_);
        return *this;
      }

      NumMat& operator=(const NumMat& C) {
        this->deallocate();
        m_ = C.m_; n_=C.n_; owndata_=C.owndata_;
        this->allocate(C.data_);
        return *this;
      }


      void Resize(Int m, Int n)  {
        if( owndata_ == false ){
          throw std::logic_error("Matrix being resized must own data.");
        }


#ifdef _NUMMAT_VECTOR_
        if(container_->size()<m*n)
        {
          container_->resize(m*n);
          data_=&(*container_)[0];
        }
        m_ = m; n_ = n;
#else
        if(m*n > bufsize_) {
          this->deallocate();
          m_ = m; n_ = n;
          this->allocate();
        }
        else{
          m_ = m; n_ = n;
        }
#endif
      }
      const F& operator()(Int i, Int j) const  { 
        if( i < 0 || i >= m_ ||
            j < 0 || j >= n_ ) {
          throw std::logic_error( "Index is out of bound." );
        }
        return data_[i+j*m_];
      }

      F& operator()(Int i, Int j)  { 
        if( i < 0 || i >= m_ ||
            j < 0 || j >= n_ ) {
          throw std::logic_error( "Index is out of bound." );
        }
        return data_[i+j*m_];
      }

#ifdef _NUMMAT_VECTOR_
      std::vector<F> * Container() const { return container_; }
#endif
      F* Data() const { return data_; }

      F* VecData(Int j)  const 
      { 
        if( j < 0 || j >= n_ ) {
          throw std::logic_error( "Index is out of bound." );
        }
        return &(data_[j*m_]); 
      }

      Int m() const { return m_; }
      Int n() const { return n_; }

      Int Size() const {return m_*n_;}
      Int ByteSize() const { return m_*n_*sizeof(F);}
#ifdef _NUMMAT_VECTOR_
      Int AllocatedSize() const {return container_->capacity();}
#else
      Int AllocatedSize() const {return bufsize_;}
#endif

//#else

      /////
      /////      Int header_size_;
      /////
      /////			/// @brief The size of the first dimension.
      /////			Int m_; 
      /////
      /////			/// @brief The size of second dimension.
      /////			Int n_;
      /////
      /////
      /////			char * container_;
      /////
      /////			/// @brief The pointer for the actual data.
      /////			F* data_;
      /////
      /////      inline void allocate(Int supidx,F* data=NULL) {
      /////        if(owndata_) {
      /////          if(m_>0 && n_>0) { data_ = new char[3*sizeof(Int)+sizeof(F)*m_*n_ ]; if( data_ == NULL ) throw std::runtime_error("Cannot allocate memory."); } else data_=NULL;
      /////          if(data!=NULL){std::copy(data,data+m_*n_,data_);}
      /////        } else {
      /////          data_ = data;
      /////        }
      /////        bufsize_ = m_*n_;
      /////      }
      /////      inline void deallocate(){
      /////				if(owndata_) {
      /////					if(bufsize_>0) { delete[] data_; data_ = NULL; }
      /////				}
      /////      }
      /////
      /////
      /////
      /////
      /////		public:
      /////			NumMat(Int m=0, Int n=0): m_(m), n_(n), owndata_(true) {
      /////        this->allocate();
      /////			}
      /////
      /////			NumMat(Int m, Int n, bool owndata, F* data): m_(m), n_(n), owndata_(owndata) {
      /////        this->allocate(data);
      /////			}
      /////
      /////			NumMat(const NumMat& C): m_(C.m_), n_(C.n_), owndata_(C.owndata_) {
      /////        this->allocate(C.data_);
      /////			}
      /////			~NumMat() {
      /////        this->deallocate();
      /////			}
      /////
      /////			NumMat& Copy(const NumMat& C) {
      /////        this->deallocate();
      /////				m_ = C.m_; n_=C.n_; owndata_=C.owndata_;
      /////        this->allocate(C.data_);
      /////				return *this;
      /////			}
      /////
      /////			NumMat& operator=(const NumMat& C) {
      /////        this->deallocate();
      /////				m_ = C.m_; n_=C.n_; owndata_=C.owndata_;
      /////        this->allocate(C.data_);
      /////				return *this;
      /////			}
      /////
      /////
      /////			void Resize(Int m, Int n)  {
      /////				if( owndata_ == false ){
      /////					throw std::logic_error("Matrix being resized must own data.");
      /////				}
      /////
      /////        
      /////#ifdef _NUMMAT_VECTOR_
      /////        if(container_->size()<m*n)
      /////        {
      /////          container_->resize(m*n);
      /////          data_=&(*container_)[0];
      /////        }
      /////				m_ = m; n_ = n;
      /////#else
      /////				if(m*n > bufsize_) {
      /////          this->deallocate();
      /////				  m_ = m; n_ = n;
      /////					this->allocate();
      /////        }
      /////        else{
      /////				  m_ = m; n_ = n;
      /////        }
      /////#endif
      /////			}
      /////			const F& operator()(Int i, Int j) const  { 
      /////				if( i < 0 || i >= m_ ||
      /////						j < 0 || j >= n_ ) {
      /////					throw std::logic_error( "Index is out of bound." );
      /////				}
      /////				return data_[i+j*m_];
      /////			}
      /////
      /////			F& operator()(Int i, Int j)  { 
      /////				if( i < 0 || i >= m_ ||
      /////						j < 0 || j >= n_ ) {
      /////					throw std::logic_error( "Index is out of bound." );
      /////				}
      /////				return data_[i+j*m_];
      /////			}
      /////
      /////
      /////
      /////
      /////
      /////
      /////
      /////			F* Data() const { return data_; }
      /////
      /////			F* VecData(Int j)  const 
      /////			{ 
      /////				if( j < 0 || j >= n_ ) {
      /////					throw std::logic_error( "Index is out of bound." );
      /////				}
      /////				return &(data_[j*m_]); 
      /////			}
      /////
      /////			Int m() const { return m_; }
      /////			Int n() const { return n_; }
      /////
      /////      Int Size() const {return m_*n_;}
      /////      Int ByteSize() const { return m_*n_*sizeof(F);}
      /////      Int AllocatedSize() const {return bufsize_;}
      /////
#endif

    };

  // Commonly used
#ifdef _NUMMAT_VECTOR_
  typedef NumMat<char>     BolNumMat;
#else
  typedef NumMat<bool>     BolNumMat;
#endif
  typedef NumMat<Int>      IntNumMat;
  typedef NumMat<Real>     DblNumMat;
  typedef NumMat<Complex>  CpxNumMat;

  // *********************************************************************
  // Utility functions
  // *********************************************************************
  /// @brief SetValue sets a numerical matrix to a constant val.
  template <class F> inline void SetValue(NumMat<F>& M, F val);

  /// @brief Energy computes the L2 norm of a matrix (treated as a vector).
  template <class F> inline Real Energy(const NumMat<F>& M);

} // namespace PEXSI

#endif // _NUMMAT_DECL_HPP_
