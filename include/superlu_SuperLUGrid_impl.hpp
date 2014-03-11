
namespace PEXSI{

inline SuperLUGrid<Real>::SuperLUGrid	( MPI_Comm comm, Int nprow, Int npcol )
{
#ifndef _RELEASE_
	PushCallStack("SuperLUGrid::SuperLUGrid");
#endif
	ptrData = new RealGridData;
	if( ptrData == NULL )
		throw std::runtime_error( "SuperLUGrid cannot be allocated." );

  ptrData->GridInit(comm, nprow, npcol);	

#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
} 		// -----  end of method SuperLUGrid::SuperLUGrid  ----- 

inline SuperLUGrid<Real>::~SuperLUGrid	(  )
{
#ifndef _RELEASE_
	PushCallStack("SuperLUGrid::~SuperLUGrid");
#endif
	// NOTE (07/21/2013): Since superlu_gridinit gets a copy of the
	// communicator, it is legal to call superlu_gridexit even if
	// grid->comm is a copy of MPI_COMM_WORLD.

  ptrData->GridExit();	
	
	delete ptrData;

#ifndef _RELEASE_
	PopCallStack();
#endif
	return ;
} 		// -----  end of method SuperLUGrid::~SuperLUGrid  ----- 

inline SuperLUGrid<Complex>::SuperLUGrid	( MPI_Comm comm, Int nprow, Int npcol )
{
#ifndef _RELEASE_
	PushCallStack("SuperLUGrid::SuperLUGrid");
#endif
	ptrData = new ComplexGridData;
	if( ptrData == NULL )
		throw std::runtime_error( "SuperLUGrid cannot be allocated." );

  ptrData->GridInit(comm, nprow, npcol);	

#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
} 		// -----  end of method SuperLUGrid::SuperLUGrid  ----- 


inline SuperLUGrid<Complex>::~SuperLUGrid	(  )
{
#ifndef _RELEASE_
	PushCallStack("SuperLUGrid::~SuperLUGrid");
#endif
	// NOTE (07/21/2013): Since superlu_gridinit gets a copy of the
	// communicator, it is legal to call superlu_gridexit even if
	// grid->comm is a copy of MPI_COMM_WORLD.

  ptrData->GridExit();	
	
	delete ptrData;

#ifndef _RELEASE_
	PopCallStack();
#endif
	return ;
} 		// -----  end of method SuperLUGrid::~SuperLUGrid  ----- 

}


